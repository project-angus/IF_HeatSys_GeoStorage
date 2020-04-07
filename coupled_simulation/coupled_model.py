import numpy as np
import datetime
import os, shutil
from logging import info, error
from pandas import DataFrame, read_csv

import powerplant as pp
from geostorage import GeoStorage
from properties import Properties


class CoupledModel:
    def __init__(self, path):
        """
         - creates power plant and storage models
        - reads input timeseries
        :param path: (string) path to *.main_ctrl.json
        """
        self.__prop = Properties(path=path[0])
        self.__gs = GeoStorage(self.__prop)
        self.__pp_info = pp.load_models(self.__prop)

        self.__directory = os.path.dirname(self.__gs.simulation_files())
        self.__basename = os.path.basename(self.__gs.simulation_files())
        info('INTERFACE powerplant models loaded')
        self.__input_ts = read_csv(os.path.join(self.__prop.working_dir, self.__prop.input_timeseries_file),
                                      delimiter=',', decimal='.')

        info('INTERFACE Input time series read')
        cols = ['time', 'Q_target', 'Q_actual', 'Q_sto', 'P_plant', 'ti_plant',
                'T_ff_sys', 'T_rf_sys', 'T_ff_sto', 'T_rf_sto','m_sto', 'pp_err']
        self.__output_ts = DataFrame(index=np.arange(0, self.__prop.t_steps_total),
                                     columns=cols)

    def prepare_timestepping(self):
        """
        - copy _HEAT_TRANSPORT.IC,  _LIQUID_FLOW.IC such that
        - initializes return temperature from storage for power plant model (takes feed-in temperture of heat network)
        pressure and temperature field are updated from primary_variables-file in time step loop
        :return: initial return temperature from storage for power plants
        """
        try:
            shutil.copy(os.path.join(self.__directory, '_HEAT_TRANSPORT.IC'),
                        os.path.join(self.__directory,
                                     self.__basename + '_HEAT_TRANSPORT_domain_primary_variables.txt'))
            # for ATES
            shutil.copy(os.path.join(self.__directory, '_LIQUID_FLOW.IC'),
                        os.path.join(self.__directory,
                                     self.__basename + '_LIQUID_FLOW_domain_primary_variables.txt'))
        except:
            pass

        info('INTERFACE time stepping prepared')

        return 40, 30  # INITIALIZATION

    def execute(self):
        """
        - prepare and execute time step loop
        :return:
        """
        T_DC, T_C = self.prepare_timestepping()
        T_rf_sto_0 = {'charging': T_C, 'discharging': T_DC, 'shutin': 0}

        for t_step in range(self.__prop.t_steps_total):

            current_time = datetime.timedelta(seconds=t_step * self.__prop.t_step_length) + self.__prop.t_start
            info('---------------------------------------------------')
            info('INTERFACE time step {} - {}'.format(t_step, current_time))

            Q_target, T_ff_sys, T_rf_sys = self.get_timestep_data(str(current_time))

            info('Target heat flow: {}'.format(Q_target))

            #storage_mode = 'charging' if Q_target > 1.e-3 elif Q_target < -1.e-3 'discharging' else 'shutin'
            if Q_target > 1.e-3:
                storage_mode = 'charging'
            elif Q_target < -1.e-3:
                storage_mode = 'discharging'
            else:
                storage_mode = 'shutin'

            #storage_mode = 'charging' if Q_target > 1.e-3 elif Q_target < -1.e-3 'discharging' else 'shutin'
            Q_sto, Q_sys, P_plant, ti_plant, T_ff_sto, T_rf_sto, m_sto, pp_err = \
                self.execute_timestep(Q_target, T_ff_sys, T_rf_sys, T_rf_sto_0[storage_mode], storage_mode)

            T_rf_sto_0[storage_mode] = T_rf_sto

            # postprocess timestep
            self.evaluate_timestep(t_step, current_time, Q_target, Q_sys,
                                   Q_sto, P_plant, ti_plant, T_ff_sys,
                                   T_rf_sys, T_ff_sto, T_rf_sto, m_sto, pp_err)

            if t_step % self.__prop.save_nth_t_step == 0:
                self.__output_ts.to_csv(
                    os.path.join(self.__prop.working_dir, self.__prop.output_timeseries_path), index=False)

        self.__output_ts.to_csv(os.path.join(self.__prop.working_dir, self.__prop.output_timeseries_path), index=False)
        print(self.__output_ts)

    def get_timestep_data(self, current_time):
        """
        - read input_timeseries
        - update ICs for geostorage
        :param current_time: (str)
        :return: (floats) target heat, temperatures and pressure to / from heat network (input to execute_timestep())
        """
        try:
            Q_target = float(self.__input_ts.loc[self.__input_ts.time == current_time, 'heat_target'].values[0])
            T_ff_sys = float(self.__input_ts.loc[self.__input_ts.time == current_time, 'temperature_feed'].values[0])
            T_rf_sys = float(self.__input_ts.loc[self.__input_ts.time == current_time, 'temperature_return'].values[0])

            shutil.copy(os.path.join(self.__directory,
                                     self.__basename + '_HEAT_TRANSPORT_domain_primary_variables.txt'),
                        os.path.join(self.__directory, 'HEAT_TRANSPORT.IC'))
            try:
                shutil.copy(os.path.join(self.__directory,
                                     self.__basename + '_LIQUID_FLOW_domain_primary_variables.txt'),
                            os.path.join(self.__directory, 'LIQUID_FLOW.IC'))
            except:
                pass
        except KeyError:
            Q_target, T_ff_sys, T_rf_sys, p_ff_sys, p_rf_sys = None, None, None, None, None
            error('INTERFACE time step not found in input data')

        return Q_target, T_ff_sys, T_rf_sys

    def execute_timestep(self, Q_target, T_ff_sys, T_rf_sys, T_rf_sto, storage_mode):
        """
        - contains interation loop
        - imitialize T_rf_sto with T_ff_sys
        :param Q_target: (float) Heat to heat network
        :param T_ff_sys: (float) temperature to heat network
        :param T_rf_sys: (float) temperature from heat network
        :param T_rf_sto: (float) Initial return temperature from geostorage
        :return: calculation results (floats) and name_plant (string)
        """
        # initialization
        Q = Q_target
        gs_belowMinumumTemperature = False # for power plant model !!!!!!!!!!

        # iteration
        for iter in range(self.__prop.iter_max):
            info('---------------------------------------------------')
            info('INTERFACE start iteration {}'.format(iter))
            try:
                # powerplant
                if abs(Q) > 1e-3:
                    if gs_belowMinumumTemperature is True:
                        Q_sto, Q_sys, P_plant, ti_plant, T_ff_sto, T_rf_sto, m_sto, pp_err = pp.calc_interface_params_limitation(
                            self.__pp_info, T_ff_sys, T_rf_sys, T_rf_sto, m_sto, storage_mode)

                    else:
                        Q_sto, Q_sys, P_plant, ti_plant, T_ff_sto, T_rf_sto, m_sto, pp_err = pp.calc_interface_params(
                            self.__pp_info, T_ff_sys, T_rf_sys, T_rf_sto, abs(Q), storage_mode)

                    if gs_belowMinumumTemperature is True:
                        Q_sto, Q_sys, P_plant, ti_plant, T_ff_sto, T_rf_sto, m_sto, pp_err = pp.calc_interface_params_limitation(
                            self.__pp_info, T_ff_sys, T_rf_sys, T_rf_sto, m_sto, storage_mode)

                    else:
                        Q_sto, Q_sys, P_plant, ti_plant, T_ff_sto, T_rf_sto, m_sto, pp_err = pp.calc_interface_params(
                            self.__pp_info, T_ff_sys, T_rf_sys, T_rf_sto, abs(Q), storage_mode)
                else:
                    Q_sto, Q_sys, P_plant, ti_plant, T_ff_sto, T_rf_sto, m_sto, pp_err = \
                        0, 0, 0, 0, 0, T_rf_sto, 0, False

                #Q_sto, Q_sys, P_plant, ti_plant, T_ff_sto, T_rf_sto, m_sto, pp_err = pp.calc_interface_params(
                #    self.__pp_info, T_ff_sys, T_rf_sys, T_rf_sto, abs(Q), storage_mode) if abs(Q) > 1e-3 else 0, 0, 0, 0, 0, T_rf_sto, 0, False

                if m_sto == 0:
                    storage_mode = 'shutin'
                    T_rf_sto_geo = T_rf_sto
                else:
                    info('POWERPLANT calculation completed')
                    # geostorage
                    T_rf_sto_geo, gs_belowMinumumTemperature, m_sto_geo = self.__gs.run_storage_simulation(T_ff_sto, m_sto, storage_mode)
                    # evaluate
                    info('GEOSTORAGE return temperature: {}'.format(T_rf_sto_geo))

                else:
                    T_rf_sto_geo = T_rf_sto

                error_T = abs(T_rf_sto_geo - T_rf_sto)
                error_m = abs(m_sto_geo - m_sto)
                info('INTERFACE coupling error: {}'.format(np.linalg.norm([error_m, error_T])))
                if (error_T < self.__prop.temperature_return_error and
                        error_m < self.__prop.mass_flow_error and
                        iter >= self.__prop.iter_min-1 or abs(Q) < 1e-3):
                    info("INTERFACE loop converged")
                    break
                # update
                T_rf_sto = T_rf_sto_geo
                m_sto = m_sto_geo
            except:
                Q_sto, P_plant, Q_actual, ti_plant, T_ff_sto, T_rf_sto_geo, m_sto, pp_err = \
                    None, None, None, None, None, None, None, False
                error("INTERFACE iteration failed")

        return Q_sto, Q_sys, P_plant, ti_plant, T_ff_sto, T_rf_sto_geo, m_sto, pp_err

    def evaluate_timestep(self, t_step, current_time, Q_target, Q_sys, Q_sto, P_plant, ti_plant,
                          T_ff_sys, T_rf_sys, T_ff_sto, T_rf_sto, m_sto, pp_err):
        """
        - write output timeseries
        :param t_step: (int) current time step
        :param current_time: (datetime.timedelta)
        :param Q_target: (float)
        :param Q_sys: (float) heat from / to system
        :param Q_sto: (float) heat from / to storage
        :param P_plant: (float) power from active plant
        :param ti_plant: (float) thermal input from active plant
        :param T_ff_sys: (float) feed flow temperature to heat network
        :param T_rf_sys: (float) return temperture from heat network
        :param T_ff_sto: (float) feed flow temperature to geostorage
        :param T_rf_sto: (float) return flow temperature from geostorage
        :param m_sto: (float) mass flow rate through heat exchanger at geostorage side
        :param pp_err: (bool) indicates whether an error occurred in power plant simulation
        :return:
        """
        try:
            os.remove(os.path.join(self.__directory, self.__basename + '0000.vtk'))

            if t_step % self.__prop.save_nth_t_step == 0:
                info('Store vtk')
                shutil.copy(os.path.join(self.__directory, self.__basename + '0001.vtk'),
                            os.path.join(self.__directory, self.__basename + '000{}'.format(t_step)[-4:] + '.vtk'))
            else:
                os.remove(os.path.join(self.__directory, self.__basename + '0001.vtk'))

            if t_step % self.__prop.save_debug_nth_t_step == 0:
                info('BREAKPOINT - store geostorage results')

                shutil.copy(os.path.join(self.__directory,
                                         self.__basename + '_HEAT_TRANSPORT_domain_primary_variables.txt'),
                            os.path.join(self.__directory, 'HEAT_TRANSPORT_{}.IC'.format(t_step)))

                if self.__gs.storage_type() == 'ATES':
                    shutil.copy(os.path.join(
                        self.__directory, self.__basename + '_LIQUID_FLOW_domain_primary_variables.txt'),
                                os.path.join(self.__directory, 'LIQUID_FLOW_{}.IC'.format(t_step)))

            for pnt in self.__gs.output_points():
                info('{} {}'.format(pnt, t_step))

                file = open(os.path.join(self.__directory, self.__basename + '_time_{}.tec'.format(pnt)), 'r')
                for line in file:
                    w = line
                file.close()

                file = open(os.path.join(self.__directory, self.__basename + '_point_{}.txt'.format(pnt)),
                            'w' if t_step == 0 else 'a')
                file.write('{}\t{}\n'.format(t_step * self.__prop.t_step_length, w.split()[-1]))
                file.close()
        except:
            pass

        if t_step % self.__prop.save_debug_nth_t_step == 0:
            self.__pp_info['power_plant_models'][self.__pp_info['discharge']['name']]['plant'].instance.save(
                os.path.join(self.__prop.working_dir, 'pp_discharging_' + str(t_step)))
            self.__pp_info['power_plant_models'][self.__pp_info['charge']['name']]['plant'].instance.save(
                os.path.join(self.__prop.working_dir, 'pp_charging_' + str(t_step)))

        self.__output_ts.loc[t_step] = np.array([current_time, Q_target, Q_sys, Q_sto, P_plant, ti_plant,
                                                 T_ff_sys, T_rf_sys, T_ff_sto, T_rf_sto, m_sto, pp_err])
