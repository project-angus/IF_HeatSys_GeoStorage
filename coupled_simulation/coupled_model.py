import numpy as np
import datetime
import os
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

        info('INTERFACE powerplant models loaded')
        self.__input_ts = read_csv(self.__prop.working_dir + self.__prop.input_timeseries_file,
                                      delimiter=',', decimal='.')
        info('INTERFACE Input time series read')

        self.__output_ts = DataFrame(index=np.arange(0, self.__prop.t_steps_total),
                                        columns=['time', 'plant', 'Q_target', 'Q_plant', 'Q_sto',
                                                 'P_plant', 'ti_actual',
                                                 'T_ff_sys', 'T_rf_sys', 'T_ff_sto', 'T_rf_sto', 'm_sto'])

    def prepare_timestepping(self):
        """
        - copy _HEAT_TRANSPORT.IC,  _LIQUID_FLOW.IC such that
        pressure and temperature field are updated from primary_variables-file in time step loop
        :return:
        """
        os.system('cp {}/_HEAT_TRANSPORT.IC {}_HEAT_TRANSPORT_domain_primary_variables.txt'.format(
            os.path.dirname(self.__gs.simulation_files()), self.__gs.simulation_files()))
        os.system('cp {}/_LIQUID_FLOW.IC {}_LIQUID_FLOW_domain_primary_variables.txt'.format(
            os.path.dirname(self.__gs.simulation_files()), self.__gs.simulation_files()))
        info('INTERFACE time stepping prepared')

    def execute(self):
        """
        - prepare and execute time step loop
        :return:
        """
        self.prepare_timestepping()
        T_rf_sto = 80

        for t_step in range(self.__prop.t_steps_total):
	
            current_time = datetime.timedelta(seconds=t_step * self.__prop.t_step_length) + self.__prop.t_start
            info('---------------------------------------------------')
            info('INTERFACE time step {} - {}'.format(t_step, current_time))


            Q_target, T_ff_sys, T_rf_sys, p_ff_sys, p_rf_sys = self.get_timestep_data(str(current_time))

            Q_sto, name_plant, Q_plant, P_plant, ti_plant, T_ff_sto, T_rf_sto, m_sto = \
                self.execute_timestep(Q_target, T_ff_sys, T_rf_sys, p_ff_sys, p_rf_sys, T_rf_sto)
            # postprocess timestep
            self.evaluate_timestep(t_step, current_time, name_plant, Q_target,  Q_plant, Q_sto, P_plant, ti_plant,
                                                     T_ff_sys, T_rf_sys,T_ff_sto, T_rf_sto, m_sto)

        self.__output_ts.to_csv(self.__prop.working_dir + self.__prop.output_timeseries_path, index=False)
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
            p_ff_sys = float(self.__input_ts.loc[self.__input_ts.time == current_time, 'pressure_feed'].values[0])
            p_rf_sys = float(self.__input_ts.loc[self.__input_ts.time == current_time, 'pressure_return'].values[0])

            os.system('cp {}_HEAT_TRANSPORT_domain_primary_variables.txt {}/HEAT_TRANSPORT.IC'.format(
                self.__gs.simulation_files(), os.path.dirname(self.__gs.simulation_files())))
            os.system('cp {}_LIQUID_FLOW_domain_primary_variables.txt {}/LIQUID_FLOW.IC'.format(
                self.__gs.simulation_files(), os.path.dirname(self.__gs.simulation_files())))
        except KeyError:
            Q_target, T_ff_sys, T_rf_sys, p_ff_sys, p_rf_sys = None, None, None, None, None
            error('INTERFACE time step not found in input data')


        return Q_target, T_ff_sys, T_rf_sys, p_ff_sys, p_rf_sys

    def execute_timestep(self, Q_target, T_ff_sys, T_rf_sys, p_ff_sys, p_rf_sys, T_rf_sto_0):
        """
        - contains interation loop
        - imitialize T_rf_sto with T_ff_sys
        :param Q_target: (float) Heat to heat network
        :param T_ff_sys: (float) temperature to heat network
        :param T_rf_sys: (float) temperature from heat network
        :param p_ff_sys: (float) pressure to heat network
        :param p_rf_sys: (float) pressure from heat network
        :return: calculation resulsts (floats) and name_plant (string)
        """
        # initialization
        storage_mode = 'charging' if Q_target > 0. else 'discharging'  # TO_DO: case Q_target = 0.
        T_rf_sto = T_rf_sto_0
        Q = Q_target
        # iteration
        for iter in range(self.__prop.iter_max):
            info('---------------------------------------------------')
            info('INTERFACE start iteration {}'.format(iter))
            try:
                # powerplant
                Q_sto, name_plant, Q_plant, P_plant, ti_plant, T_ff_sto, T_rf_sto, m_sto = pp.calc_interface_params(
                    self.__pp_info, T_ff_sys, T_rf_sys, T_rf_sto, p_ff_sys, p_rf_sys, abs(Q), storage_mode, 0)

                info('POWERPLANT calculation completed')
                # geostorage
                T_rf_sto_geo = self.__gs.run_storage_simulation(T_ff_sto, m_sto, storage_mode)
                # evaluate
                info('T_rf_sto_geo: {}'.format(T_rf_sto_geo))

                error = abs(T_rf_sto_geo - T_rf_sto)
                info('INTERFACE coupling error: {}'.format(error))
                if error < self.__prop.temperature_return_error and iter >= self.__prop.iter_min-1:
                    info("INTERFACE loop converged")
                    break
                # update
                T_rf_sto = T_rf_sto_geo
            except:
                Q_sto, name_plant, Q_plant, P_plant, ti_plant, T_ff_sto, T_rf_sto_geo, m_sto = \
                    None, None, None, None, None, None, None, None
                error("INTERFACE iteration failed")

        return Q_sto, name_plant, Q_plant, P_plant, ti_plant, T_ff_sto, T_rf_sto_geo, m_sto

    def evaluate_timestep(self, t_step, current_time, name_plant, Q_target,  Q_plant, Q_sto, P_plant, ti_plant,
                          T_ff_sys, T_rf_sys,T_ff_sto, T_rf_sto, m_sto):
        """
        - write output timeseries
        :param t_step: (int) current time step
        :param current_time: (datetime.timedelta)
        :param name_plant: (str) name of active power plant
        :param Q_target: (float)
        :param Q_plant: (float) heat from active plant
        :param Q_sto: (float) heat from / to storage
        :param P_plant: (float) power from active plant
        :param ti_plant: (float) thermal input from active plant
        :param T_ff_sys: (float) feed flow temperature to heat network
        :param T_rf_sys: (float) return temperture from heat network
        :param T_ff_sto: (float) feed flow temperature to geostorage
        :param T_rf_sto: (float) return flow temperature from geostorage
        :param m_sto: (float) mass flow rate through heat exchanger at geostorage side
        :return:
        """
        self.__output_ts.loc[t_step] = np.array([current_time, name_plant, Q_target, Q_plant, Q_sto, P_plant, ti_plant,
                                                 T_ff_sys, T_rf_sys, T_ff_sto, T_rf_sto, m_sto])
