import os
from logging import debug, info, error
from json import load
from subprocess import call
from abc import ABC, abstractmethod
from tespy.tools.helpers import modify_path_os
from shutil import copy
from pathlib import Path

import fileinput

def replace(filename, text_to_search, replacement_text):
    with fileinput.FileInput(filename, inplace=True, backup='.bak') as file:
        for line in file:
            print(line.replace(text_to_search, replacement_text), end='')


class GeoStorageSimulator(ABC):
    @abstractmethod
    def preprocess(self):
        pass
    @abstractmethod
    def run(self):
        pass
    @abstractmethod
    def postprocess(self):
        pass


class OgsKb1(GeoStorageSimulator):
    def __init__(self, data):
        self.__data = data

        self.__directory = os.path.dirname(data['simulation_files'])
        self.__basename = os.path.basename(data['simulation_files'])

        self.__distribution = data['distribution']
        self.__factor = float(data['factor'])
        self.__density = float(data['density'])
        # print(data)

    def preprocess(self, T_ff_sto, m_sto, storage_mode):
        """
        - write input files for geostorage simulation
        - Requirements:
            - _*.st prepared with keywords $WARM, $COLD for input / output source terms
            -  *.bc prepared with keywords $TYPE, $VALUE for boundary condition at inlet
        :param T_ff_sto: feed flow temperature to geostorage  (in *.bc )
        :param m_sto: (float) mass flow rate through heat exchanger at geostorage side (in *.st)
        :param storage_mode: (str) 'charging' or 'discharging'
        :return:
        """

        info('GEOSTORAGE inflow temperature: {}'.format(T_ff_sto))

        for i in range(len(self.__distribution)):
            flow_rate = max(1.e-6,  # required for eskilson model
                            self.__factor * float(self.__distribution[i]) *
                            m_sto / self.__density)
            info('GEOSTORAGE flow rate {}: {}'.format(i, flow_rate))

            # ST-FILE
            st_file = os.path.join(self.__directory, self.__basename + '.st')

            copy(os.path.join(self.__directory, '_' + self.__basename + '.st'), st_file)

            replace(st_file, "$INFLOW_TEMPERATURE", str(T_ff_sto))
            # ATES
            replace(st_file, "$WARM_{}".format(i), str(flow_rate))
            replace(st_file, "$COLD_{}".format(i), str(-flow_rate))
            # BTES
            replace(st_file, "$FLOW_RATE_{}".format(i), str(flow_rate))

            # BC-FILE
            # ATES
            try:
                bc_file = os.path.join(self.__directory, self.__basename + '.bc')
                copy(os.path.join(self.__directory, '_' + self.__basename + '.bc'), bc_file)

                replace(bc_file, "$INFLOW_TEMPERATURE", str(T_ff_sto + 273.15))
                if storage_mode == 'charging':
                    replace(bc_file, "$INFLOW_POSITION_{0}".format(i), "WARM_{0}".format(i))
                elif storage_mode == 'discharging':
                    replace(bc_file, "$INFLOW_POSITION_{0}".format(i), "COLD_{0}".format(i))
                else:
                    raise RuntimeError("Preprocess - Storage operation type unknown")
            except:
                pass

    def run(self):
        """
        - call geostorage simulator after file preparation
        :return:
        """
        Path(os.path.join(self.__directory, 'out.txt')).touch()
        Path(os.path.join(self.__directory, 'error.txt')).touch()
        try:
            call([self.__data['simulator_file'], self.__data['simulation_files']],
                stdout=open(os.path.join(self.__directory, 'out.txt')),
                stderr=open(os.path.join(self.__directory, 'error.txt')))
        except:
            error("Call of storage simulator failed")
        info('GEOSTORAGE calculation completed')

    def postprocess(self, storage_mode):
        """
        evaluate result
        :param storage_mode: (string) 'charging' or 'discharging'
        :return: return flow temperature from geostorage
        """

        for i in range(len(self.__distribution)):
            try:
                with open(os.path.join(self.__directory,
                                       self.__basename + '_HEAT_TRANSPORT_Contraflow_{}.tec'.format(i))) as file:
                    line = file.readline()
                    line = file.readlines()[1]
                    t_rf_sto = float(line.split()[2])
            except:
                position = 'PNT_COLD' if storage_mode == 'charging' else 'PNT_WARM'
                try:
                    with open(os.path.join(self.__directory, self.__basename + '_time_{}_{}.tec'.format(position, i))) as file:
                        line = file.readlines()[4]
                        t_rf_sto = float(line.split()[1]) - 273.15
                except:
                    t_fr_sto = None
                    raise RuntimeError("Postprocess - No output for storage {}".format(i))

        return t_rf_sto


class GeoStorage:
    def __init__(self, cd):
        info('GEOSTORAGE Reading input file .geostorage_ctr.json')
        base_path = os.path.join(cd.working_dir, cd.geostorage_path)
        path = os.path.join(base_path, cd.scenario + '.geostorage_ctrl.json')
        # print("PATH: {}".format(path))
        self.__specification = dict()
        with open(path) as file:
            self.__specification.update(load(file))

        self.__specification['simulator_file'] = modify_path_os(self.__specification['simulator_file'])
        self.__specification['simulation_files'] = modify_path_os(self.__specification['simulation_files'])
        if self.__specification['simulator_file'][0] == '.':
            self.__specification['simulator_file'] = base_path + self.__specification['simulator_file']
        if self.__specification['simulation_files'][0] == '.':
            self.__specification['simulation_files'] = base_path + self.__specification['simulation_files']

        if self.__specification['simulator_name'] == 'ogs_kb1':
            info('GEOSTORAGE simulator is OGS_kb1')
            self.__simulator = OgsKb1({'simulator_file': self.__specification['simulator_file'],
                                        'simulation_files': self.__specification['simulation_files'],
                                       'factor': self.__specification['factor'],
                                       'distribution': self.__specification['distribution'],
                                       'density': self.__specification['density']
                                       })
        else:
            error('GEOSTORAGE simulator not supported')
            self.__simulator = None
        # print(self.__specification)

        self.__flag_belowMinimumTemperature = False
        self.__T_min = float(self.__specification['minimum_discharge_temperature'])


    def simulation_files(self):
        return self.__specification['simulation_files']

    def output_points(self):
        return self.__specification['output_points']

    def storage_type(self):
        return self.__specification['storage_type']

    def run_storage_simulation(self, T_ff_sto, m_sto, storage_mode):
        """
        :param T_ff_sto: (float) feed flow temperature to geostorage
        :param m_sto: (float) mass flow rate through heat exchanger at geostorage side
        :param storage_mode: (str) 'charging' or 'discharging'
        :return: return flow temperature from geostorage
        """
        self.__simulator.preprocess(T_ff_sto, m_sto, storage_mode)
        self.__simulator.run()

        T_rf_sto = self.__simulator.postprocess(storage_mode)

        if T_rf_sto < self.__T_min and self.__flag_belowMinimumTemperature == False:
            debug("GEOSTORAGE Discharge below minimum temperature")
            self.__flag_belowMinimumTemperature = True

        if self.__flag_belowMinimumTemperature:
            T_diff = T_rf_sto - self.__T_min
            debug("T_diff: {} - T_rf_sto: {} - m_sto: {}".format(T_diff, T_rf_sto, m_sto))

            if T_diff < -0.01:
                self.run_storage_simulation(T_ff_sto, m_sto/(1 - 0.2 * T_diff), storage_mode)
            elif T_diff > 0.01:
                self.run_storage_simulation(T_ff_sto, m_sto*(1 + 0.2 * T_diff), storage_mode)
            else:
                info("GEOSTORAGE Inner iteration converged")

        return T_rf_sto, self.__flag_belowMinimumTemperature, m_sto
