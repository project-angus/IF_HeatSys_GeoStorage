#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 15:17:46 2018

__author__ = "witte, wtp"

"""

import sys
import getopt
import pandas as pd
import numpy as np
from coupled_simulation import powerplant as pp, geostorage as gs
import json
import datetime
import os

def __main__(path):
    """
    main function to initialise the calculation

    - creates power plant and storage models
    - reads input timeseries
    - starts the loop
    - writes results to .csv-file

    :param md: object containing the basic model data
    :type md: model_data object
    :returns: no return value
    """

#    if len(argv) == 0:
#        return
#
#    #read main input file and set control variables, e.g. paths, identifiers, ...
#    #path = (r'D:\Simulations\if_testcase\testcase.main_ctrl.json')
#    path = ''
#
#    try:
#        opts, args = getopt.getopt(argv,"hi:o:",["ipath="])
#    except getopt.GetoptError:
#        print('test.py -i <inputpath>')
#        sys.exit(2)
#    for opt, arg in opts:
#        if opt == '-h':
#            print('test.py -i <inputpath>')
#            sys.exit()
#        elif opt in ("-i", "--ipath"):
#            path = arg
#
#    if path[0] == "r":
#        path = path[1:]
#
#    path_log = path[:-15]
#    path_log += ".log"
#    #check if file exists and delete if necessary
#    if os.path.isfile(path_log):
#        os.remove(path_log)
#    #save screen output in file
#    sys.stdout = Logger(path_log)
#
#    print('Input file is:')
#    print(path)
#
#    print('######################################################################')
#    print('Assembling model data...')

    cd = coupling_data(path=path)
#
#    # create instances for power plant and storage
#    geostorage = gs.geo_sto(cd)

    power_plant_path = cd.working_dir + cd.powerplant_path + cd.scenario + '.powerplant_ctrl.json'
    power_plant_info = {}
    with open(power_plant_path) as f:
        power_plant_info.update(json.load(f))

    for key, val in power_plant_info['power_plant_models'].items():
        power_plant_info['power_plant_models'][key] = pp.model(cd, val)
        print(val, key)

    print('######################################################################')
    print('Reading input time series...')

    input_ts = read_series(cd.working_dir + cd.input_timeseries_path)

    #prepare data structures
    print('######################################################################')
    print('Preparing output data structures...')
    variable_list = []
    if cd.auto_eval_output == True:
        variable_list = ['time', 'power_target', 'massflow_target', 'power_actual', 'massflow_actual','storage_pressure', 'Tstep_accepted', 'delta_power', 'delta_massflow' ]
    else:
        variable_list = ['time', 'power_target', 'massflow_target', 'power_actual', 'massflow_actual', 'storage_pressure' ]
    #one output line per timestep...
    output_ts = pd.DataFrame(index=np.arange(0, cd.t_steps_total),columns=variable_list)

    #print(output_ts)
    '''debug values from here onwards'''
    #data = [0.0, 0.0]
    #data = geostorage.CallStorageSimulation(1.15741, 1, cd, 'charging')
    #data = geostorage.CallStorageSimulation(0.0, 2, cd, 'shut-in')
    #data = geostorage.CallStorageSimulation(-1.15741, 3, cd, 'discharging')
    '''end of debug values'''

    print('######################################################################')
    p0 = 0.0 #old pressure (from last time step / iter)
    # get initial pressure before the time loop
    p0, dummy_flow = geostorage.CallStorageSimulation(0.0, 0, 0, cd, 'init')
    print('Simulation initialzation completed.')
    print('######################################################################')

    last_time = cd.t_start
    for t_step in range(cd.t_steps_total):

        current_time = datetime.timedelta(seconds=t_step * cd.t_step_length) + cd.t_start

#        try:
#            power_target = input_ts.loc[current_time].power / 100 * 1e6
#            last_time = current_time
#        except KeyError:
#            power_target = input_ts.loc[last_time].power / 100 * 1e6
#
#        # calculate pressure, mass flow and power
#        p_actual, m_target, m_actual, power_actual, success,  = calc_timestep(
#                powerplant, geostorage, power_target, p0, cd, t_step)
#
#        # save last pressure (p1) for next time step as p0
#        p0 = p_actual
#        #deleting old files
#        geostorage.deleteFFile(t_step)
#
#        # write pressure, mass flow and power to .csv
#        if cd.auto_eval_output == True:
#            delta_power = abs(power_actual) - abs(power_target)
#            delta_massflow = abs(m_actual) - abs(m_target)
#
#            output_ts.loc[t_step] = np.array([current_time, power_target, m_target, power_actual, m_actual,
#                                                    p_actual, success, delta_power, delta_massflow])
#        else:
#            output_ts.loc[t_step] = np.array([current_time, power_target, m_target, power_actual, m_actual,
#                                                    p_actual])

        #Logger.flush()

        #sys.stdout.flush() #force flush of output

        #if t_step % cd.save_nth_t_step == 0:
        output_ts.to_csv(cd.working_dir + cd.output_timeseries_path, index=False)

#def calc_timestep(powerplant, geostorage, power, p0, md, tstep):
#    """
#    calculates one timestep of coupled power plant - storage simulation
#
#    :param powerplant: powerplant model
#    :type powerplant: powerplant.model object
#    :param storage: storage model
#    :type storage: storage.model object
#    :param power: scheduled power for timestep
#    :type power: float
#    :param p0: initual pressure at timestep
#    :type p0: float
#    :param md: object containing the basic model data
#    :type md: model_data object
#    :returns: - p1 (*float*) - interface pressure at the end of the timestep
#              - m_corr (*float*) - mass flow for this timestep
#              - power (*float*) - power plant's input/output power for this
#                timestep
#    """
#    tstep_accepted = False
#    storage_mode = ''
#
#    if power == 0:
#        m = 0
#        storage_mode = 'shut-in'
#    elif power < 0:
#        storage_mode = 'discharging'
#        m, power = powerplant.get_mass_flow(power, p0, storage_mode)
#    else:
#        storage_mode = 'charging'
#        m, power = powerplant.get_mass_flow(power, p0, storage_mode)
#
#    if m == 0:
#        storage_mode = "shut-in"
#
#
#    #moved inner iteration into timestep function,
#    #iterate until timestep is accepted
#    p0_temp = p0
#    print('######################################################################')
#    print('######################################################################')
#    print('######################################################################')
#    print('Advancing to timestep:\t', tstep, 'Operational mode is: ', storage_mode)
#
#    for iter_step in range(md.max_iter): #do time-specific iterations
#
#
#        if tstep_accepted:
#            print('Message: Timestep accepted after iteration ', iter_step - 1)
#            break
#        print('----------------------------------------------------------------------')
#        print('----------------------------------------------------------------------')
#        print('Current iteration:\t', iter_step)
#        print('----------------------------------------------------------------------')
#
#        #get pressure for the given target rate and the actually achieved flow rate from storage simulation
#        p1, m_corr = geostorage.CallStorageSimulation(m, tstep, iter_step, md, storage_mode )
#
#
#        if storage_mode == 'charging' or storage_mode == 'discharging':
#            # pressure check
#            if abs((p0_temp - p1) / p1) > md.pressure_diff_rel or abs(p0_temp - p1) > md.pressure_diff_abs:
#                m, power = powerplant.get_mass_flow(power, p1, storage_mode)
#                print('Adjusting mass flow rate.')
#                print('m / m_corr\t\t', '%.6f'%m, '/', '%.6f'%m_corr, '[kg/s]')
#                print('p0_new / p1\t\t', '%.6f'%p0_temp, '/', '%.6f'%p1, '[bars]')
#                if m == 0:
#                    print('Forcing shut-in mode as m is zero.')
#                    storage_mode = 'shut-in'
#
#            # if pressure check is successful, mass flow check:
#            # check for difference due to pressure limitations
#            elif abs(m_corr) <= 1E-5:
#                power = 0
#                tstep_accepted = True
#                print('Adjusting power to ZERO')
#                print('m / m_corr\t\t', '%.6f'%m, '/', '%.6f'%m_corr, '[kg/s]')
#                print('p0_new / p1\t\t', '%.6f'%p0_temp, '/', '%.6f'%p1, '[bars]')
#
#            elif abs((m - m_corr) / m_corr) > md.flow_diff_rel or abs(m - m_corr) > md.flow_diff_abs:
#                m, power = powerplant.get_power(m_corr, p1, storage_mode)
#                tstep_accepted = True
#                print('Adjusting power to ', power)
#                print('m / m_corr\t\t', '%.6f'%m, '/', '%.6f'%m_corr, '[kg/s]')
#                print('p0_new / p1\t\t', '%.6f'%p0_temp, '/', '%.6f'%p1, '[bars]')
#
#            else:
#                #return p1, m_corr, power
#                tstep_accepted = True
#                m = m_corr
#
#        elif storage_mode == "shut-in":
#            print('Force accepting timestep b/c storage shut-in')
#            tstep_accepted = True
#        else:
#            print('Problem: Storage mode not understood')
#            tstep_accepted = True
#        #saving old pressure
#        p0_temp = p1
#
#    if not tstep_accepted:
#        print('----------------------------------------------------------------------')
#        print('----------------------------------------------------------------------')
#        print('Problem: Results in timestep ', tstep, 'did not converge, accepting last iteration result.')
#
#    return p1, m, m_corr, power, tstep_accepted
#
#
#def read_series(path):
#    """
#    reads the input time series
#
#    :param path: path to input time series
#    :type path: str
#    :returns: ts (*pandas.DataFrame*) - dataframe containing the time series
#    """
#    ts = pd.read_csv(path, delimiter=',', decimal='.')
#    ts = ts.set_index(pd.to_datetime(ts.index, unit='h'))
#    ts['power'] = ts['input'] - ts['output']
#    return ts
#
#
#class coupling_data:
#    """
#    creates a data container with the main model parameters
#
#    :returns: no return value
#    """
#
#    def __init__(self, path):
#
#        # load data.json information into objects dictionary (= attributes of
#        # the object)
#        self.path = path
#
#        with open(path) as f:
#            self.__dict__.update(json.load(f))
#
#        self.auto_eval_output = False
#
#        if self.eval_output == "True":
#            self.auto_eval_output = True
#
#        self.coupled_simulation()
#
#    def coupled_simulation(self):
#        """
#        Function to set all required default data, e.g. well names, paths, ...
#
#        :returns: no return value
#        """
#
#        str_tmp = self.path.strip('.main_ctrl.json')
#        self.scenario = ""
#        self.working_dir = ""
#
#        i = 0
#        key = ""
#        if os.name == 'nt':
#            key = "\\"
#            self.powerplant_path = self.powerplant_path.replace('/', '\\')
#            self.geostorage_path = self.geostorage_path.replace('/', '\\')
#            self.input_timeseries_path = self.input_timeseries_path.replace('/', '\\')
#            self.output_timeseries_path = self.output_timeseries_path.replace('/', '\\')
#        elif os.name == 'posix':
#            key = "/"
#            self.powerplant_path = self.powerplant_path.replace('\\', '/')
#            self.geostorage_path = self.geostorage_path.replace('\\', '/')
#            self.input_timeseries_path = self.input_timeseries_path.replace('\\', '/')
#            self.output_timeseries_path = self.output_timeseries_path.replace('\\', '/')
#        else:
#            print('Error: OS not supported')
#
#
#        for c in str_tmp[::-1]:
#            if c == key:
#                self.working_dir = str_tmp[:-i]
#                break
#            self.scenario += c
#            i += 1
#
#        self.scenario = self.scenario[::-1]
#        self.debug = bool(self.debug)
#        date_format = '%Y-%m-%d %H:%M:%S'
#        self.t_start = datetime.datetime.strptime(self.t_start, date_format)
#
#        print('Reading inputile \"' + self.scenario + '.main_ctrl.json\" '
#              'in working directory \"' + self.working_dir + '\"')
#
#'''        if self.debug:
#            print('DEBUG-OUTPUT for main control data')
#            print('Time series path:\t' + self.input_timeseries_path)
#            print('Start time:\t' + str(self.t_start))
#            print('Time step length:\t' + str(self.t_step_length))
#            print('Number of time steps:\t' + str(self.t_steps_total))
#            print('Iteration limits:\t' + str(self.min_iter) + '\t' +
#                  str(self.max_iter))
#            print('Pressure convergence criteria:\t' +
#                  str(self.pressure_diff_abs) +
#                  ' bars\t' + str(self.pressure_diff_rel * 100) + ' %')
#            print('END of DEBUG-OUTPUT for main control data')'''
#
#class Logger(object):
#
#    def __init__(self, a_string):
#        self.terminal = sys.stdout
#        self.log = open(a_string, "a")
#
#    def write(self, message):
#        self.terminal.write(message)
#        self.log.write(message)
#
#    def flush(self):
#        #this flush method is needed for python 3 compatibility.
#        #this handles the flush command by doing nothing.
#        #you might want to specify some extra behavior here.
#        pass
#
##__main__(sys.argv[1:])
