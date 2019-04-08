#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 15:17:46 2018

@author: witte
"""

# %% imports

import pandas as pd
import numpy as np
import logging
from tespy import nwkr, logger, con, hlp
from CoolProp import AbstractState as AS
import CoolProp as CP
import json

global water
water = AS('HEOS', 'water')

logger.define_logging(log_path=True, log_version=True, screen_level=logging.WARNING)
# %% power plant model class


def load_models(cd):

    power_plant_path = cd.working_dir + cd.powerplant_path + cd.scenario + '.powerplant_ctrl.json'
    power_plant_models = {}
    with open(power_plant_path) as f:
        power_plant_models.update(json.load(f))

    for key, val in power_plant_models['power_plant_models'].items():
        power_plant_models['power_plant_models'][key]['model'] = model(cd, val, key)

    return power_plant_models


class model:
    """
    Creates the model for the power plant. Parameters are loaded from
    coupling data object cd.

    Parameters
    ----------
    cd : coupling_data
        Generel data for the interface handling.
    """

    def __init__(self, cd, data, name):

        # load data.json information into objects dictionary (= attributes of
        # the object)
        self.wdir = cd.working_dir + cd.powerplant_path
        self.sc = cd.scenario
        self.model_data = data

        self.load_tespy_model()
        self.instance.solve('offdesign', design_path=self.wdir + self.model_data['path'],
                            init_path=self.wdir + self.model_data['path'], path_abs=True)
        msg = 'Successfully loaded TESPy model ' + name + '.'
        logging.debug(msg)

    def load_tespy_model(self):

        # load tespy models with the network_reader module
        self.instance = nwkr.load_nwk(self.wdir + self.model_data['path'])
        self.instance.set_printoptions(print_level='none')

def calc_interface_params(ppinfo, T_ff_sys, T_rf_sys, T_ff_sto, T_rf_sto, p_ff, p_rf, Q, mode, it):
    """
    Calculates the interface parameters at the heat exchanger of the thermal
    energy storage. The heat exchanger is always part of the network, it may
    be altered, if requirements of the heating system cannot be met, depending
    on the calculation mode.

    **Charging**

    If the thermal energy storage is charged (at T_feed), the return flow temperature
    and the mass flow through the heat exchanger on the hot side are results of the simulation.

    **Discharging**

    If the thermal energy storage is discharged, the mass flow through the heat exchanger of
    the thermal energy storage is determined by the discharging heat flow and the temperature requirements.
    The actual feed flow temperature will be calculated by the heat exchanger providing the
    calculated mass flow. If the target temperature cannot be met, a power plant
    will be operated in order to meet the heat flow target. The power plant to
    operate is predefined in the configuration files:

    Different options are defined, one of which must be a flexible power plant, that can
    meet any part load. This option will be used, if all of the other options cannot
    deliver the required heat flow.

    Parameters
    ----------
    ppinfo : dict
        Dictionary with information on TESPy power plant models.

    T_ff_sys : float
        Heating system feed flow temperature (to heating system).

    T_rf_sys : float
        Heating system return flow temperature (from heating system).

    T_ff_sto : float
        Storage feed flow temperature (to storage).

    T_rf_sto : float
        Storage return flow temperature (from storage).

    Q : float
        Total heat flow of the heat exchanger.

    mode : str
        Charging or discharging mode.

    p_ff : float
        Heating system feed flow pressure.

    p_rf : float
        Heating system return flow pressure.

    it : int
        Internal iteration count.

    Returns
    -------
    Q_sto : float
        Storage heat exchanger transferred heat.

    plant : str
        Name of additional power plant.

    Q_plant : float
        Transferred heat from additional power plant.

    P_plant : float
        Power input/output of additional power plant.

    ti_plant : float
        Thermal input of additional power plant.

    T_ff_sto :
        Storage feed flow temperature (to storage).

    T_rf_sto :
        Storage return flow temperature (from storage).

    m_sto : float
        Mass flow at storage side of heat exchanger.
    """

    if mode == 'charging':

        'do something'

    elif mode == 'discharging':

        ttd_min = 3

        # check if the storage return flow temperature is above system feed flow temperature
        if T_rf_sto > T_ff_sys + ttd_min:
            msg = ('Storage return flow temperature is above heating system feed flow temperature: '
                   'Setting system feed flow temperature and transferred heat at heat exchanger. '
                   'Storage feed flow temperature is result of calculation!')
            logging.debug(msg)
            # calculate mass flow for given return and feed flow temperature and transferred heat
            # this is for convergence stability purposes only!
            water.update(CP.PT_INPUTS, p_rf * 1e5, T_rf_sys + 273.15)
            h_return = water.hmass()
            water.update(CP.PT_INPUTS, p_ff * 1e5, T_ff_sys + 273.15)
            h_feed = water.hmass()
            mass_flow = Q / (h_feed - h_return)

            # load heat exchanger model and specify heat exchanger parameters
            heat_ex = ppinfo['power_plant_models']['he_discharge']['model']
            he = heat_ex.instance
#            he.iterinfo = True
            # Heating system feed flow temperature is a result of the calculation, storage temperatures are set.
            he.imp_conns[heat_ex.model_data['ff_sys']].set_attr(m=mass_flow, T=T_ff_sys, design=[])
            he.imp_busses[heat_ex.model_data['heat_bus']].set_attr(P=np.nan)
            he.imp_conns[heat_ex.model_data['rf_sys']].set_attr(T=T_rf_sys)
            he.imp_conns[heat_ex.model_data['ff_sto']].set_attr(T=np.nan)
            he.imp_conns[heat_ex.model_data['rf_sto']].set_attr(T=T_rf_sto)

            # solving
            design = heat_ex.wdir + heat_ex.model_data['path']
            if Q < heat_ex.model_data['Q_design'] * heat_ex.model_data['Q_low']:
                # initialise low heat transfer cases with init_path_low_Q
                init = heat_ex.wdir + heat_ex.model_data['init_path_low_Q']
                he.solve('offdesign', design_path=design, init_path=init, path_abs=True)
            else:
                he.solve('offdesign', design_path=design, init_path=design, path_abs=True)

            # solve again, this time with specified heat transfer
            he.imp_conns[heat_ex.model_data['ff_sys']].set_attr(m=np.nan, T=T_ff_sys)
            he.imp_busses[heat_ex.model_data['heat_bus']].set_attr(P=Q)
            # use values from previous calculation for initialisation
            he.solve('offdesign', design_path=design, path_abs=True)

            Q_sto = he.imp_busses[heat_ex.model_data['heat_bus']].P.val
            T_sys = he.imp_conns[heat_ex.model_data['ff_sys']].T.val
            T_ff_sto = he.imp_conns[heat_ex.model_data['ff_sto']].T.val
            T_rf_sto = he.imp_conns[heat_ex.model_data['rf_sto']].T.val
            m_sto = he.imp_conns[heat_ex.model_data['ff_sto']].m.val

            msg = ('Heat extraction from storage successful, feed flow temperature level met by storage heat exchanger only: '
                   'Transferred heat in ' + mode + ' mode is ' + str(round(he.imp_busses[heat_ex.model_data['heat_bus']].P.val, 0)) + ' W, target value was ' + str(round(Q, 0)) + ' W. '
                   'System feed flow temperature is at ' + str(round(T_sys, 1)) + ' C with storage feed flow temperature at ' + str(round(T_ff_sto, 1)) + ' C.')
            logging.info(msg)
            return Q_sto, None, 0, 0, 0, T_ff_sto, T_rf_sto, m_sto

        elif T_rf_sto <= T_rf_sys + ttd_min:

            if ppinfo['heat_balance_only']:
                # storage operation impossible, mass flow through heat exchanger on storage side is zero
                # temperature at feed and return flow are identical, "T_rf_sto, T_rf_sto" is not a typo!
                return 0, None, 0, 0, 0, T_rf_sto, T_rf_sto, 0

            else:
                plant, P, Q, ti, p_ff = sim_power_plant_operation(ppinfo, T_ff_sys, T_sys, p_ff, mass_flow, Q)
                if plant is None:
                   return 0, None, 0, 0, 0, T_rf_sto, T_rf_sto, 0

        else:
            # calculate mass flow for given return and feed flow temperature and transferred heat
            water.update(CP.PT_INPUTS, p_rf * 1e5, T_rf_sys + 273.15)
            h_return = water.hmass()
            water.update(CP.PT_INPUTS, p_ff * 1e5, T_ff_sys + 273.15)
            h_feed = water.hmass()
            mass_flow = Q / (h_feed - h_return)

            # load heat exchanger model and specify heat exchanger parameters
            heat_ex = ppinfo['power_plant_models']['he_discharge']['model']
            he = heat_ex.instance

            # Heating system feed flow temperature is a result of the calculation, storage temperatures are set.
            he.imp_conns[heat_ex.model_data['ff_sys']].set_attr(m=mass_flow, T=T_rf_sto - ttd_min)
            he.imp_busses[heat_ex.model_data['heat_bus']].set_attr(P=np.nan)
            he.imp_conns[heat_ex.model_data['rf_sys']].set_attr(T=T_rf_sys)
            he.imp_conns[heat_ex.model_data['ff_sto']].set_attr(T=np.nan)
            he.imp_conns[heat_ex.model_data['rf_sto']].set_attr(T=T_rf_sto)

            # solve system
            design = heat_ex.wdir + heat_ex.model_data['path']
            if Q < heat_ex.model_data['Q_design'] * heat_ex.model_data['Q_low']:
                init = heat_ex.wdir + heat_ex.model_data['init_path_low_Q']
                he.solve('offdesign', design_path=design, init_path=init, path_abs=True)
            else:
               he.solve('offdesign', design_path=design, init_path=design, path_abs=True)

            Q_sto = he.imp_busses[heat_ex.model_data['heat_bus']].P.val
            T_sys = he.imp_conns[heat_ex.model_data['ff_sys']].T.val
            T_ff_sto = he.imp_conns[heat_ex.model_data['ff_sto']].T.val
            T_rf_sto = he.imp_conns[heat_ex.model_data['rf_sto']].T.val
            m_sto = he.imp_conns[heat_ex.model_data['ff_sto']].m.val

            msg = ('Heat extraction from storage successful, feed flow temperature level NOT met by storage heat exchanger only: '
                   'Transferred heat in ' + mode + ' mode is ' + str(round(he.imp_busses[heat_ex.model_data['heat_bus']].P.val, 0)) + ' W, target value was ' + str(round(Q, 0)) + ' W. '
                   'System feed flow temperature is at ' + str(round(T_sys, 1)) + ' C with target temperature at ' + str(round(T_ff_sys, 1)) + ' C.')
            logging.info(msg)

            if ppinfo['heat_balance_only']:
                # return storage heat exchanger parameters, no operation calculation for power plant
                return Q_sto, None, 0, 0, 0, T_ff_sto, T_rf_sto, m_sto

            else:
            # calculate residual for transferred heat
                Q_res = Q - Q_sto

                p_sys = he.imp_conns[heat_ex.model_data['ff_sys']].p.val
                plant, P, Q_plant, ti, p_ff = sim_power_plant_operation(ppinfo, T_ff_sys, T_sys, p_sys, mass_flow, Q_res)

                if plant is None:
                   return Q_sto, None, 0, 0, 0, T_ff_sto, T_rf_sto, m_sto

                Q_diff = Q_res - Q_plant
                if abs(Q_diff) < ppinfo['Q_res_abs'] and abs(Q_diff / Q) < ppinfo['Q_res_rel']:
                    # calculation complete
                    msg = 'Calculation successful, residual value of heat flow is ' + str(round(Q_diff, 0)) + ' W.'
                    logging.info(msg)
                    return Q_sto, plant, P, Q_plant, ti, T_ff_sto, T_rf_sto, m_sto
                else:
                    if it < ppinfo['max_iter']:
                        it += 1
                        msg = 'Start iteration loop, as residual value of heat flow (' + str(round(Q_diff, 0)) + ' W) is higher than allowed (absolute deviation ' + str(ppinfo['Q_res_abs']) + ' W, relative deviation ' + str(ppinfo['Q_res_rel']) + ').'
                        logging.info(msg)
                        return calc_interface_params(ppinfo, T_ff_sys, T_rf_sys, T_ff_sto, T_rf_sto, p_ff, he.imp_conns[heat_ex.model_data['rf_sys']].p.val, Q, mode, it)
                    else:
                        msg = 'Reached maximum recursion depth, aborting calculation with residual heat flow of ' + str(round(Q_diff, 0)) + ' W.'
                        return Q_sto, plant, P, Q_plant, ti, T_ff_sto, T_rf_sto, m_sto

    else:
        msg = 'Mode for interface parameter calculation at coupled_simulation.powerplant must be \'charging\' or \'discharging\'.'
        logging.error(msg)
        raise ValueError(msg)


def sim_power_plant_operation(ppinfo, T_ff_sys, T_rf_sys, p_rf, mass_flow, Q_res):
    """
    Calculate

    **Charging**

    If the thermal energy storage is charged (at T_feed), the return flow temperature
    and the mass flow through the heat exchanger on the hot side are results of the simulation.

    **Discharging**

    If the thermal energy storage is discharged, the mass flow through the heat exchanger of
    the thermal energy storage is determined by the discharging heat flow and the temperature requirements.
    The actual feed flow temperature will be calculated by the heat exchanger providing the
    calculated mass flow. If the target temperature cannot be met, a power plant
    will be operated in order to meet the heat flow target. The power plant to
    operate is predefined in the configuration files:

    Different options are defined, one of which must be a flexible power plant, that can
    meet any part load. This option will be used, if all of the other options cannot
    deliver the required heat flow.

    Parameters
    ----------
    ppinfo : dict
        Dictionary with information on TESPy power plant models.

    T_ff_sys : float
        Heating system feed flow temperature (to heating system).

    T_rf_sys : float
        Heating system return flow temperature (from heat exchanger!!).

    p_rf : float
        Heating system return flow pressure.

    m : float
        District heating water mass flow through the plant.

    Returns
    -------
    plant : str
        Name of additional power plant.

    Q_plant : float
        Transferred heat from additional power plant.

    P_plant : float
        Power input/output of additional power plant.

    ti_plant : float
        Thermal input of additional power plant.

    p_ff : float
        Feed flow pressure of disctrict heating water.
    """

    models = ppinfo['power_plant_models']
    order = ppinfo['storage_boost_order']

    # go through the order of plants defined in the powerplant control file
    for model in order:
        # check if minimum heat flow of the plant is larger than required residual heat flow
        # jump to next power plant model if this is not the case
        if Q_res > models[model]['Q_design'] * models[model]['Q_min']:
            msg = 'Heat flow required (' + str(round(Q_res, 0)) + ' W) will be delivered by power plant \'' + model + '\'.'
            logging.info(msg)
            msg = 'power plant model: \'' + model + '\'.'
            logging.debug(msg)
            plant = models[model]['model']
            nw = plant.instance
            nw.imp_conns[plant.model_data['ff_sys']].set_attr(T=T_ff_sys)
            nw.imp_conns[plant.model_data['rf_sys']].set_attr(m=mass_flow, p=p_rf, T=T_rf_sys)
            nw.imp_busses[plant.model_data['heat_bus']].set_attr(P=np.nan)
            design = plant.wdir + plant.model_data['path']

            # solving
            design = plant.wdir + plant.model_data['path']
            if Q_res < plant.model_data['Q_design'] * plant.model_data['Q_low']:
                # initialise low heat transfer cases with init_path_low_Q
                init = plant.wdir + plant.model_data['init_path_low_Q']
                try:
                    nw.solve('offdesign', design_path=design, init_path=init, path_abs=True)
                except ValueError:
                    msg = 'Encountered ValueError in (most likely fluid property) calculation at plant \'' + model + '\'. Skipping to next plant in order. If this happens frequently, the power plant might not be suited for the respective task.'
                    logging.error(msg)
                    continue
            else:
                try:
                    nw.solve('offdesign', design_path=design, init_path=design, path_abs=True)
                except ValueError:
                    msg = 'Encountered ValueError in (most likely fluid property) calculation at plant \'' + model + '\'. Skipping to next plant in order. If this happens frequently, the power plant might not be suited for the respective task.'
                    logging.error(msg)
                    continue

            if nw.res[-1] < 1e-3 and not nw.lin_dep:
                P = nw.imp_busses[plant.model_data['power_bus']].P.val
                Q = nw.imp_busses[plant.model_data['heat_bus']].P.val
                ti = nw.imp_busses[plant.model_data['ti_bus']].P.val
                p_ff = nw.imp_conns[plant.model_data['ff_sys']].p.val
                msg = 'Calculation for power plant model \'' + model + '\' successful, heat flow is ' + str(round(Q, 0)) + ' W.'
                logging.info(msg)
                return model, P, Q, ti, p_ff
            else:
                msg = 'Error during simulation of power plant model \'' + model + '\', trying next plant.'
                logging.error(msg)

        else:
            msg = 'Heat flow required (' + str(round(Q_res, 0)) + ' W) below minimum possible heat flow for power plant \'' + model + '\' (' + str(round(models[model]['Q_design'] * models[model]['Q_min'])) + ' W).'
            logging.info(msg)

    msg = 'Calculation not successful, could not find a solution for any power plant model regarding the heat flow requirement.'
    logging.error(msg)
    return None, 0, 0, 0, 0
