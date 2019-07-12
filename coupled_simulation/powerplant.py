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

logger.define_logging(log_path=True, log_version=True, screen_level=logging.INFO)
# %% power plant model class


def load_models(cd):

    power_plant_path = cd.working_dir + cd.powerplant_path + cd.scenario + '.powerplant_ctrl.json'
    power_plant_models = {}
    with open(power_plant_path) as f:
        power_plant_models.update(json.load(f))

    for key, val in power_plant_models['power_plant_models'].items():
        power_plant_models['power_plant_models'][key]['plant'] = model(cd, val, key)

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
        self.name = name

        self.load_tespy_model()
        self.instance.solve('offdesign', design_path=self.wdir + self.model_data['path'],
                            init_path=self.wdir + self.model_data['path'])
        msg = 'Successfully loaded TESPy model ' + name + '.'
        logging.debug(msg)

    def load_tespy_model(self):

        # load tespy models with the network_reader module
        self.instance = nwkr.load_nwk(self.wdir + self.model_data['path'])
        self.instance.set_printoptions(print_level='none')

def calc_interface_params(ppinfo, T_ff_sys, T_rf_sys, T_rf_sto, p_ff, p_rf, Q, mode, it):
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

        # calculate mass flow for given return and feed flow temperature and
        # transferred heat, this is for convergence stability purposes only!
        water.update(CP.PT_INPUTS, p_rf * 1e5, T_rf_sys + 273.15)
        h_return = water.hmass()
        water.update(CP.PT_INPUTS, p_ff * 1e5, T_ff_sys + 273.15)
        h_feed = water.hmass()
        mass_flow = Q / (h_feed - h_return)

        # specify storage interface properties according to temperature levels

        if T_rf_sto <= T_rf_sys + ttd_min:
            # heat extraction from storage impossible
            msg = 'Storage heat extraction impossible as storage temperature is below heating system return flow temperature.'
            logging.error(msg)
            Q_res = Q
            Q_sto = 0
            m_sto = 0
            T_ff_sto = T_rf_sto
            T_sys = T_rf_sys
            p_sys = p_rf

        else:
            # load heat exchanger model and specify heat exchanger parameters
            plant = ppinfo['power_plant_models']['he_discharge']['plant']
            model = plant.instance

            # specify parameters
            model.imp_busses[plant.model_data['heat_bus']].set_attr(P=np.nan)
            model.imp_conns[plant.model_data['rf_sys']].set_attr(T=T_rf_sys)
            model.imp_conns[plant.model_data['ff_sto']].set_attr(T=np.nan)
            model.imp_conns[plant.model_data['rf_sto']].set_attr(T=T_rf_sto)

            if T_rf_sto > T_ff_sys + ttd_min:
                # storage temperature above system feed flow temperature
                T = T_ff_sys

            else:
                # storage temperature below system feed flow temperature
                T = T_rf_sto - ttd_min

            # specify temperature value for system feed flow
            model.imp_conns[plant.model_data['ff_sys']].set_attr(m=mass_flow, T=T, design=[])

            # solving
            design = plant.wdir + plant.model_data['path']
            if Q < plant.model_data['Q_design'] * plant.model_data['Q_low']:
                # initialise low heat transfer cases with init_path_low_Q
                init = plant.wdir + plant.model_data['init_path_low_Q']
                model.solve('offdesign', design_path=design, init_path=init)
            else:
                model.solve('offdesign', design_path=design, init_path=design)

            if T_rf_sto > T_ff_sys + ttd_min:
                # system feed flow temperature requirement met by heat exchanger
                # solve again, this time with specified heat transfer
                model.imp_conns[plant.model_data['ff_sys']].set_attr(m=np.nan, T=T_ff_sys)
                model.imp_busses[plant.model_data['heat_bus']].set_attr(P=Q)
                # use values from previous calculation for initialisation
                model.solve('offdesign', design_path=design)

                # get storage heat exchanger parameters
                Q_sto = model.imp_busses[plant.model_data['heat_bus']].P.val
                T_sys = model.imp_conns[plant.model_data['ff_sys']].T.val
                T_ff_sto = model.imp_conns[plant.model_data['ff_sto']].T.val
                T_rf_sto = model.imp_conns[plant.model_data['rf_sto']].T.val
                m_sto = model.imp_conns[plant.model_data['ff_sto']].m.val

                msg = ('Heat extraction from storage successful, feed flow temperature level met by storage heat exchanger only: '
                       'Transferred heat in ' + mode + ' mode is ' + str(round(model.imp_busses[plant.model_data['heat_bus']].P.val, 0)) + ' W, target value was ' + str(round(Q, 0)) + ' W. '
                       'System feed flow temperature is at ' + str(round(T_sys, 1)) + ' C with storage feed flow temperature at ' + str(round(T_ff_sto, 1)) + ' C.')
                logging.info(msg)
                return Q_sto, None, 0, 0, 0, T_ff_sto, T_rf_sto, m_sto

            # get storage heat exchanger parameters
            Q_sto = model.imp_busses[plant.model_data['heat_bus']].P.val
            T_sys = model.imp_conns[plant.model_data['ff_sys']].T.val
            T_ff_sto = model.imp_conns[plant.model_data['ff_sto']].T.val
            T_rf_sto = model.imp_conns[plant.model_data['rf_sto']].T.val
            m_sto = model.imp_conns[plant.model_data['ff_sto']].m.val
            Q_res = Q - Q_sto
            # heating system pressure for storage booster is pressure after
            # heat exchanger
            p_sys = model.imp_conns[plant.model_data['ff_sys']].p.val

            # system feed flow temperature not met by heat exchanger
            # calculate storage booster operation
            msg = ('Heat extraction from storage successful, feed flow temperature level NOT met by storage heat exchanger only: '
                   'Transferred heat in ' + mode + ' mode is ' + str(round(model.imp_busses[plant.model_data['heat_bus']].P.val, 0)) + ' W, target value was ' + str(round(Q, 0)) + ' W. '
                   'System feed flow temperature is at ' + str(round(T_sys, 1)) + ' C with target temperature at ' + str(round(T_ff_sys, 1)) + ' C.')
            logging.info(msg)

        # storage booster calculation (both, for failed and partial heat extraction)
        if ppinfo['heat_balance_only']:
            # storage operation impossible, mass flow through heat exchanger on storage side is zero
            # temperature at feed and return flow are identical, "T_rf_sto, T_rf_sto" is not a typo!
            return 0, None, 0, 0, 0, T_ff_sto, T_rf_sto, m_sto

        else:
            # calculate power plant operation to replace heat extraction
            for sto_booster in ppinfo['storage_boost_order']:
                plant = ppinfo['power_plant_models'][sto_booster]['plant']

                data = sim_power_plant_operation(plant, T_ff_sys, T_sys, p_sys, mass_flow, Q_res)

                if data[0]:
                    # power plant operation possible
                    P_plant = data[1]
                    Q_plant = data[2]
                    ti_plant = data[3]
                    Q_diff = Q_res - Q_plant


                    if (abs(Q_diff / Q) > ppinfo['heat_max_dev_rel'] and
                            it < ppinfo['max_iter']):
                        msg = ('Repeat powerplant calculation, residual value of heat flow is ' + str(round(Q_diff, 0)) +
                               ' W, relative deviation of actual heat flow to target is too high: ' + str(round(Q_diff / Q, 4)) + '.')
                        logging.info(msg)
                        it += 1
                        p_ff = plant.instance.imp_conns[plant.model_data['ff_sys']].p.val
                        return calc_interface_params(ppinfo, T_ff_sys, T_rf_sys, T_rf_sto, p_ff, p_rf, Q, mode, it)

                    else:
                        msg = ('Calculation successful, residual value of heat flow is ' + str(round(Q_diff, 0)) + ' W, relative deviation is ' + str(round(Q_diff / Q, 4)) + '. '
                               'Power plant booster is \'' + sto_booster + '\'.')
                        logging.info(msg)

                        return Q_sto, sto_booster, P_plant, Q_plant, ti_plant, T_ff_sto, T_rf_sto, m_sto

            # no power plant operation possible
            msg = 'No power plant operation possible, heating system feed flow temperature could not be reached.'
            logging.error(msg)
            return Q_sto, None, 0, 0, 0, T_ff_sto, T_rf_sto, m_sto

    else:
        msg = 'Mode for interface parameter calculation at coupled_simulation.powerplant must be \'charging\' or \'discharging\'.'
        logging.error(msg)
        raise ValueError(msg)


def sim_power_plant_operation(plant, T_ff_sys, T_rf_sys, p_rf, mass_flow, Q_res):
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
    plant : coupled_simulation.powerplant.model
        Power plant model object, holding power plant information.

    T_ff_sys : float
        Heating system feed flow temperature (to heating system).

    T_rf_sys : float
        Heating system return flow temperature (from heat exchanger!!).

    p_rf : float
        Heating system return flow pressure.

    mass_flow : float
        District heating water mass flow through the plant.

    Q_res : float
        Value of residual heat flow.

    Returns
    -------
    success : bool
        Calculation successful.

    Q_plant : float
        Transferred heat from additional power plant.

    P_plant : float
        Power input/output of additional power plant.

    ti_plant : float
        Thermal input of additional power plant.
    """

    model = plant.instance
    Q_min = plant.model_data['Q_design'] * plant.model_data['Q_min']

    if Q_res < Q_min:
        msg = ('Power plant operation of plant \'' + plant.name + '\' not possible due to load restriction: '
               'Heat flow required (' + str(round(Q_res, 0)) + ' W) below minimum possible heat flow for power plant \'' + plant.name + '\' (' + str(round(Q_min, 0)) + ' W).')
        logging.debug(msg)
        return False, 0, 0, 0

    msg = 'Heat flow required (' + str(round(Q_res, 0)) + ' W) will be delivered by power plant \'' + plant.name + '\'.'
    logging.debug(msg)

    # parameter specification
    model.imp_conns[plant.model_data['ff_sys']].set_attr(T=T_ff_sys)
    model.imp_conns[plant.model_data['rf_sys']].set_attr(m=mass_flow, p=p_rf, T=T_rf_sys)
    model.imp_busses[plant.model_data['heat_bus']].set_attr(P=np.nan)

    # solving
    design = plant.wdir + plant.model_data['path']
    if Q_res < plant.model_data['Q_design'] * plant.model_data['Q_low']:
        # initialise low heat transfer cases with init_path_low_Q
        init = plant.wdir + plant.model_data['init_path_low_Q']
        try:
            model.solve('offdesign', design_path=design, init_path=init)
        except ValueError:
            msg = 'Encountered ValueError in (most likely fluid property) calculation at plant \'' + plant.name + '\'. Skipping to next plant in order. If this happens frequently, the power plant might not be suited for the respective task.'
            logging.error(msg)
            return False, 0, 0, 0
    else:
        try:
            model.solve('offdesign', design_path=design, init_path=design)
        except ValueError:
            msg = 'Encountered ValueError in (most likely fluid property) calculation at plant \'' + plant.name + '\'. Skipping to next plant in order. If this happens frequently, the power plant might not be suited for the respective task.'
            logging.error(msg)
            return False, 0, 0, 0

    if model.res[-1] < 1e-3 and model.lin_dep is False:
        P = model.imp_busses[plant.model_data['power_bus']].P.val
        Q = model.imp_busses[plant.model_data['heat_bus']].P.val
        ti = model.imp_busses[plant.model_data['ti_bus']].P.val
        msg = 'Calculation for power plant model \'' + plant.name + '\' successful, heat flow is ' + str(round(Q, 0)) + ' W.'
        logging.debug(msg)
        return True, P, Q, ti
    else:
        msg = 'Could not find a steady state solution for power plant model \'' + plant + '\' operation, trying next plant.'
        logging.warning(msg)
        return False, 0, 0, 0
