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

logger.define_logging(
    log_path=True, log_version=True, screen_level=logging.INFO
)
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
        msg = 'Successfully loaded TESPy model ' + name + '.'
        logging.debug(msg)

    def load_tespy_model(self):

        # load tespy models with the network_reader module
        self.instance = nwkr.load_nwk(self.wdir + self.model_data['path'])
        self.instance.set_printoptions(print_level='none')

def calc_interface_params(models, order, T_ff_sys, T_rf_sys, T_ff_sto, T_rf_sto, Q, mode):
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

    2 Options are defined, one of which must be a flexible power plant, that can
    meet any part load. This option will be used, if the primary option cannot
    deliver the required heat flow.

    Parameters
    ----------
    models : dict
        Dictionary with power plant models of type coupled_simulation.powerplant.model
        identified by name, e. g.: he, ice_he, hp_he, eb_he, ....

    T_return : float
        Return flow temperature.

    T_feed : float
        Feed flow temperature.

    Q : float
        Total heat flow of the heat exchanger.

    mode : str
        Charging or discharging mode.
    """
    if mode == 'charging':

        'do something'

    elif mode == 'discharging':

        p = 10e5
        water.update(CP.PT_INPUTS, p, T_rf_sys + 273.15)
        h_return = water.hmass()
        water.update(CP.PT_INPUTS, p, T_ff_sys + 273.15)
        h_feed = water.hmass()
        mass_flow = Q / (h_feed - h_return)

        heat_ex = models['he_discharge']['model']
        he = heat_ex.instance
        he.imp_conns[heat_ex.model_data['ff_sys']].set_attr(m=mass_flow)
        he.imp_conns[heat_ex.model_data['rf_sys']].set_attr(T=T_rf_sys)
        he.imp_conns[heat_ex.model_data['ff_sto']].set_attr(T=T_ff_sto)
        he.imp_conns[heat_ex.model_data['rf_sto']].set_attr(T=T_rf_sto)
        he.imp_busses[heat_ex.model_data['heat_bus']].set_attr(P=np.nan)
        design = heat_ex.wdir + heat_ex.model_data['path']
        he.solve('offdesign', design_path=design, path_abs=True)
        Q_res = Q - he.imp_busses[heat_ex.model_data['heat_bus']].P.val
        T = he.imp_conns[heat_ex.model_data['ff_sys']].T.val
        p = he.imp_conns[heat_ex.model_data['ff_sys']].p.val

        if Q_res > 1e2:

            for model in order:
                if Q_res > models[model]['Q_design'] * models[model]['Q_low']:
                    plant = models[model]['model']
                    nw = plant.instance
                    nw.imp_conns[plant.model_data['ff_sys']].set_attr(T=T_ff_sys)
                    nw.imp_conns[plant.model_data['rf_sys']].set_attr(m=mass_flow, T=T, p=p)
                    nw.imp_busses[plant.model_data['heat_bus']].set_attr(P=np.nan)
                    design = plant.wdir + plant.model_data['path']
                    nw.set_printoptions(print_level='info')
                    nw.solve('offdesign', design_path=design, path_abs=True)
                    nw.print_results()

        else:

            he.imp_conns[heat_ex['ff_sys']].set_attr(m=np.nan)
            he.imp_busses[heat_ex['heat_bus']].set_attr(P=Q)
            he.solve('offdesign', design_path=design, path_abs=True)


#            print(models)
#            print(order)
#            print(model)
#        if Q_res > model

#        he.imp_conns[heat_ex['ff_sys']].T.val_SI


    else:
        msg = 'Mode for interface parameter calculation at coupled_simulation.powerplant must be \'charging\' or \'discharging\'.'
        logging.error(msg)
        ValueError(msg)


def calc_power_plant_operation(models, T_return, T_feed, Q):
    """

    """
