#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 15:17:46 2018

@author: witte
"""

# %% imports

import numpy as np
import logging
from tespy.networks import load_network
from tespy.tools import logger
from tespy.tools.helpers import TESPyNetworkError

import json

logger.define_logging(log_version=True, file_level=logging.DEBUG,
                      screen_level=logging.INFO)
# %% power plant model class


def load_models(cd):

    power_plant_path = (cd.working_dir + cd.powerplant_path +
                        cd.scenario + '.powerplant_ctrl.json')
    power_plant_models = {}
    with open(power_plant_path) as f:
        power_plant_models.update(json.load(f))

    for key, val in power_plant_models['power_plant_models'].items():
        power_plant_models['power_plant_models'][key]['plant'] = model(
                cd, val, key)

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
        msg = 'Successfully loaded TESPy model ' + name + '.'
        logging.debug(msg)

    def load_tespy_model(self):

        # load tespy models with the network_reader module
        self.instance = load_network(self.wdir + self.model_data['path'])
        self.instance.set_attr(
            m_range=self.model_data['m_range'], iterinfo=False)
        if self.model_data['debug'] is True:
            self.instance.set_attr(iterinfo=True)

        # design heat flow
        self.instance.imp_busses[self.model_data['heat_bus_sys']].set_attr(
            P=self.model_data['Q_design'])

        # design system temperatures
        self.instance.imp_conns[self.model_data['ff_sys']].set_attr(
            T=self.model_data['T_ff_sys_design'])

        self.instance.imp_conns[self.model_data['rf_sys']].set_attr(
            T=self.model_data['T_rf_sys_design'])

        # design storage temperatures
        self.instance.imp_conns[self.model_data['ff_sto']].set_attr(
            T=self.model_data['T_ff_sto_design'])

        self.instance.imp_conns[self.model_data['rf_sto']].set_attr(
            T=self.model_data['T_rf_sto_design'])

        # solve design case and save to path
        self.instance.solve('design')
        self.new_design = self.wdir + self.model_data['path'] + '_design'
        self.instance.save(self.new_design)

        # offdesign test
        self.instance.solve('offdesign', design_path=self.new_design,
                            init_path=self.new_design)

        # create low load data
        Q_range = np.linspace(
            self.model_data['Q_low'], 1, 3)[::-1] * self.model_data['Q_design']
        for Q in Q_range:
            self.instance.imp_busses[
                self.model_data['heat_bus_sys']].set_attr(P=Q)

            self.instance.solve('offdesign', design_path=self.new_design)
        self.instance.save(self.new_design + '_low_Q')


def calc_interface_params(ppinfo, T_ff_sys, T_rf_sys, T_rf_sto, Q, mode):
    """
    Calculates the interface parameters at the heat exchanger of the thermal
    energy storage. The heat exchanger is always part of the network, it may
    be altered, if requirements of the heating system cannot be met, depending
    on the calculation mode.

    **Charging**

    If the thermal energy storage is charged (at T_feed), the return flow
    temperature and the mass flow through the heat exchanger on the hot side
    are results of the simulation.

    **Discharging**

    If the thermal energy storage is discharged, the mass flow through the heat
    exchanger of the thermal energy storage is determined by the discharging
    heat flow and the temperature requirements. The actual feed flow
    temperature will be calculated by the heat exchanger providing the
    calculated mass flow. If the target temperature cannot be met, a power
    plant will be operated in order to meet the heat flow target. The power
    plant to operate is predefined in the configuration files:

    Different options are defined, one of which must be a flexible power plant,
    that can meet any part load. This option will be used, if all of the other
    options cannot deliver the required heat flow.

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

    Returns
    -------
    IF_DATA : list
        List of exchange parameters.
    """

    if mode == 'charging':

        # read interface model information
        charge = ppinfo['charge']['name']
        ttd_restriction = ppinfo['charge']['restricted']
        plant = ppinfo['power_plant_models'][charge]['plant']

        # check for minimum heat transfer
        Q_min = plant.model_data['Q_min'] * plant.model_data['Q_design']
        if Q < Q_min:
            # storage interface plant operation impossible
            msg = ('Target heat extraction rate of ' + str(round(Q, 0)) + ' W '
                   'below minimum possible heat extraction rate ofF ' +
                   str(round(Q_min, 0)) + ' W. '
                   'Interface operation impossible.')
            logging.warning(msg)
            T_ff_sto = T_rf_sto
            return 0, 0, 0, 0, T_ff_sto, T_rf_sto, 0, False

        else:

            # temperature restriction for storage feed flow temperature
            ttd = ppinfo['charge']['ttd_min']
            T_ff_sto_max = ppinfo['charge']['T_ff_sto_max']

            # system feed flow temperature as temperature at interface
            # inlet
            IF_data = sim_IF_charge(
                plant, T_ff_sys, T_rf_sto, Q, ttd, T_ff_sto_max)

        return IF_data

    elif mode == 'discharging':

        # read interface model information
        discharge = ppinfo['discharge']['name']
        plant = ppinfo['power_plant_models'][discharge]['plant']

        # check for minimum heat transfer
        Q_min = plant.model_data['Q_min'] * plant.model_data['Q_design']
        if Q < Q_min:
            # storage interface plant operation impossible
            msg = ('Target heat extraction rate of ' + str(round(Q, 0)) + ' W '
                   'below minimum possible heat extraction rate of ' +
                   str(round(Q_min, 0)) + ' W. '
                   'Interface operation impossible.')
            logging.warning(msg)
            T_ff_sto = T_rf_sto
            IF_data = 0, 0, 0, 0, T_ff_sto, T_rf_sto, 0, False

        else:
            IF_data = sim_IF_discharge(plant, T_ff_sys, T_rf_sys, T_rf_sto, Q)

        return IF_data

    else:
        msg = ('Mode for interface parameter calculation at '
               'coupled_simulation.powerplant must be \'charging\' or '
               '\'discharging\'.')
        logging.error(msg)
        raise ValueError(msg)


def calc_interface_params_limitation(ppinfo, T_ff_sys, T_rf_sys, T_rf_sto, m, mode):
    """

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

    m : float
        Mass flow.

    mode : str
        Charging or discharging mode.

    Returns
    -------
    IF_DATA : list
        List of exchange parameters.
    """

    if mode == 'charging':

        # read interface model information
        charge = ppinfo['charge']['name']
        ttd_restriction = ppinfo['charge']['restricted']
        plant = ppinfo['power_plant_models'][charge]['plant']

        # check for minimum heat transfer
        Q_min = plant.model_data['Q_min'] * plant.model_data['Q_design']
        if Q < Q_min:
            # storage interface plant operation impossible
            msg = ('Target heat extraction rate of ' + str(round(Q, 0)) + ' W '
                   'below minimum possible heat extraction rate ofF ' +
                   str(round(Q_min, 0)) + ' W. '
                   'Interface operation impossible.')
            logging.warning(msg)
            T_ff_sto = T_rf_sto
            return 0, 0, 0, 0, T_ff_sto, T_rf_sto, 0, False

        else:

            # temperature restriction for storage feed flow temperature
            ttd = ppinfo['charge']['ttd_min']
            T_ff_sto_max = ppinfo['charge']['T_ff_sto_max']

            # system feed flow temperature as temperature at interface
            # inlet
            IF_data = sim_IF_charge(
                plant, T_ff_sys, T_rf_sto, Q, ttd, T_ff_sto_max)

        return IF_data

    elif mode == 'discharging':

        # read interface model information
        discharge = ppinfo['discharge']['name']
        plant = ppinfo['power_plant_models'][discharge]['plant']
        IF_data = sim_IF_discharge_limitation(plant, T_ff_sys, T_rf_sys, T_rf_sto, m)
        return IF_data

    else:
        msg = ('Mode for interface parameter calculation at '
               'coupled_simulation.powerplant must be \'charging\' or '
               '\'discharging\'.')
        logging.error(msg)
        raise ValueError(msg)


def sim_IF_discharge(plant, T_ff_sys, T_rf_sys, T_rf_sto, Q):
    """
    Calculate the operation of the interface plant (e. g. heat exchanger or
    heat pump) without temperature restirctions.

    Parameters
    ----------
    plant : coupled_simulation.powerplant.model
        Power plant model object, holding power plant information.

    T_ff_sys : float
        Heating system feed flow temperature (to heating system).

    T_rf_sys : float
        Heating system return flow temperature (to interface).

    T_rf_sto : float
        Storage return flow temperature (from storage).

    Q : float
        Value of heat flow.

    Returns
    -------
    Q_sto : float
        Heat transferred from storage.

    Q_sys : float
        Heat transferred to system.

    P_IF : float
        Power input/output of storage interface.

    TI_IF : float
        Thermal input of storage interface.

    T_ff_sto : float
        Storage feed flow temperature (to storage).

    T_rf_sto : float
        Storage return flow temperature (from storage).

    mass_flow : float
        District heating water mass flow through the plant.

    err : bool
        Indicates whether an error occurred in calculation.
    """
    model = plant.instance

    # specify system parameters
    model.imp_busses[plant.model_data['heat_bus_sys']].set_attr(P=np.nan)
    model.imp_conns[plant.model_data['rf_sys']].set_attr(T=T_rf_sys)
    model.imp_conns[plant.model_data['rf_sto']].set_attr(T=T_rf_sto)
    model.imp_conns[plant.model_data['ff_sto']].set_attr(T=np.nan)
    model.imp_conns[plant.model_data['ff_sys']].set_attr(
            m=np.nan, T=T_ff_sys, design=[])
    model.imp_busses[plant.model_data['heat_bus_sys']].set_attr(P=Q)

    # solving
    try:
        design = plant.new_design
        if Q < plant.model_data['Q_design'] * plant.model_data['Q_low']:
            # initialise low heat transfer cases with init_path_low_Q
            try:
                model.solve('offdesign', design_path=design)
                if model.lin_dep or model.res[-1] > 1e-3:
                    raise TESPyNetworkError

            except (TESPyNetworkError, ValueError):
                init = plant.new_design + '_low_Q'
                model.solve('offdesign', design_path=design, init_path=init)
                if model.lin_dep or model.res[-1] > 1e-3:
                    raise TESPyNetworkError

        else:
            try:
                model.solve('offdesign', design_path=design)
                if model.lin_dep or model.res[-1] > 1e-3:
                    raise TESPyNetworkError

            except (TESPyNetworkError, ValueError):
                model.solve('offdesign', design_path=design, init_path=design)
                if model.lin_dep or model.res[-1] > 1e-3:
                    raise TESPyNetworkError

    except (TESPyNetworkError, ValueError):
        model.lin_dep = True

    if model.lin_dep or model.res[-1] > 1e-3:
        return 0, 0, 0, 0, 0, T_rf_sto, 0, True

    for conn_id, limits in plant.model_data['limiting_mass_flow'].items():
        conn = model.imp_conns[conn_id]
        m_max = conn.m.design * limits[1]
        m_min = conn.m.design * limits[0]
        m = conn.m.val_SI

        if m > m_max:
            model.imp_busses[plant.model_data['heat_bus_sys']].set_attr(
                P=np.nan)
            m_range = np.linspace(m_max, m, num=3, endpoint=False)
            for m_val in m_range[::-1]:
                conn.set_attr(m=m_val)
                model.solve('offdesign', design_path=design)
            msg = ('Limiting heat flow due to mass flow restriction in '
                   'extraction plant: mass flow: ' + str(round(m, 2)) +
                   'kg/s; maximum mass flow: ' + str(round(m_max, 2)) +
                    'kg/s.')
            logging.warning(msg)

        elif m < m_min:
            msg = ('Shutting off plant due to mass flow restriction in '
                   'extraction plant: mass flow: ' + str(round(m, 2)) +
                   'kg/s; minimum mass flow: ' + str(round(m_min, 2)) +
                    'kg/s.')
            logging.warning(msg)
            return 0, 0, 0, 0, 0, T_rf_sto, 0, False

        conn.set_attr(m=np.nan)

    # storage interface temperatures
    T_ff_sys = model.imp_conns[plant.model_data['ff_sys']].T.val
    T_ff_sto = model.imp_conns[plant.model_data['ff_sto']].T.val
    T_rf_sto = model.imp_conns[plant.model_data['rf_sto']].T.val

    # storage mass flow
    m_sto = model.imp_conns[plant.model_data['ff_sto']].m.val

    # interface transferred energy params
    Q_sto = model.imp_busses[plant.model_data['heat_bus_sto']].P.val
    Q_sys = model.imp_busses[plant.model_data['heat_bus_sys']].P.val
    P_IF = model.imp_busses[plant.model_data['power_bus']].P.val
    TI_IF = model.imp_busses[plant.model_data['ti_bus']].P.val

    return Q_sto, Q_sys, P_IF, TI_IF, T_ff_sto, T_rf_sto, m_sto, False


def sim_IF_charge(plant, T_rf_sys, T_rf_sto, Q, ttd, T_ff_sto_max):
    """
    Calculate the operation of the interface plant for storage charging.

    Parameters
    ----------
    plant : coupled_simulation.powerplant.model
        Power plant model object, holding power plant information.

    T_rf_sys : float
        Heating system return flow temperature (to interface).

    T_rf_sto : float
        Storage return flow temperature (from storage).

    Q : float
        Value of heat flow.

    ttd : float
        Minimum terminal temperature difference.

    T_ff_sto_max : float
        Maximum temperature to storage.

    Returns
    -------
    Q_sto : float
        Heat transferred to storage.

    Q_sys : float
        Heat transferred from system.

    P_IF : float
        Power input/output of storage interface.

    TI_IF : float
        Thermal input of storage interface.

    T_ff_sto : float
        Storage feed flow temperature (to storage).

    T_rf_sto : float
        Storage return flow temperature (from storage).

    mass_flow : float
        District heating water mass flow through the plant.

    err : bool
        Indicates whether an error occurred in calculation.
    """
    model = plant.instance
    model.imp_busses[plant.model_data['heat_bus_sys']].set_attr(P=np.nan)

    if T_ff_sto_max < T_rf_sys - ttd:
        T_ff_sto = T_ff_sto_max
    else:
        T_ff_sto = T_rf_sys - ttd

    # specify system parameters
    model.imp_busses[plant.model_data['heat_bus_sys']].set_attr(P=Q)
    model.imp_conns[plant.model_data['rf_sys']].set_attr(T=T_rf_sys)
    model.imp_conns[plant.model_data['rf_sto']].set_attr(T=T_rf_sto)
    model.imp_conns[plant.model_data['ff_sto']].set_attr(T=T_ff_sto)
    model.imp_conns[plant.model_data['ff_sys']].set_attr(
        m=np.nan, T=np.nan, design=[])

    # solving
    try:

        design = plant.new_design
        if Q < plant.model_data['Q_design'] * plant.model_data['Q_low']:
            # initialise low heat transfer cases with init_path_low_Q
            try:
                model.solve('offdesign', design_path=design)
                if model.lin_dep or model.res[-1] > 1e-3:
                    raise TESPyNetworkError

            except (TESPyNetworkError, ValueError):
                init = plant.new_design + '_low_Q'
                model.solve('offdesign', design_path=design, init_path=init)
        else:
            try:
                model.solve('offdesign', design_path=design)
                if model.lin_dep or model.res[-1] > 1e-3:
                    raise TESPyNetworkError

            except (TESPyNetworkError, ValueError):
                model.solve('offdesign', design_path=design, init_path=design)

        if model.lin_dep or model.res[-1] > 1e-3:
            raise TESPyNetworkError
    except (TESPyNetworkError, ValueError):
        model.lin_dep = True

    if model.lin_dep or model.res[-1] > 1e-3:
        return 0, 0, 0, 0, 0, T_rf_sto, 0, True

    for conn_id, limits in plant.model_data['limiting_mass_flow'].items():
        conn = model.imp_conns[conn_id]
        m_max = conn.m.design * limits[1]
        m_min = conn.m.design * limits[0]
        m = conn.m.val_SI

        if m > m_max:
            model.imp_busses[plant.model_data['heat_bus_sys']].set_attr(
                P=np.nan)
            m_range = np.linspace(m_max, m, num=3, endpoint=False)
            for m_val in m_range[::-1]:
                conn.set_attr(m=m_val)
                model.solve('offdesign', design_path=design)
            msg = ('Limiting heat flow due to mass flow restriction in '
                   'injection plant: mass flow: ' + str(round(m, 2)) +
                   'kg/s; maximum mass flow: ' + str(round(m_max, 2)) +
                    'kg/s.')
            logging.warning(msg)

        elif m < m_min:
            msg = ('Shutting off plant due to mass flow restriction in '
                   'injection plant: mass flow: ' + str(round(m, 2)) +
                   'kg/s; minimum mass flow: ' + str(round(m_min, 2)) +
                    'kg/s.')
            logging.warning(msg)
            return 0, 0, 0, 0, 0, T_rf_sto, 0, False

        conn.set_attr(m=np.nan)

    # storage interface temperatures
    T_ff_sys = model.imp_conns[plant.model_data['ff_sys']].T.val
    T_ff_sto = model.imp_conns[plant.model_data['ff_sto']].T.val
    T_rf_sto = model.imp_conns[plant.model_data['rf_sto']].T.val

    # storage mass flow
    m_sto = model.imp_conns[plant.model_data['ff_sto']].m.val

    # interface transferred energy params
    Q_sto = model.imp_busses[plant.model_data['heat_bus_sto']].P.val
    Q_sys = model.imp_busses[plant.model_data['heat_bus_sys']].P.val
    P_IF = model.imp_busses[plant.model_data['power_bus']].P.val
    TI_IF = model.imp_busses[plant.model_data['ti_bus']].P.val

    return Q_sto, Q_sys, P_IF, TI_IF, T_ff_sto, T_rf_sto, m_sto, False


def sim_IF_discharge(plant, T_ff_sys, T_rf_sys, T_rf_sto, m):
    """
    Calculate the operation of the interface plant (e. g. heat exchanger or
    heat pump) without temperature restirctions.

    Parameters
    ----------
    plant : coupled_simulation.powerplant.model
        Power plant model object, holding power plant information.

    T_ff_sys : float
        Heating system feed flow temperature (to heating system).

    T_rf_sys : float
        Heating system return flow temperature (to interface).

    T_rf_sto : float
        Storage return flow temperature (from storage).

    Q : float
        Value of heat flow.

    Returns
    -------
    Q_sto : float
        Heat transferred from storage.

    Q_sys : float
        Heat transferred to system.

    P_IF : float
        Power input/output of storage interface.

    TI_IF : float
        Thermal input of storage interface.

    T_ff_sto : float
        Storage feed flow temperature (to storage).

    T_rf_sto : float
        Storage return flow temperature (from storage).

    mass_flow : float
        District heating water mass flow through the plant.

    err : bool
        Indicates whether an error occurred in calculation.
    """
    model = plant.instance

    # specify system parameters
    model.imp_busses[plant.model_data['heat_bus_sys']].set_attr(P=np.nan)
    model.imp_conns[plant.model_data['rf_sys']].set_attr(T=T_rf_sys)
    model.imp_conns[plant.model_data['rf_sto']].set_attr(T=T_rf_sto, m=m)
    model.imp_conns[plant.model_data['ff_sto']].set_attr(T=np.nan)
    model.imp_conns[plant.model_data['ff_sys']].set_attr(
            m=np.nan, T=T_ff_sys, design=[])

    # solving
    try:
        try:
            model.solve('offdesign', design_path=design)
            if model.lin_dep or model.res[-1] > 1e-3:
                raise TESPyNetworkError

        except (TESPyNetworkError, ValueError):
            model.solve('offdesign', design_path=design, init_path=design)
            if model.lin_dep or model.res[-1] > 1e-3:
                raise TESPyNetworkError

    except (TESPyNetworkError, ValueError):
        model.lin_dep = True

    if model.lin_dep or model.res[-1] > 1e-3:
        return 0, 0, 0, 0, 0, T_rf_sto, 0, True

    for conn_id, limits in plant.model_data['limiting_mass_flow'].items():
        conn = model.imp_conns[conn_id]
        m_max = conn.m.design * limits[1]
        m_min = conn.m.design * limits[0]
        m = conn.m.val_SI

        if m > m_max:
            model.imp_conns[plant.model_data['rf_sto']].set_attr(m=np.nan)
            m_range = np.linspace(m_max, m, num=3, endpoint=False)
            for m_val in m_range[::-1]:
                conn.set_attr(m=m_val)
                model.solve('offdesign', design_path=design)
            msg = ('Limiting heat flow due to mass flow restriction in '
                   'extraction plant: mass flow: ' + str(round(m, 2)) +
                   'kg/s; maximum mass flow: ' + str(round(m_max, 2)) +
                    'kg/s.')
            logging.warning(msg)

        elif m < m_min:
            msg = ('Shutting off plant due to mass flow restriction in '
                   'extraction plant: mass flow: ' + str(round(m, 2)) +
                   'kg/s; minimum mass flow: ' + str(round(m_min, 2)) +
                    'kg/s.')
            logging.warning(msg)
            return 0, 0, 0, 0, 0, T_rf_sto, 0, False

        conn.set_attr(m=np.nan)

    # storage interface temperatures
    T_ff_sys = model.imp_conns[plant.model_data['ff_sys']].T.val
    T_ff_sto = model.imp_conns[plant.model_data['ff_sto']].T.val
    T_rf_sto = model.imp_conns[plant.model_data['rf_sto']].T.val

    # storage mass flow
    m_sto = model.imp_conns[plant.model_data['ff_sto']].m.val

    # interface transferred energy params
    Q_sto = model.imp_busses[plant.model_data['heat_bus_sto']].P.val
    Q_sys = model.imp_busses[plant.model_data['heat_bus_sys']].P.val
    P_IF = model.imp_busses[plant.model_data['power_bus']].P.val
    TI_IF = model.imp_busses[plant.model_data['ti_bus']].P.val

    return Q_sto, Q_sys, P_IF, TI_IF, T_ff_sto, T_rf_sto, m_sto, False
