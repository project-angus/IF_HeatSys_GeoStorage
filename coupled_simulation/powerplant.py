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

logger.define_logging(
    log_path=True, log_version=True, screen_level=logging.WARNING
)
# %% power plant model class


class model:
    """
    Creates the model for the power plant. Parameters are loaded from
    coupling data object cd.

    Parameters
    ----------
    cd : coupling_data
        Generel data for the interface handling.
    """

    def __init__(self, cd, data):

        # load data.json information into objects dictionary (= attributes of
        # the object)
        self.wdir = cd.working_dir + cd.powerplant_path
        self.sc = cd.scenario
        self.model_data = data

        self.load_tespy_model()

        print('test')

    def load_tespy_model(self):

        # load tespy models with the network_reader module
        self.instance = nwkr.load_nwk(self.wdir + self.model_data['path'])
        self.instance.set_printoptions(print_level='none')
