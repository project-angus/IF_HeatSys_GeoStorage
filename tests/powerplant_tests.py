#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# %% imports

from tespy import nwk, cmp, con, hlp, logger
from coupled_simulation import cp, pp
from nose.tools import eq_, raises
import numpy as np
import os

# %% component tests

cp.__main__(os.getcwd() + '/testdata/testcase.main_ctrl.json')
