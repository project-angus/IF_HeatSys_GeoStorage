#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 15:17:46 2018

__author__ = "witte, wtp, jod"

"""
import sys
from coupled_model import CoupledModel


def run(path):
    """
    main function to initialise the calculation

    """
    coupled_model = CoupledModel(path)
    coupled_model.execute()


if __name__ == '__main__':
    run(sys.argv[1:])
