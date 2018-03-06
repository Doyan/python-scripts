#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 13:36:53 2018

@author: gabgus
"""
import sys
sys.path.append('/chalmers/users/gabgus/miniconda2/envs/p3/lib/python3.5/site-packages/cantera/')


import numpy as np, matplotlib.pyplot as plt
import cantera as ct
#------------------------------------------------------------------------------
g = ct.Solution('gri30.xml')
g.TPX = 300.0, ct.one_atm, 'CH4:0.95,O2:2,N2:7.52'
g.equilibrate('TP')
