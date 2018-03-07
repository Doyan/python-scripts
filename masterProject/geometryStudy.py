#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 13:36:53 2018

@author: gabgus
"""
import sys
sys.path.append('/chalmers/users/gabgus/miniconda2/envs/p3/lib/python3.5/site-packages/cantera/')


import numpy as np
#import matplotlib.pyplot as plt
#import cantera as ct

#--------Key variables----------------------------------------------------------------------
alpha_w = 0.1 # Moisture content
SB = 1        # Steam-Biomass ratio
Y_ch = 0.1    # Char-yield


# Temperatures, in deg C
Ts = 400    # -- check sensitivity
Tbed = 850
Tfeed = 80  # total guess -- check sensitivity

# Composition of biomass C,H,O,Ash,H2O
Ydaf = np.array([54.55, 5.95, 39.5, 0, 0]) / 100 # from data for pine bark (N,S etc lumped in O)
Ywet = ((1-alpha_w)*Ydaf + (alpha_w)*np.array([0,0,0,0,1]))

#--- Table values ---------------------------------------------------------------------------
MW=np.array([12.0107, 1.00794, 15.999, -1000, 18.01528]) # g/mol -- Molar mass, Ash omitted

#--- Relations ------------------------------------------------------------------------------
def HHV_bm(Y):  # Changdong Sheng 2005 for dry biomass with lumped N and O
    return -1.3675 + 0.3137*Y[0] + 0.7009*Y[1] + 0.0318*Y[2]

def HHV_ch(Y): # IGT 1978  used for char (made for coal) 
    return 0.341*Y[0] + 1.323*Y[1] + 0.0685 - 0.0153*Y[3] - 0.1194*Y[2]

# Function for calculating enthalpy of combustion products.
def h_post(Y): # kJ/kg
    hmCO2=-393.51 # kJ/mol
    hmH2O=-285.83 # kJ/mol
    
    XCO2=Y[0]*1000/MW[0]
    XH2O=Y[1]*1000/(2*MW[1])
    return XCO2*hmCO2 + XH2O*hmH2O  

# Model for estimating pyrolysis products 
def pyrolysis(Y,T,Y_ch): # Neves, Thunman et al 2011.
    # Elements in tar
    tarC=1.005 + 1.9e-4 * T
    tarH=0.93 + 3.8e-4 * T
    tarO=0.92 - 2.2e-4 * T
    Ytar=np.array([tarC, tarH, tarO, 0, 0])*Ydaf
    
    # Elements in char
    charC=0.93 - 0.92*np.exp(-0.42e-2*T)
    charH=-0.41e-2 + 0.1*np.exp(-0.24e-2*T)
    charO=0.07 - 0.85*np.exp(-0.48e-2*T)
    Ych=np.array([charC, charH, charO, 0, 0])
    
    # Elements in gaseous species
    Ycxhy=np.array([24/28, 4/28, 0, 0, 0]) # longer hc's
    Ych4=np.array([MW[0], 4*MW[1], 0, 0, 0]) / (MW[0]+4*MW[1])
    Yco=np.array([MW[0], 0, MW[2], 0, 0]) / (MW[0]+MW[2])
    Yco2=np.array([MW[0], 0, 2*MW[2], 0, 0]) / (MW[0]+2*MW[2])
    Yh20=np.array([0, 2*MW[1], MW[2], 0, 0]) / (2*MW[1] + MW[2])
    Yh2=np.array([0,2*MW[1], 0, 0, 0]) / (2*MW[1])
    
    # Relation for H2/CO
    omega1=3e-4 + 0.0429/(1 + (T/632)**(-7.23))
    
    # Lower heating values of components
    
    






