#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 13:36:53 2018

@author: gabgus NOTE: PYTHON2
"""
import sys
sys.path.append('/chalmers/users/gabgus/miniconda2/envs/p3/lib/python3.5/site-packages/cantera/')


import numpy as np
#import matplotlib.pyplot as plt
import cantera as ct

gas = ct.Solution('gri30.xml')
iCH4  = gas.species_index('CH4')
iCO  = gas.species_index('CO')
iCO2 = gas.species_index('CO2')
iH2O = gas.species_index('H2O')
iH2 = gas.species_index('H2')
iO2 = gas.species_index('O2')
indices = [iCH4,iCO,iCO2,iH2O,iH2]


#--------Key variables----------------------------------------------------------------------
alpha_w = 0.1 # Moisture content
SB = 0.63        # Steam-Biomass ratio

# Temperatures, in deg C
Ts = 400    # -- check sensitivity
Tbed = 850
Tfeed = 80  # total guess -- check sensitivity

# for conversion to Kelvin
TrefK=298.15
TK=273.15

# ------------------------ Composition and energy in gas -----------------------
# Set state of product gas at reference , Xi based on hofbauer article
gas.X = {'CH4': 0.09377, 'CO': 0.2036, 'CO2': 0.1225, 'H2': 0.2323, 'H2O': 0.3478}

gas.TP =  TrefK, ct.one_atm

# Save enthalpy of product gas at reference
hgasref=gas.enthalpy_mass/1e6

# Save hot gas enthalpy
gas.TP = Tbed+TK, ct.one_atm
hgas = gas.enthalpy_mass/1e6
 
# calculate needed Oxygen for combustion
O2= gas.X[iCO]*0.5 + gas.X[iH2]*0.5 + gas.X[iCH4]*2

gas.X = {'CH4': 0.09377, 'CO': 0.2036, 'CO2': 0.1225, 'H2': 0.2323, 'H2O': 0.3478, 'O2':O2}
gas.TP =  TrefK, ct.one_atm

gas.equilibrate('TP')
hfgas=gas.enthalpy_mass/1e6

HHV = -(hfgas-hgasref)
LHV = HHV - gas[iH2O].Y[0]*2264.705/1e3

# Char yield
Ych = 0.106+2.43*np.exp(-0.66e-2*Tbed)

# Mass flows
             # MW
mGas = Pgas/LHV         # kg/s
mBM = mGas/(1-Ych+SB)   # kg/s
mS = mBM*SB             # kg/s
mCh = mBM*Ych           # kg/s

#----------------------- Composition of char and BM --------------------------
compBM = np.array([50.8, 41.67, 7.51])/100 # c o h dry basis

MW = gas['C','O','H','H2O'].molecular_weights
# From Thunman article, T in deg C
charC=0.93 - 0.92*np.exp(-0.42e-2*Tbed)
charH=-0.41e-2 + 0.1*np.exp(-0.24e-2*Tbed)
charO=0.07 - 0.85*np.exp(-0.48e-2*Tbed)
compCh=np.array([charC, charO, charH])


# Boie correlation for estimation of HHVs for char and BM 
# could choose more specialised correlations but this suffices for now
def HHV_boie(Y): # Y in mass percent
    Y=Y*100
    return (351.7*Y[0] + 1419.33*Y[2] -110.95*Y[1])/1e3


# Function for calculating enthalpy of combustion products.
def h_post(Y): # kJ/kg
    hmCO2=-393.51 # kJ/mol
    hmH2O=-285.83 # kJ/mol
    
    XCO2=Y[0]*1000/MW[0]
    XH2O=Y[2]*1000/(2*MW[2])
    
    return (XCO2*hmCO2 + XH2O*hmH2O)/1e3, XH2O*MW[3]/1000 

# Heating values and post combustion values for BM and char
HHVch = HHV_boie(compCh)
HHVbm = HHV_boie(compBM)

hfbm, h2obm = h_post(compBM)
hfch, h2och = h_post(compCh)

# enthalpy of dry bm and char at reference
hbmref = hfbm + HHVbm
hchref = hfch + HHVch

# sensible heat of bm and char
cpDrybm=0.00486*(TrefK + Tfeed+TK)/2 - 0.21293 # kJ/kg/K - Peduzzi Boisssonnet 2016

hbm=hbmref + cpDrybm*(Tfeed-25)/1e3
    
x=np.linspace(TrefK,Tbed+TK)
hchar = hchref + np.trapz(1.4*x + 688,x)/1e6

hH2Ofeed=335.05/1e3

hbm=(1-alpha_w)*hbm + alpha_w*hH2Ofeed

hs400 = 3278/1e3    # J/kg steam -- from Elliot Lira book


Qbed = mGas*hgas + mCh*hchar - mS*hs400 - mBM*hbm

print(Qbed/mBM)



gas['H2', 'CO', 'CO2','CH4','H2O'].standard_enthalpies_RT













