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

#--------Key variables----------------------------------------------------------------------
alpha_w = 0.1 # Moisture content
SB = 0.63        # Steam-Biomass ratio


# Temperatures, in deg C
Ts = 400    # -- check sensitivity
Tbed = 850
Tfeed = 80  # total guess -- check sensitivity
TrefK=298.15
TK=273.15

# Composition of biomass C,O,H,Ash,H2O
Ydaf = np.array([49.5, 44.6, 5.90, 0, 0]) / 100 # from data for pine bark (N,S etc lumped in O)
Ywet = ((1-alpha_w)*Ydaf + (alpha_w)*np.array([0,0,0,0,1]))

#--- Nice Table values ---------------------------------------------------------------------------
MW=np.array([12.0107, 15.999, 1.00794, 0, 18.01528]) # g/mol -- Molar mass, Ash omitted
dhVap=2264.705 # kj/kg latent heat of water

hs400 = 3278    # kj/kg steam -- from Elliot Lira book
hs850 = 4278.2  # kj/kg steam -- from Elliot Lira book



# Lazy two-way index dict
i = {'tar': 0, 'cxhy': 1, 'ch4': 2, 'co': 3, 'co2': 4, 'h2o': 5, 'h2': 6}
temp=dict(i)
for item in enumerate(temp):
    i.update({i[item[1]]:item[1]})
del temp

#--- Relations ------------------------------------------------------------------------------


def HHV_bm(Y):  # Changdong Sheng 2005 for dry biomass with lumped N and O
    Y=Y*100
    return (-1.3675 + 0.3137*Y[0] + 0.7009*Y[2] + 0.0318*Y[1])*1e3

def HHV_boie(Y):
    Y=Y*100
    return 351.7*Y[0] + 1419.33*Y[2] -110.95*Y[1]

def HHV_ch(Y): # IGT 1978  used for char (made for coal) 
    Y=Y*100
    return (0.341*Y[0] + 1.323*Y[2] + 0.0685 - 0.0153*Y[3] - 0.1194*Y[1])*1e3

def cpBM(alpha,T): # Temperature in K, alpha is moisture content
    cpdry=0.00486*T - 0.21293 # kJ/kg/K - Peduzzi Boisssonnet 2016
    cph2o=4.18
    cpwet=(1-alpha)*cpdry + alpha*cph2o
    return cpwet

def cpCh(T): # T in K
    return 1.4*T + 688

# Function for calculating enthalpy of combustion products.
def h_post(Y): # kJ/kg
    hmCO2=-393.51 # kJ/mol
    hmH2O=-285.83 # kJ/mol
    
    XCO2=Y[0]*1000/MW[0]
    XH2O=Y[2]*1000/(2*MW[2])
    
    return XCO2*hmCO2 + XH2O*hmH2O, XH2O*MW[4]/1000 

def YchF(T):
    return 0.106+2.43*np.exp(-0.66e-2*T)

def compPch(T):
    charC=0.93 - 0.92*np.exp(-0.42e-2*T)
    charH=-0.41e-2 + 0.1*np.exp(-0.24e-2*T)
    charO=0.07 - 0.85*np.exp(-0.48e-2*T)
    Ych=np.array([charC, charO, charH, 0, 0])
    return Ych

def compPtar(T,Y):
    tarC=1.05 + 1.9e-4 * T
    tarO=0.92 - 2.2e-4 * T
    tarH=0.93 + 3.8e-4 * T
    Ytar=np.array([tarC, tarO, tarH, 0, 0])*Y
    return Ytar

# function for elemental mass composition from molecular formula as list [C,O,H] 
def elemComp(formula):
    formula=np.array(formula)
    Mcomp=np.sum(formula*MW[0:3])
    Ycomp=(formula*MW[0:3])/Mcomp
    return Mcomp, Ycomp
    
# Elements in gaseous species
Mcxhy, Ycxhy = elemComp([2,0,4])
Mch4, Ych4 = elemComp([1,0,4])
Mco, Yco = elemComp([1,1,0])
Mco2, Yco2 = elemComp([1,2,0])
Mh2o, Yh2o = elemComp([0,1,2])
Mh2, Yh2 = elemComp([0,0,2])

Mprod=np.array([Mch4, Mco, Mco2, Mh2o, Mh2])
elemProd=np.array([Ych4, Yco, Yco2, Yh2o, Yh2]) 



# Model for estimating pyrolysis products 
def pyrolysis(Y,T,Y_ch): 
    # T in deg C Neves, Thunman et al 2011.  
    
    Ytar=compPtar(T,Y)  # Elements in tar    
    Ych=compPch(T)      # Elements in char   
    
    # Relation for H2/CO
    omega1=3.0e-4 + 0.0429/(1 + (T/632.0)**(-7.23))
    
    # Lower heating value of gas and components
    LHVg = -6.23 + 2.47e-2*T
    LHVcxhy = 47.2
    LHVco = 10.1
    LHVch4 = 50.0
    LHVh2 = 120     
    
    # Relation for hydrogen vs temperature
    omega2=1.145*(1-np.exp(-0.11e-2*T))**9.384
    
    # columns for matrix equation
    col1 = np.append(Ytar[0:3],[0, 0, LHVg, 0])
    col2 = np.append(Ycxhy[0:3],[0, 0, LHVcxhy, 0])
    col3 = np.append(Ych4[0:3],[0, -1, LHVch4, 0])
    col4 = np.append(Yco[0:3],[-omega1, 0.146, LHVco, 0])
    col5 = np.append(Yco2[0:3],[0, 0, 0, 0])
    col6 = np.append(Yh2o[0:3],[0, 0, LHVg, 0])
    col7 = np.append(Yh2[0:3],[1, 0, LHVh2, 1])
    
    A=np.column_stack((col1,col2,col3,col4,col5,col6,col7))
    
    Y_gas=np.sum(Y[0:3])-Y_ch*np.sum(Ych[0:3]) # total gas yield
    
    b=np.append(np.resize(Y-Ych*Y_ch,3),[0, 2.18e-4, Y_gas*LHVg, omega2])
    Yprod=np.linalg.solve(A,b)
    
    return Ytar, Ych, LHVg, Yprod # LHV in MJ/kg

#------------ Main calculation ------------------------------------------------
# First iteration - chosen product composition
    # ch4 co co2 h2o h2
Xprod=np.array([13.4, 29.1, 17.5, 0, 33.2])/100 # volume fraction
alpha_w_p=0.332 # moisture content, v/v
Xprod_wet = Xprod*(1-alpha_w_p) + alpha_w_p*np.array([0,0,0,1,0])

# Cantera gas initialisation
gas = ct.Solution('gri30.xml')
iCH4  = gas.species_index('CH4')
iCO  = gas.species_index('CO')
iCO2 = gas.species_index('CO2')
iH2O = gas.species_index('H2O')
iH2 = gas.species_index('H2')
indices = [iCH4,iCO,iCO2,iH2O,iH2]

# Set state of gas
X=np.zeros_like(gas.X)
np.put(X,indices,Xprod_wet)
gas.TPX = Tbed+TK, ct.one_atm, X

compGas=np.zeros(3)
compGas[0]=gas.elemental_mass_fraction('C')
compGas[1]=gas.elemental_mass_fraction('O')
compGas[2]=gas.elemental_mass_fraction('H')

hfgas, pwgas=h_post(compGas)

# Fetch char composition
compCh=compPch(Tbed)

# Enthalpy of combustion products and amount of  combustion water
#hfgas, pwgas = h_post(compProd)     # kj/kg, kg -- for gas 
hfbm, pwbm = h_post(Ydaf)           # kj/kg, kg -- for biomass
hfch, pwch = h_post(compCh)         # kj/kg, kg -- for char

# Enthalpy of unburnt compounds at standard state (kj/kg)
hbm = hfbm + HHV_bm(Ydaf)
hch = hfch + HHV_ch(compCh)
hgas=gas.enthalpy_mass/1000


# adding the sensible heat to char and bm
x=np.linspace(TrefK,Tbed+TK)
hchar=hch+np.trapz(cpCh(x),x)/1000

x=np.linspace(TrefK,Tfeed+TK)
hbiomass=hbm+np.trapz(cpBM(alpha_w,x),x)

# Energy balance
qbed=(1-YchF(Tbed)+SB)*hgas + YchF(Tbed)*hchar - SB*hs400 -hbiomass











