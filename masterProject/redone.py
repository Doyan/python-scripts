# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 18:45:19 2018

@author: Gabriel
"""
import numpy as np, matplotlib.pyplot as plt

import pandas as pd

from thermopy import nasa9polynomials as nasa9
db=nasa9.Database()

# -------------------------------------------------------------------------
# Table values
# -------------------------------------------------------------------------

# thermodynamic data for compounds through thermopy:
h2o = db.set_compound('h2o')
h2ol = db.set_compound('h2o(l)')
h2 = db.set_compound('h2') 
co = db.set_compound('carbon monoxide')
co2 = db.set_compound('carbon dioxide')
ch4 = db.set_compound('ch4')
c2h2 = db.set_compound('C2H2,acetylene')
c2h4 = db.set_compound('c2h4')
c2h6 = db.set_compound('c2h6')
c3h6 = db.set_compound('c3h6,propylene')
n2 = db.set_compound('n2')
o2 = db.set_compound('o2')
C = db.set_compound('c(gr)')

# component lists for different gas mixtures
comps_gas = [h2, co, co2, ch4, c2h2, c2h4, c2h6, c3h6]
comps_air = [o2, n2]
comps_flue = [co2, h2o, o2, n2]

# 0 deg C in Kelvin, for conversion purposes
TK= 273.15 # K

# molar composition of air 
x_air = np.array([1, 3.76])/4.76

# enthalpy of vaporisation for water.
dh_vap = h2o.enthalpy(373.15) - h2ol.enthalpy(373.15) # J/mol

# Heating value of anthracite from "Data och Diagram"
lhv_ch = 32.7                 # MJ/kg to represent char

# heating values for gas components
# from NREL answer to query:
# http://www.ieatask33.org/app/webroot/files/file/publications/HeatingValue.pdf

names_gas = ['h2', 'co', 'co2', 'ch4', 'c2h2', 'c2h4', 'c2h6', 'c3h6']
lhv_gas = [241.79, 282.9, 0, 802.71, 1256.9, 1323.2, 1428.83, 1926.1]     # kj/mol

# add lhv to thermopy objects for convenience
for i,comp in enumerate(comps_gas):
    comp.lhv = lhv_gas[i]

# molar masses through thermopy
mw_gas =[]
for comp in comps_gas:
    mw_gas.append(comp.molecular_weight)

# kg / mol molar mass of water
mw_h2o = h2ol.molecular_weight  

# -------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------

# Specific enthalpy 
# (thermopy compound, T in Celsius) -> kj / kg

def enthalpy(comp,T): 
    return comp.enthalpy(T+273.15)/comp.molecular_weight/1e3

# Specific enthalpy of mixture 
# (thermopy compound list, molar fraction list, T in Celsius) -> kJ/kg

def enthalpyMix(comps,X,T):
    h = []
    X= X / sum(X)
    for i, comp in enumerate(comps):
        h.append(X[i] * enthalpy(comp,T))
    return sum(h)

# Molar mass of mixture
# (thermopy compound list, molar fraction list,) -> kg/mol    

def mixtureMW(comps,X):
    MWmid = []
    X= X / sum(X)
    for i, comp in enumerate(comps):
        MWmid.append(X[i] * comp.molecular_weight)
    return sum(MWmid)


# Lower heating value of mixture
# (thermopy compound list, molar fraction list,) -> MJ/kg    

def mixtureLHV(comps,X):
    lhv = []
    X= X / sum(X)
    for i, comp in enumerate(comps):
        lhv.append(X[i] * comp.lhv )
    return sum(lhv) / mixtureMW(comps,X)/1000

# Mean specific heat capacity of thermopy compund over interval 
# (thermopy compound, start temp in deg C, end temp in deg C) -> kj/kgK    

def meanCp(comp,T1,T2):
    cp = []
    n = 20
    Tspan = np.linspace(min([T1,T2])+TK,max([T2,T1])+TK,n)
    
    for temp in Tspan:
        cp.append(comp.heat_capacity(temp))
    return sum(cp)/n/1e3/comp.molecular_weight


# Mean specific heat capacity of mixture over interval
# (thermopy compound list, molar fraction list, ...
#  ... start temp in deg C, end temp in deg C) -> kj/kgK

def meanCpmix(comps, X, T1, T2):
    cp = []
    X = X/X.sum()
    for i, comp in enumerate(comps):
        cp.append(X[i] * meanCp(comp,T1, T2))
    return sum(cp)    


# -------------------------------------------------------------------------
# Definition of base case
# -------------------------------------------------------------------------
class Case(object):
    """Base Case class to start from"""
    
    # Collection of given gas compositions in mol % 
    gas_order = ['h2', 'co','co2','ch4','c2h2','c2h4','c2h6','c3h6']
    gas_comps = {'Anton_high':[40.58, 16.03, 31.54, 8.31, 0.38, 2.88, 0.13, 0.16], 
              'Anton_mid': [36.32, 20.48, 28.65, 9.94, 0.5, 3.62, 0.21, 0.27],
              'Anton_low': [27.51, 25.92, 28.25, 12.16, 0.73, 4.45, 0.26, 0.73]
    }
    
    # Corresponding gas yields in kg/kgdaf      
    gas_yields = {'Anton_high': 0.88, 'Anton_mid': 0.7, 'Anton_low': 0.53}
    
    # -------------------------------------------------------------------------
    # Inputs 
    gascombo='Anton_mid'    # namestring for chosen gas comp and yield
    
    X_gas = np.array(gas_comps[gascombo])/100    # molar fraction, gas composition
    yield_gas = gas_yields[gascombo]                # kg gas/kgdaf
    
    # Fuel properties, change all at once!
    y_char = 0.2622875      # yield char kg/kgdaf
    y_vol = 0.703075        # yield volatiles kg/kgdaf	
    y_ash = 0.0391125       # yield ash kg/kgdaf 
    
    fuel_HHV = 20.50125     # MJ/kg
    fuel_LHV = 19.25125     # MJ/kg
    
    CHO = np.array([52.675, 5.7375, 37.625]) / 100 # mass %, major elements     
    
    # Temperatures
    TG = 850        # deg C
    TC = 850        # deg C
    Tair = 200      # deg C
    Tsteam = 200    # deg C
    
    # Loads
    P_gas = 10      # MW
    P_heat = 85     # MW
    
    # Velocities
    u_air = 1.2     # m/s
    u_steam = 1.2   # m/s
    
    # Dimensions
    L_chamber = 3   # m, Length
    W_chamber = 1.5 # m, Width
    
    # Moisture content, wet basis
    xh2o_g = 0.55
    xh2o_c = 0.55

    # Initialisation method to allow passing dict or keyword arguments 
    # to change individual attributes in order to create different cases
    def __init__(self, name='base', *change_data, **kwargs):
        self.name=name
        for dictionary in change_data:
            for key in dictionary:
                setattr(self, key, dictionary[key])
        for key in kwargs:
            setattr(self, key, kwargs[key])

# --------------------------------------------------------------------------
# Source term calculation (enclose in loop or make function?)
# --------------------------------------------------------------------------
case = Case()

lhv_cg = mixtureLHV(comps_gas,case.X_gas)   #  MJ/kg, lhv of produced gas

m_cg = case.P_gas / lhv_cg                  # kg/s, needed massflow cold gas    

m_gfuel = m_cg / case.yield_gas         # kg/s, needed massflow gasifier fuel







# qfuel_g = ygas*cp_gas(Tbed-Tref) + ych*cp_ch(T_bed - Tref) - cp_f*(Tf-Tref) 
# qh2o_g = xh2o_g*(enthalpy(h2o,TG) - enthalpy(h2o,Tf)) 

# sink_g = mfuel_g*qsink_g 

# sink_c = mh2o_c*qh2o + mfuel_c*qfuel + mchar_c*qchar
# swhole_c = mchar_c * LHVchar + (1-Xchar)*LHVchar 
# swhole_g = mchar_g * Xchar*dHchar


