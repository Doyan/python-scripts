#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  3 13:52:51 2018

@author: gabgus
"""
# ---------------------------------------------------------------------
import numpy as np
from thermopy import nasa9polynomials as nasa9
db=nasa9.Database()
p_atm = 0.101325 # MPa

# read fuel - gas variables from file
def readCase(fnam):
    data = np.load(fnam).flat[0]
    header = ''
    delim = ' '
    for k, v in data.items():
        header = header + delim + k
    print('Case: ' + data['name'] + '. -- Loaded variables:\n')
    print(header)
    return data, header

# ---------------------------------------------------------------------
# Input variables
xH2O_G = 0.1     # Moisture content in Gasification feedstock
xH2O_C = 0.2     # Moisture content in Combustion fuel

Xch=0.75      # char conversion

ER = 1.2    # Equivalence ratio

# Effects in the gasifier and combustor
Pgas = 15 # MW
Pheat = 85 # MW
nboil= 0.89 # boiler efficiency

# Feed composition, dry ash free basis, phyllis id #1269
#CHONS = np.array([50.7, 6.12, 42.93, 0.20, 0.03])
pa_f = np.array([0.86, 0.14])  # daf, Vol, fC -- ash 0.4% dry basis
LHV_f = 17.62                 # MJ/kg

# Temperatures, in deg C
TG = 850
TC = 850     # Same temp in both chambers ... Unreasonable...?

Ts = 200     # steam temp                    -- check sensitivity
Tf = 80      # fuel temp, total guess        -- check sensitivity
Tair = 200   # air temp,                     -- check sensitivity

Tevap=100   # Boiling point of water
Tpyro=400   # Temp at which fuel dissociates


Tref = 25   # Reference temp
TK=273.15   # for conversion to Kelvin


case,header = readCase('case_chalmers.npy')

# -------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Table values

h2o = db.set_compound('h2o')
h2ol = db.set_compound('h2o(l)')
h2 = db.set_compound('h2') 
co = db.set_compound('carbon monoxide')
co2 = db.set_compound('carbon dioxide')
ch4 = db.set_compound('ch4')
c2h4 = db.set_compound('c2h4')
n2 = db.set_compound('n2')
o2 = db.set_compound('o2')
C = db.set_compound('c(gr)')

comps_gas = [co, h2, co2, ch4, c2h4, h2o]
comps_air = [o2, n2]
comps_flue = [co2, h2o, o2, n2]

X_air = np.array([1, 3.76])/4.76

dHvap = h2o.enthalpy(373.15) - h2ol.enthalpy(373.15) # J/mol
MW_h2o = h2ol.molecular_weight  # kg / mol

# Heating values
LHV_ch = 32.0                 # MJ/kg Heating value of char
LHVi = [283.0, 286.0, 0, 891.0, 1411.0, 0]     # kj/mol gas components

#--------------functions for dealing with composition and conversions----------

#  mol /mol -> mol/(mol gas), mol gas/mol fuel, mol 
def fg_comp(X_fuel,ER):
    co2 = X_fuel[0]
    h2o = X_fuel[1]/2
    o2 = co2 + h2o/2 - X_fuel[2]/2
    X = np.array([co2, h2o, o2*(ER-1), o2*ER*3.76])
    yFg = X.sum()
    return X/yFg, yFg, sum([o2*ER, o2*ER*3.76])

# 
def mixtureMW(X,gcomps):
    MWmid = []
    for i, comps in enumerate(gcomps):
        MWmid.append(X[i] * comps.molecular_weight)
    return sum(MWmid)

# molar fraction => MJ/kg
def LHVgas(Xdry):
    LHV = Xdry * LHVi
    MW = mixtureMW(Xdry, comps_gas)
    return LHV.sum()/1e3/MW

def makeDryBasis(xh2o):
    return xh2o/(1-xh2o)

# --------- Sensible heat and Cp for compounds and mixtures ------------------
# comp refers to compunds with polynomials in the nasa9 database

def meanCp(comp,T1,T2):
    cpsum=0
    n = 20
    Tspan = np.linspace(min([T1,T2])+TK,max([T2,T1])+TK,n)
    
    for i, temp in enumerate(Tspan):
        cpsum += comp.heat_capacity(temp) /comp.molecular_weight / n
    return cpsum /1e3

def meanCpmix(comps, X, T1, T2):
    cp = []
    X = X/X.sum()
    for i, comp in enumerate(comps):
        cp.append(X[i] * meanCp(comp,T1, T2))
    return sum(cp)    

def sensHeat(comp,T):
    h = (comp.enthalpy(T+TK) - comp.enthalpy(Tref+TK)) / comp.molecular_weight
    return h / 1e3

def senseHmix(comps,X,T):
    h = []
    X= X / sum(X)
    for i, comp in enumerate(comps):
        h.append(X[i] * sensHeat(comp,T))
    return sum(h)

# -------------------------------------------------------------------------    
# mass flows 

Pfuel_c = (Pgas+Pheat)/nboil

m_cg = Pgas / case['LHV_cg_mass']

mf_g = m_cg / case['yGas_mass']
m_steam = mf_g*case['SB']


mf_c = Pfuel_c/case['LHV_fuel']
m_air = mf_c*case['AFR']

# Convert wet basis moisture content to dry basis
xH2O_C = makeDryBasis(xH2O_C)
xH2O_G = makeDryBasis(xH2O_G)

mm_c = xH2O_C*mf_c
mm_g = xH2O_G*mf_g

mchar_c = mf_c*case['yield_char'] 
mchar_g = mf_g*case['yield_char'] 




# -------------------------------------------------------------------------
# Heat sinks 

# Cp for fuel
Cp_f = (5.46*(Tf+Tref+2*TK)/2 - 524.77)/1000 # kj/kgK fuel Gupta 2003

# Cp for gas
Cp_g = meanCpmix(comps_gas,case['X_cg'],Tref,TG)

# Cp for char j/kgK fuel Gupta 2003
Cp_ch = -0.0038*((Tpyro + TG +2*TK)/2)**2 + 5.98*(Tpyro + TG + 2*TK)/2 - 795.28
Cp_ch = Cp_ch/1000 # kj/kgK

# needed to heat 1kg of daf fuel to dissociation
qfuel = Cp_f*(Tf - Tref) - Cp_f * (Tpyro -Tref) #kj/kg

# needed heat to heat 1kg char to Bed temp
qchar = Cp_ch*(Tpyro - Tref) - Cp_ch*(TG - Tref) #kj/kg

# needed to heat and evaporate 1kg moisture 
qh2o=(sensHeat(h2ol,Tf) - sensHeat(h2ol,Tevap) - dHvap/MW_h2o/1e3) #kj/kg



















