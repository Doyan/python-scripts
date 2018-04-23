# -*- coding: utf-8 -*-
"""
Created on Fri Mar 23 13:32:55 2018

@author: Gabriel
"""

import numpy as np
from thermopy import nasa9polynomials as nasa9
db=nasa9.Database()
p_atm = 0.101325 # MPa

# ---------------------------------------------------------------------
# Input variables
xH2O_G = 0.1     # Moisture content in Gasification feedstock
xH2O_C = 0.2     # Moisture content in Combustion fuel

Xch=0.75      # char conversion
Xvol=0.1299     # fraction of heat from combusting volatiles

SB = 0.6    # Steam-Biomass ratio
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

Tref = 25   # Reference temp
TK=273.15   # for conversion to Kelvin

yield_gas = 0.91 # kg cold gas per kg fuel


# Fuel CHONS chalmers dfb

# Assume ni based on chalmers DFB at SB=0.7
# co h2 co2 ch4 c2h4  
ni_cg=np.array([12.1, 7.5, 4.05, 4.1, 1.5]) # mol/kg daf fuel

# Given product composition
# [co, h2, co2, ch4, c2h4, h2o]
X_gas_dry = np.array([0.358, 0.304, 0.166, 0.137, 0.03, 0])
X_gas_dry = X_gas_dry/X_gas_dry.sum()

# Elemental composition of fuel
Xc_fuel = np.array([1.0, 1.37, 0.61]) # normalised for carbon
X_fuel = Xc_fuel/Xc_fuel.sum()     # normalised to 1 mol fuel

# Molar weight of elements and fuel
MW_CHO = np.array([12.01, 1.007, 15.998])/1e3 # kg/mol
MW_fuel = np.sum(X_fuel * MW_CHO) # kg/mol


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
    
# -----------------------------------------------------------------------------
# 
    
# readCase


# Mass flows 
    
# get LHV of cold gas and from that massflow of cold gas
LHV_gas = LHVgas(X_gas_dry)
m_cg = Pgas/LHV_gas

# massflow of fuel needed to generate said gas
mf_g = m_cg/yield_gas

# massflow of fuel needed to generate the total needed heat, 
# boiler efficiency of 0.89
mf_c = (Pheat + Pgas)/nboil/LHV_f
nf_c  = mf_c/MW_fuel

# Convert wet basis moisture content to dry basis
xH2O_C = makeDryBasis(xH2O_C)
xH2O_G = makeDryBasis(xH2O_G)

# Fetch total amount of water from moisture in each chamber
mm_c = xH2O_C*mf_c
mm_g = xH2O_G*mf_g


# amount and comp of dry fluegas per mol fuel
X_fg, yfg, yAir = fg_comp(X_fuel,ER)
MW_fg = mixtureMW(X_fg,comps_flue)

MW_air = mixtureMW(X_air,comps_air)

m_fg = nf_c*yfg*MW_fg 
m_air = nf_c*yAir*MW_air

# Mass steam leaving with gas
steam_reacted =  6.5*MW_h2o/(96*MW_fuel) # kg steam reacted per kg fuel
ms_g = mf_g * (SB - steam_reacted)

# -----------------------------------------------------------------------------
# Heat calculations

# Heat needed due to 1kg moisture
qh2o=(sensHeat(h2ol,Tf) - sensHeat(h2ol,Tevap) - dHvap/MW_h2o/1e3) # kj/ kg h2o(l)

# Heat entering with 1 kg air
qair = senseHmix(comps_air,X_air,Tair) # kj/kg air

# Heat entering with 1 kg steam
qsteam = sensHeat(h2o, Ts)  # kj/kg steam

# Fetch good cp for gas and flue
Cp_g = meanCpmix(comps_gas,X_gas_dry,Tref,TG) # kj/kgK
Cp_fg = meanCpmix(comps_flue,X_fg,Tref,TC)
Cp_m = meanCp(h2o,Tref,TC)

# Cp for fuel
Cp_f = (5.46*(Tf+Tref+2*TK)/2 - 524.77)/1000 # kj/kgK fuel Gupta 2003

# Reaction enthalpy per kg daf fuel
dHr = 131.48 # kJ/mol  reaction: C(s) + H2O -> CO + H2

qrchar= pa_f[1]/C.molecular_weight *Xch*dHr # kj/kg

# Energy from combustion of char per kg daf fuel
qcchar = pa_f[1] * LHV_ch*1e3

# energy from combusting volatiles
qcvola = (LHV_f*1e3 - qcchar ) * Xvol

# energy from gasification char
bonus = pa_f[1]*(1-Xch)*LHV_ch

## -------------------- Actual source terms ------------------------------------

# Boundary conditions
Q_air = qair*m_air # kW
Q_steam = qsteam*mf_g*SB # kW
mcp_C = (m_fg*Cp_fg + mm_c*Cp_m) # kW/K 
mcp_G = (m_cg*Cp_g + mm_g*Cp_m) # kW/K

# Volumetric terms
Q_C = mf_c * Cp_f*(Tf-Tref) + mf_c*qcchar - mm_c*qh2o + mf_g*bonus + mf_c*qcvola # kW
Q_G = mf_g * Cp_f*(Tf-Tref) - mf_g*qrchar - mm_g*qh2o  # kW

# total balance to see if more radiation is needed
Qsum = Q_air + Q_steam + Q_C + Q_G - mcp_C*(TC-Tref) - mcp_G*(TG-Tref)

# ---- Printing ------------
print('\n\n===================================================================\n')
print('Final input terms:\n')

print('Q_G = {:.5} MW'.format(Q_G/1e3))
print('Q_C = {:.5} MW\n'.format(Q_C/1e3))

print('qair = {:.5} MW'.format(Q_air/1e3))
print('qsteam = {:.5} MW\n'.format(Q_steam/1e3))
print('mCP_flue = {:.5} kW/K'.format(mcp_C))
print('mCP_gas = {:.5} kW/K\n'.format(mcp_G))

print(Xvol)
print(Qsum)









