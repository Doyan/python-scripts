# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 11:31:16 2018

@author: Gabriel
"""

import numpy as np
import cantera as ct
from pyXSteam.XSteam import XSteam
st = XSteam(XSteam.UNIT_SYSTEM_MKS)

gas = ct.Solution('gri30.xml')
iCH4  = gas.species_index('CH4')
iCO  = gas.species_index('CO')
iCO2 = gas.species_index('CO2')
iH2O = gas.species_index('H2O')
iH2 = gas.species_index('H2')
iO2 = gas.species_index('O2')
iC2H4 = gas.species_index('C2H4')

#import matplotlib.pyplot as plt
# -----------------------------------------------------------------------------

# Key variables
xH2O_G = 0.1     # Moisture content in Gasification feedstock
xH2O_C = 0.2     # Moisture content in Combustion fuel

Xch=1      # char conversion
Xvol=0      # fraction of heat from combusting volatiles

SB = 0.7    # Steam-Biomass ratio
ER = 1.2
# Effects in the gasifier and combustor
Pgas = 15 # MW
Pheat = 85 # MW

# Feed composition, dry ash free basis, phyllis id #1269
CHONS = np.array([50.7, 6.12, 42.93, 0.20, 0.03])
paf = np.array([81.83, 18.17])  # daf, Vol, fC -- ash 0.4% dry basis
LHV_f = 17.62                   # MJ/kg

# Temperatures, in deg C

TG = 850
TC = 850     # Same temp in both chambers ... Unreasonable...?

Ts = 400     # steam temp                    -- check sensitivity
Tf = 80      # fuel temp, total guess        -- check sensitivity
Tair = 375   # air temp,                     -- check sensitivity
Tstack = 200 # Stack temp for combustion    -- check sensitivity
Troom = 15   # room temp for fractional Fg loss --||--

Tevap=100   # Boiling point of water

Tref = 25   # Reference temp
TK=273.15   # for conversion to Kelvin

# Pressure for easy steamtable
p=1.0       # bar 

# Assume ni based on chalmers DFB at SB=0.7
# co h2 co2 ch4 c2h4  
ni=np.array([12.1, 7.5, 4.05, 4.1, 1.5]) # mol/kg daf fuel
# -----------------------------------------------------------------------------
# Table values
LHVi = [283.0, 286.0, 0, 891.0, 1411.0]     # kj/mol
MWi  = [28.01, 2.016, 44.01, 16.043, 28.053] # g/mol

MWH2O = 18.015

dHvap=st.hV_t(Tevap)-st.hL_t(Tevap)

cpH2OL = (st.Cp_pt(p,Tref) + st.Cp_pt(p,Tevap-0.5))/2
cpH2OV = np.array([st.Cp_pt(p,Tevap), st.Cp_pt(p,Ts), st.Cp_pt(p,TG)])

cpf = (5.46*(Tf+Tref+2*TK)/2 - 524.77)/1000 # kj/kg K Gupta 2003

dHrchar=131.48 # heat of reaction for c + h2o -> co + h2

def cpGasCt(T,gas): # degC, ni/ntot -> kj/kg
    gas.TP=T+TK,ct.one_atm
    return gas.cp_mass/1e3

# ----------------------------------------------------------------------------- 
# Mass Calculations

# First, deal with the char conversion
Ychar=paf[1] # amount of char/kgdaf, assume = fixed carbon

nch= Ychar*10/12.01 # coarse assumption: all char = C_(s) 

# more coarse assumptions: char only turns into co and h2
nco= nch*Xch
nh2= nch*Xch

ni[0]+=nco
ni[1]+=nco
 
# Final yield and composition of gas after char conversion done
Ygas = np.sum(ni*MWi/1000) 
Xgas = ni/ni.sum()

# Calculate Heating value of gas
LHVgas=np.sum(LHVi*ni)/Ygas/1000 # MJ/kg cold gas

# Mass flows
mgas=Pgas/LHVgas    # kg cold gas
mf_G=mgas/Ygas        # kg fuel

# -----------------------------------------------------------------------------
# Heat calculations

# Inclusion of pyrolytic water in gas composition for cp calculation
ni_pw = np.append(ni,0.2e3/18.015) # assume worst case 0.2kg/kgdaf
Ygas_pw = np.sum(ni_pw*np.append(MWi,MWH2O)/1000)

gas['CO','H2','CO2','CH4','C2H4','H2O'].X = ni_pw
cpGas = (cpGasCt(Tref,gas) + cpGasCt(TG,gas))/2

# Individual heats for gasifier

qfuel = (cpGas*(TG-Tref) - cpf*(Tf-Tref)) # kj/kgdaf

# heat and evaporate moisture
cpS =(cpH2OV[0] + cpH2OV[2])/2
qmoist = xH2O_G*(cpH2OL*(Tevap-Tref) - cpH2OL*(Tf-Tref) + dHvap + cpS*(TG-Tevap)) 

# heat steam from Ts to TG
cpS2 = (cpH2OV[0] + cpH2OV[1])/2
qsteam = (cpS*(TG - Tref) - cpS2*(Ts-Tref))*SB

# heat needed for char conversion reaction
qrchar=nch*Xch*dHrchar


qGtot=qfuel + qmoist + qsteam + qrchar
QG = qGtot*mf_G


# --- Printing-----------------------------------------------------------------
print('---------- Gasifier Heating demands in MJ/kg daf fuel ---------- \n')
print('fuel: {:.5}  moisture: {:.5}  steam: {:.5} reaction: {:.5}'.format(qfuel/1e3,qmoist/1e3,qsteam/1e3,qrchar/1e3))
print('---------------------------------------------------------------')
print('total: {:.5}\n'.format(qGtot/1e3))
print('m_fuel = {:.5} kg/s'.format(mf_G))
print('LHV_gas = {:.5} MJ/kg cold gas\n'.format(LHVgas))

print('Final Gasifier Heat demand = {:.5} MW\n'.format(QG/1e3))


# =============================================================================
# =============================================================================
# Combustion chamber -- Most calculations on kg daf fuel basis

LHVchar = 32.0 # MJ/kg daf

# Have to estimate flue gas composition
MWchons = gas['C','H','O','N'].molecular_weights
MWchons = np.append(MWchons,32.065)

Xfuel = CHONS*10/MWchons # g/(g/mol) = mol

# neglect sulfur present
# c -> co2, h -> h2o, n -> n2
# co2 h20 n2 - comp dry fluegas
ni_fg_dry = np.array([Xfuel[0], (Xfuel[1])/2.0, Xfuel[3]/2.0]) # mol/kg daf

# needed o2
o2 = ni_fg_dry[0] + ni_fg_dry[1]/2.0 - Xfuel[2]/2.0 # mol/kg daf

# o2 n2 - comp air
ni_air_stoich = np.array([o2, o2*3.72])
ni_air = ni_air_stoich * ER 
MW_air = gas['O2', 'N2'].molecular_weights 

yAir = np.sum(ni_air*MW_air/1000.0)

# co2 h2o n2 o2
ni_fg = ni_fg_dry + np.array([0, 0, ni_air[1]])
ni_fg = np.append(ni_fg,ni_air[0] - ni_air_stoich[0])

ni_fg_wet = ni_fg  + np.array([0,xH2O_C*1e3/MWH2O, 0, 0])

MW_fg = gas['CO2','H2O','N2','O2'].molecular_weights
yFlue = np.sum(ni_fg*MW_fg/1000.0)
yFlue_wet = np.sum(ni_fg_wet*MW_fg/1000.0)

#------------------------------------------------------------------------------

# Cp for flue gas
gas['CO2','H2O','N2','O2'].X = ni_fg
cp_fg = (cpGasCt(Tref,gas) + cpGasCt(TC,gas))/2

gas['CO2','H2O','N2','O2'].X = ni_fg_wet
cp_fg_wet = (cpGasCt(Tref,gas) + cpGasCt(TC,gas))/2

# cp for air
gas['N2', 'O2'].X = ni_air
cp_air = (cpGasCt(Tref,gas) + cpGasCt(TC,gas))/2

# mass calculation

gas['CO2','H2O','N2','O2'].X = ni_fg_wet
gas.TP = Tstack + TK, p
hstack = gas.enthalpy_mass/1e3

gas.TP = Tref+TK, p
href = gas.enthalpy_mass/1e3

gas['N2', 'O2'].X = ni_air
gas.TP = Troom + TK, p

hl = gas.enthalpy_mass/1e3

gas.TP = Tref+TK, p
hlref = gas.enthalpy_mass/1e3
 
Fg = (yFlue_wet*(hstack - href) - yAir*(hl - hlref))/(LHV_f*1e3)


Fr = 0.0065
nu_b = 1 - Fr - Fg
Qf_C = (Pheat + Pgas)/nu_b
mf_C = Qf_C/LHV_f 

# mf_C = P_C/(yGas*(hgi -ho))
# Hi + yAir*(hl - hlref) = yGas*(hgi -hgref)
# Hi + yAir *cp_air*(Tair - Tref) = yGas*cp(Tgi- Tref)
# Tgi = VL/(yGas*cpgas) + Tref
# mf_C = QHeat/(yGas*Cpgas*(Tgi -Tstack))




# -----------------------------------------------------------------------------
# Heat demand calculations

qC_char = Ychar/100*(1 - Xch) * cp_fg*(TC-TG) # kJ/kg daf GASIFICATION fuel - heat needed to heat char from gsification section

qC_air = yAir * (cp_air*(TC - Tref) - cp_air*(Tair - Tref))
qC_moist = xH2O_C*(cpH2OL*(Tevap-Tref) - cpH2OL*(Tf-Tref) + dHvap + cpS*(TC-Tevap))
qC_fuel = (cp_fg*(TG-Tref) - cpf*(Tf-Tref))

heatC_char = Ychar/100 * LHVchar

qCtot = qC_air + qC_moist + qC_fuel
net_C = heatC_char - qCtot/1e3

r_GC= mf_G/mf_C
bonus = Ychar/100*(1 - Xch)*r_GC*(LHVchar - cp_fg*(TC-TG)/1e3)
QC = mf_C*(net_C + bonus)


# --- Printing-----------------------------------------------------------------
print('========================================================================')
print('\n\n---------- Combustor Heating demands in MJ/kg daf fuel ---------- \n')
print('fuel: {:.5}  moisture: {:.5}  air: {:.5}'.format(qC_fuel/1e3,qC_moist/1e3,qC_air/1e3))
print('---------------------------------------------------------------')
print('total: {:.5}\n'.format(qCtot/1e3))

print('heat from char: {:.5} MJ/kg daf'.format(heatC_char))
print('net extra from gasification char: {:.5} MJ/kg daf'.format(bonus))
print('total amount of fuel: {:.5} kg/s\n'.format(mf_C))

print('Net Excess heat: {:.5} MW'.format(QC))


mCP_flue= yFlue_wet*mf_C*cp_fg_wet
qair = yAir*mf_C*cp_air*(Tair - Tref)

# need new Ygas and cpGas with included moisture

ni_wet = ni_pw + np.array([0,0,0,0,0,xH2O_G/MWH2O])
Ygas_wet = np.sum(ni_wet*np.append(MWi,MWH2O)/1000)

gas['CO','H2','CO2','CH4','C2H4','H2O'].X = ni_wet
cpGas = (cpGasCt(Tref,gas) + cpGasCt(TG,gas))/2

mCP_gas= Ygas_wet*mf_G*cpGas

qsteam = SB*mf_G*cpS2*(Ts - Tref)

vf_sand=0.40
cpbed = 0.78*vf_sand + 1*(1-vf_sand)
density= 2650*vf_sand + 0.72*(1-vf_sand)


print('\n\n===================================================================\n')
print('Final input terms:\n')

print('QG = {:.5} MW'.format(QG/1e3))
print('QC = {:.5} MW\n'.format(QC))

print('mCP_flue = {:.5} kW/K'.format(mCP_flue))
print('qair = {:.5} kW'.format(qair))
print('mCP_gas = {:.5} kW/K'.format(mCP_gas))
print('qsteam = {:.5} kW\n'.format(qsteam))

print('bed cp = {:.5} kJ/kgK'.format(cpbed))
print('bed density = {:.5} kg/m3'.format(density))


