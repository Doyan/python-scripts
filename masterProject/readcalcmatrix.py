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

# Test matrix formulation
N = 2**3
Fukt = np.resize([0,1],N)
char = np.resize([0,0,1,1],N)
D = np.resize([0,0,0,0,1,1,1,1],N)

imatrix = np.column_stack((Fukt,char,D))

Fukt = np.choose(Fukt,[0.08,0.7])
char = np.choose(char,[0,1])
D = np.choose(D,[20,100])

matrix = np.column_stack((Fukt,char,D))


    
# ---------------------------------------------------------------------
# Input variables
xH2O_G = 0.08     # Moisture content in Gasification feedstock
xH2O_C = 0.5     # Moisture content in Combustion fuel

Xch=1.0      # char conversion

partAir = 0.6 #  60 percent of total air used as primary air

# Effects in the gasifier and combustor
Pgas = 15 # MW
Pheat = 85 # MW
gasheat = 3 # MW heat needed per kgdaf gasifier fuel. 

compensate = True
bonus = True

Fudgecorr = -1.78

fnam = "test.csv"

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
Tsteam = 200 # steam temp

Tevap=100   # Boiling point of water
Tpyro=500   # Temp at which fuel dissociates


Tref = 25   # Reference temp
TK = 273.15   # for conversion to Kelvin


case,header = readCase('case_chalmers.npy')

ychar = case["yield_char"] 

# -------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Table values

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

comps_gas = [co, h2, co2, ch4, c2h4, h2o]
comps_air = [o2, n2]
comps_flue = [co2, h2o, o2, n2]

X_air = np.array([1, 3.76])/4.76

dHvap = h2o.enthalpy(373.15) - h2ol.enthalpy(373.15) # J/mol
MW_h2o = h2ol.molecular_weight  # kg / mol

# Heating values
LHV_ch = 32.0                 # MJ/kg Heating value of char
LHVi = [283.0, 286.0, 0, 891.0, 1411.0, 0]     # kj/mol gas components

# Molar masses
MWi  = [28.01, 2.016, 44.01, 16.043, 28.053, 18.015] # g/mol
MWi = np.array(MWi) / 1000 # kg / mol


#--------------functions for dealing with composition and conversions----------
def charGas(mchar,Xchar):
    nchar = mchar/C.molecular_weight
    nco = nchar*Xchar
    nh2 = nchar*Xchar
    ni_char = np.array([nco, nh2, 0, 0, 0, 0])
    return ni_char


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
directory = '/scratch/gabgus/geometry/present/sourcefiles/'
print('\n============================= Start of run ===============================\n')



D_array = np.array([5,10,15,20,30,40,60,80,100]) # matrix[:,2]
moist_array = np.ones_like(D_array)*0.6
char_array = np.ones_like(D_array)*0.1 # matrix[:,1]

name_array=[]
for i in np.arange(N):
    name_array.append(str(imatrix[i][0]) + str(imatrix[i][1]) + str(imatrix[i][2]))
    
name_array=['5','10','15','20','30','40','60','80','100']

N = len(name_array)
for i in np.arange(N):
    xH2O_G = moist_array[i]
    xH2O_C = 0.5
    
    Xch = char_array[i]
    height = 0.5
    
    casestring = name_array[i]
    fnam = directory + 'source' + casestring + '.csv'
    # Convert wet basis moisture content to dry basis
    xH2O_C = makeDryBasis(xH2O_C)
    xH2O_G = makeDryBasis(xH2O_G)
    
    # char conversion effects on gas
    ni_char = charGas(case['yield_char'],Xch)
    ni_cg = case['ni'] + ni_char
    
    X_cg = ni_cg/ni_cg.sum()    
    
    LHV_cg = LHVgas(ni_cg/ni_cg.sum())
    
    yGas_mass = np.sum(ni_cg * MWi) 
        
    
    # mass flows 
    
    m_cg = Pgas / LHV_cg
    
    mf_g = m_cg / yGas_mass
    m_steam = mf_g*case['SB']
    mchar_g = mf_g*case['yield_char'] 
    mm_g = xH2O_G*mf_g
      
    # ------------------------------------------------------------------------
    # Heat sinks 
    
    # Cp for fuel 
    Cp_f = (5.46*(Tf+Tref+2*TK)/2 - 524.77)/1000 # kj/kgK fuel Gupta 2003 https://doi.org/10.1016/S0016-2361(02)00398-8
    
    # Cp for gas
    Cp_g = meanCpmix(comps_gas,X_cg,Tref,TG)
    
    # Cp for char j/kgK fuel Gupta 2003
    Cp_ch = -0.0038*((Tpyro + TG +2*TK)/2)**2 + 5.98*(Tpyro + TG + 2*TK)/2 - 795.28
    Cp_ch = Cp_ch/1000 # kj/kgK
    
    # needed to heat 1kg of daf fuel to dissociation
    qfuel = Cp_f*(Tf - Tref) - Cp_f * (Tpyro -Tref) #kj/kg
    
    # needed heat to heat 1kg char to Bed temp
    qchar = Cp_ch*(Tpyro - Tref) - Cp_ch*(TG - Tref) #kj/kg
    
    # needed to heat and evaporate 1kg moisture 
    qh2o= (sensHeat(h2ol,Tf) - sensHeat(h2ol,Tevap) - dHvap/MW_h2o/1e3) #kj/kg
    
    qsteam= (sensHeat(h2o,Tsteam) - sensHeat(h2o,TG))
    
    
    # Source terms gas
    
    Sm_g = mchar_g*Xch
    
    n_added =  mchar_g*Xch/C.molecular_weight
    nR = n_added * 8.314 # * T/P / Ag
    
    # Reaction enthalpy per kg daf fuel
    dHr = 131.48 # kJ/mol  reaction: C(s) + H2O -> CO + H2
    Sq_g = -n_added * dHr / 1e3 # MW 
    
    # heating and drying
    qsink_g = qfuel + xH2O_G*qh2o + ychar*qchar
    
    sink_g = mf_g*qsink_g/1e3  

    # fuel energy needed - bonus energy from nonconverted char. assumes 100% efficient heat transfer to gasifier... 
    qsink_c = qfuel + xH2O_C*qh2o + ychar*qchar
    qair = senseHmix(comps_air,X_air,Tair) - senseHmix(comps_air,X_air,TC) 
    
    Pfuel_c = (Pheat)/nboil
    mfc0 = Pfuel_c/LHV_f
    
    bonusheat = 0
    bonusmass = 0
    gasheat = 0
    if compensate:
        gasheat = -(sink_g + m_steam*qsteam/1e3 + Sq_g)
    if bonus:
        bonusmass = mchar_g*(1-Xch)
        bonusair = bonusmass*(32/12) + 3.76*bonusmass*(28/12)
        bonusheat = bonusmass*LHV_ch
    
    net_correction =  (gasheat - bonusmass*LHV_ch)
    
    if net_correction < 0:
        mf_c = mfc0
    else:
        mf_c = mfc0 + net_correction/LHV_ch/ychar
        
    
    m_air = partAir*mf_c*case['AFR']
    
    mm_c = xH2O_C*mf_c
    
    mchar_c = mf_c*case['yield_char'] 
    
    
    # Source terms Combustor
    Sm_c = mchar_c
    
    Sq_c = mchar_c * LHV_ch + bonusheat + Fudgecorr # MW
    
    sink_c = mf_c * qsink_c/1e3 
    
    Stot = Sq_c+Sq_g+sink_c+sink_g
    # BC values
    
    # density of air
    rho_air = p_atm*1e6*mixtureMW(X_air,comps_air)/(8.314*(Tair+TK))
    rho_s = 8.3
    
    
    
    data = [1, int(casestring), Sm_g, Sm_c, nR, Sq_c*1e6, Sq_g*1e6, sink_c*1e6, sink_g*1e6, m_air/rho_air, m_steam*2.2,Tsteam+TK,Tair+TK, height,mf_c,mf_g, LHV_f,D_array[i]]
    head = 'row, case, Sm_g, Sm_c, nR, Sq_c, Sq_g, sink_c, sink_g, Vdot_air, Vdot_steam, Tsteam, Tair, height, mf_c, mf_g, LHV_fuel,D'

    with open(fnam, 'w') as f:
        f.write(head)
        f.write('\n')
        for i,val in enumerate(data):    
            f.write(str(val) + ',')
        f.close()
    
    print('\n\n===================================================================\n')
    print('Final input terms ' + casestring + ':\n')

    print('QG = {:.5} MW'.format(Sq_g))
    print('QC = {:.5} MW\n'.format(Sq_c))
    print('sink_c = {:.5} MW'.format(sink_c))
    print('sink_g = {:.5} MW\n'.format(sink_g))
    
    print('Correction = {:.5}'.format(net_correction))
    print('Total = {:.5} MW\n'.format(Stot))

    print('Total fuel = {:.4} + {:.4} = {:.5} kg'.format(mf_c,mf_g,mf_c + mf_g))
    print('Total fuel energy for heat = {:.5} MW'.format(mf_c*LHV_f*nboil))

    print('Total mass char = {:.5} kg'.format(Sm_c))

    print('Air speed : {:.4} m/s'.format(m_air/rho_air/39))
    print('Steam speed : {:.4} m/s'.format(m_steam*2.2/4.5))



