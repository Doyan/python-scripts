#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  3 11:39:15 2018

@author: gabgus
"""
import numpy as np
# -------------------------
# Collection of fuels and corresponding gas yields

LHVi = np.array([283.0, 286.0, 0, 891.0, 1411.0, 0])     # kj/mol gas components
MWi  = [28.01, 2.016, 44.01, 16.043, 28.053, 18.015] # g/mol
MWi = np.array(MWi) / 1000 # kg / mol

MW_CHO = np.array([12.01, 1.007, 15.998])/1e3 # kg/mol


def saveCase(casedict,fnam):
    np.save(fnam,casedict)
    return

def readCase(fnam):
    data = np.load(fnam).flat[0]
    header = ''
    delim = ' '
    for k, v in data.items():
        header = header + delim + k
    print('Case: ' + data['name'] + '. -- Loaded variables:\n')
    print(header)
    return data, header
    
def getHeader(case):
    header = ''
    delim = ' '
    data = []
    for k, v in case.items():
        header = header + delim + k
        data.append(v)
    print(header)
    print(data)  


######## CHALMERS DFB Gasifier ############################3

case_chalmers = {}

# Input

# [co, h2, co2, ch4, c2h4, h2o]
ni_cg =  np.array([12.1, 7.5, 4.05, 4.1, 1.6, 0])
CHO = np.array([50.55,6.13,43.22]) / 100
SB = 0.7
ER = 1.2
LHV_fuel = 19.03 # MJ/kg daf
Ych = 0.16


# derived

CHO = CHO/CHO.sum()
 
yGas_mass = np.sum(ni_cg * MWi) # kg / kg daf
X_cg = ni_cg/ni_cg.sum()


MW_cg = np.sum(X_cg * MWi) # kg / mol gas
LHV_cg_mole = np.sum(X_cg*LHVi)/1000 # MJ/mol

LHV_cg_mass = LHV_cg_mole/MW_cg


nCHO_fuel = CHO[0:3] / MW_CHO
xCHO_fuel = nCHO_fuel/nCHO_fuel.sum()
MW_fuel = np.sum(xCHO_fuel*MW_CHO)


# flue gas and air

ni_fuel = CHO / MW_CHO

co2 = ni_fuel[0]
h2o = ni_fuel[1] / 2
o2 = co2 + h2o/2 - ni_fuel[2]/2

ni_fg = np.array([co2, h2o, o2*(ER-1), o2*ER*3.76])
ni_air = o2*ER + o2*ER*3.76

MW_O2 = 15.998*2/1000
MW_N2 = 14.0067*2/1000

m_air = o2*ER*MW_O2 + o2*ER*3.76*MW_N2




case_chalmers['name'] = 'Chalmers'
case_chalmers['ni'] = ni_cg
case_chalmers['X_cg'] = X_cg
case_chalmers['yGas_mass'] = yGas_mass
case_chalmers['yGas_mole'] = ni_cg.sum()
case_chalmers['CHO'] = CHO
case_chalmers['yield_char'] = Ych 
case_chalmers['SB'] = SB
case_chalmers['ER'] = ER
case_chalmers['AFR'] = m_air
case_chalmers['LHV_fuel'] = LHV_fuel
case_chalmers['LHV_cg_mass'] = LHV_cg_mass
case_chalmers['LHV_cg_mole'] = LHV_cg_mole   

saveCase(case_chalmers,'case_chalmers')




############### Anton bark values #######################################

case_anton_high = {}
case_anton_mid = {}
case_anton_low = {}

# volume percent
# high act high T
# low aght high T
# low act low T

# H2	CO	CO2	CH4	C2H2	C2H4	C2H6	C3H6
Xi = [[40.58, 16.03, 31.54, 8.31, 0.38, 2.88, 0.13, 0.16],
      [36.32, 20.48, 28.65, 9.94, 0.5, 3.62, 0.21, 0.27],
      [27.51, 25.92, 28.25, 12.16, 0.73, 4.45, 0.26, 0.73]
]


# kg gas (no tar)/kg fuel
yield_gas = [0.88, 0.7, 0.53] 

# elements CHO
CHO =[52.675, 5.7375, 37.625]


ychar= 0.2622875 # yield char kg/kgdaf
yvol = 0.703075 # yield volatiles kg/kgdaf	
yash = 0.0391125 # yield ash kg/kgdaf

HHVfuel = 20.50125 # MJ/kg
LHVfuel = 19.25125 # MJ/kg

# shoehorn in the heavier compounds as ethene
xi=[[],[],[]]
for i in range(3):
    xi[i] = Xi[i][0:4]
    xi[i].append(sum(Xi[i][4:])) 
    xi[i].append(0.0) # add spot for water
    xi[i] = np.array(xi[i]) # make array
    xi[i] = xi[i]/xi[i].sum() #  normalise











case_anton_high['name'] = 'Anton_high'








