# -*- coding: utf-8 -*-
"""
Created on Thu Sep  6 11:43:33 2018

@author: Gabriel
"""
import numpy as np

from thermopy import nasa9polynomials as nasa9
db=nasa9.Database()

# -----------------------------------------------------------------------------
# scrap calculation to get loads for empiric data

# 0 deg C in Kelvin, for conversion purposes
TK= 273.15 # K
Tref = 25# deg C

# Pressure at one atm
p_atm = 101325 # Pa

# molar Gas constant 
R=8.3144


CHONS = np.array([0.035, 0.518, 0.06, 0.381, 0.0054, 0.0009])
CHONS_daf = CHONS[1:]/CHONS[1:].sum()

lhv_fuel = 20.456 # MW/kgdaf

MW_CHONS=np.array([0.01201 , 0.001007, 0.015998, 0.0140067, 0.032])



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
H = db.set_compound('H')
O = db.set_compound('O')
naphtalene=db.set_compound('c10h8')
naphtalene.elements[0] = ('C',10) # thermopy has parsed the nasa polys wrong for 2 digit subscripts it seems.

# component lists for different gas mixtures
comps_gas = [h2, co, co2, ch4, c2h2, c2h4, c2h6, c3h6]
comps_air = [co2, h2o, o2, n2]
comps_flue = [co2, h2o, o2, n2]

# -----------------------------------------------------------------------------

X_CHONS = CHONS_daf/MW_CHONS # mol element/kgdaf 

n_co2 = X_CHONS[0] # mol co2 /kgdaf

n_h2o = X_CHONS[1]/2 # mol h2o/kgdaf

n_so2 = X_CHONS[4]

n_o2 = (X_CHONS[0] + X_CHONS[1]/4 - X_CHONS[2]/2 + X_CHONS[4]) # mol o2 /kgdaf

n_n2 = n_o2*3.76 + X_CHONS[3]/2    # mol n2/kgdaf

lt = (n_o2 + n_n2)* R*(TK)/p_atm # theoretical Nm3 air/kgdaf


#SEEMS WRONG? need o2conc in dry fluegas -> need total amount of flue gas
o2conc=np.array([0.058, 0.048]) # vol % o2 in moist flue gas?


xh2o_wet = 0.55
xh2o_daf = xh2o_wet /(1- xh2o_wet) 

n_moisture = xh2o_daf/ h2o.molecular_weight

gt = (n_co2 + n_n2 + n_so2)* R*(TK)/p_atm # theoretical Nm3 gas / kgdaf

# need to remake with dry o2conc
m = 1.0 + gt/lt * o2conc/(0.21 - o2conc) # air factor 

lv = lt*m  # real Nm3 air/kgdaf 

ltot = np.array([39.2241, 69.4282])/3.6 # Nm3 total air 

mfuel= ltot/lv # kgdaf 

Qfuel = mfuel * lhv_fuel


print(m)
print(Qfuel)
#%%
# Scap calc just for fun to calculate Ar number 
# and see if we're in the B or D particle region

g= 9.82
dp=0.95e-3
rho_g = 1.2
rho_p = 2610
mu_g =  44.0e-6

Ar = g*rho_g*(rho_p - rho_g)*dp**3/mu_g**2
print((Ar)**(1/3))




