# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 18:45:19 2018

@author: Gabriel
"""

# (gen input)

# read input

# calc quantities

# calc source terms


# Specific enthalpy (thermopy compound, T in Celsius) -> kj / kg
def enthalpy(comp,T): 
    return comp.enthalpy(T+273.15)/comp.molecular_weight/1e3

# Specific enthalpy of mixture (thermopy compound list, molar fraction list, T in Celsius) -> kJ/kg
def enthalpyMix(comps,X,T):
    h = []
    X= X / sum(X)
    for i, comp in enumerate(comps):
        h.append(X[i] * enthalpy(comp,T))
    return sum(h)



# qfuel_g = ygas*cp_gas(Tbed-Tref) + ych*cp_ch(T_bed - Tref) - cp_f*(Tf-Tref) 
# qh2o_g = xh2o_g*( 



# sink_g = mfuel_g*qsink_g 

# sink_c = mh2o_c*qh2o + mfuel_c*qfuel + mchar_c*qchar
# swhole_c = mchar_c * LHVchar + (1-Xchar)*LHVchar 
# swhole_g = mchar_g * Xchar*dHchar


