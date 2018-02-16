# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 13:46:10 2017

@author: Gabriel

Investigation of flame temperature in crossflow configuration 
"""

import numpy as np, matplotlib as mpl, cantera as ct
#from scipy.optimize import fsolve
plt = mpl.pyplot
colors = mpl.colors




# -----------------------------------------------------------------------------
# Input Data

# Hypothetical area
D = 10.0e-3     # m
A=np.pi*D**2/4

p=ct.one_atm
tin_f= 473.0
tin_o= 300.0
mdot_f0=3.24631218245e-06 /A
mdot_o=mdot_f0*3.9
comp_o='O2:1.0, N2:3.76' 
comp_f='CH4:1'
Lf=0.413
width= Lf*3  # 3 times flame length in a)

gas = ct.Solution('gri30.xml', 'gri30_mix')
gas.TP = gas.T, p

iN2  = gas.species_index('N2')
iO2  = gas.species_index('O2')
iH2O = gas.species_index('H2O')
iF = gas.species_index('CH4')
iCO2= gas.species_index('CO2')


# -----------------------------------------------------------------------------
# Functions

def counterFlow(p,tin_f,tin_o,mdot_o,mdot_f,comp_o,comp_f,width,gas):
    # Create an object representing the counterflow flame configuration,
    # which consists of a fuel inlet on the left, the flow in the middle,
    # and the oxidizer inlet on the right.
    f = ct.CounterflowDiffusionFlame(gas, width=width)
    # Set the state of the two inlets
    f.fuel_inlet.mdot = mdot_f
    f.fuel_inlet.X = comp_f
    f.fuel_inlet.T = tin_f
    f.oxidizer_inlet.mdot = mdot_o
    f.oxidizer_inlet.X = comp_o
    f.oxidizer_inlet.T = tin_o
    #Solving setups
    loglevel = 1  # amount of diagnostic output (0 to 5)
    f.set_boundary_emissivities(0.0, 0.0)
    f.radiation_enabled = False
    f.set_refine_criteria(ratio=4, slope=0.2, curve=0.3, prune=0.04)  
    #Solve counterflow
    f.solve(loglevel, auto=True)
    return f

def printDataSheet(fnam,f,param,value):
    file = open(fnam,"w")
    file.write("# Studied params: " + param + "=" + str(value))
    file.write("\n")
    #file.write("\n# Temp: " + str(t))
    
    Z=f.mixture_fraction('C')           
    z = len(Z)     
    
    for i in range(0,z):
        file.write(str("{:011.8f}".format(Z[i])).ljust(16))
        file.write(str("{:011.8f}".format(f.T[i])).ljust(16))
        file.write(str("{:011.8f}".format(f.Y[iF,i])).ljust(16))
        file.write(str("{:011.8f}".format(f.Y[iO2,i])).ljust(16))
        file.write(str("{:011.8f}".format(f.Y[iCO2,i])).ljust(16))
        file.write(str("{:011.8f}".format(f.Y[iH2O,i])).ljust(16))
        file.write(str("{:011.8f}".format(f.Y[iN2,i])).ljust(16))
        file write(str("{:011.8f}".format(f.).ljust(16))
        file.write("\n")
    
    file.close()   
    return 

# -----------------------------------------------------------------------------
# Main calculation

f = counterFlow(p,tin_f,tin_o,mdot_o,mdot_f0,comp_o,comp_f,width,gas)
printDataSheet('../ref.txt',f,'afr:3.9, width/Lf=',3)


## parameter variation width
#w=np.array([0.009, 0.01, 0.1, 1, 10.0, 100.0, 1000.0])
##
##for val in w:
##    print 'w= ' + str(val)
##    width=Lf*val
##    f = counterFlow(p,tin_f,tin_o,mdot_o,mdot_f,comp_o,comp_f,width,gas)
##    fnam= 'width' + str(val) + '.txt'    
##    printDataSheet(fnam,f,'afr=3.9   width/Lf',val)
#
## parameter variation air to fuel ratio
#
#afr=np.array([0.01, 0.1, 1, 10, 100])
#width=Lf*3
#for val in afr:
#    print 'afr= ' + str(val)
#    mdot_o=mdot_f0*val
#    f = counterFlow(p,tin_f,tin_o,mdot_o,mdot_f0,comp_o,comp_f,width,gas)
#    fnam= 'afr' + str(val) + '.txt'    
#    printDataSheet(fnam,f,'width/Lf=3 afr= ',val)
#
## parameter variation massflow fuel
#afr=5
#m=np.array([0.0001, 5.7, 5.8])
#for val in m:
#    print 'm= ' + str(val)
#    mdot_f=mdot_f0*val
#    mdot_o=mdot_f*afr
#    f = counterFlow(p,tin_f,tin_o,mdot_o,mdot_f,comp_o,comp_f,width,gas)
#    fnam= 'm' + str(val) + '.txt'    
#    printDataSheet(fnam,f,'m= ',val)
#  
#
#



















