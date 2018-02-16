# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 13:10:41 2017

@author: fistler
"""

# Counterflow function
import cantera as ct

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

# Input parameters
p =   # pressure
...

# Create the gas object used to evaluate all thermodynamic, kinetic, and
# transport properties.
gas = ct.Solution('gri30.xml', 'gri30_mix')
gas.TP = gas.T, p

f = counterFlow(p,tin_f,tin_o,mdot_o,mdot_f,comp_o,comp_f,width,gas)

Z = f.mixture_fraction('C')
T = f.T
grid = f.flame.grid


#Print data in data sheet
file = open("???.txt","w")
file.write("# Header "+ str(x))
file.write("\n# Temp: " + str(T))

for i in range(0,z):
    file.write(str("{:011.8f}".format(Z[i])).ljust(16))
    file.write("\n")
 
file.close()   

#Read data sheet
data  = np.genfromtxt('?????.txt') #Input 
x = data[:,0]
y = data[:,1]

