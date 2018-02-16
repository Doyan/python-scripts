# -*- coding: utf-8 -*-
"""
Created on Mon Aug 28 11:16:05 2017

@author: fistler
"""
import matplotlib.pyplot as plt
import cantera as ct
import numpy as np
from scipy.integrate import ode

def func(t,Z,a,b):
    T = Z[0]
    Y = Z[1:]
    ...
    dYdt = 
    dTdt = 
    dZdt = np.append(dTdt, dYdt)
    return dZdt

    
Y0 = 
T0 = 
Z0 = np.append(T0,Y0)
...
t0 = 
t1 = 

solver = ode(constP_ODE).set_integrator('vode')
solver.set_initial_value(Z0,t0).set_f_params(a, b)
    
tsol = []
Tsol = []
NOsol = []

while solver.successful() and np.any(solver.t < t1):
    solver.integrate(t1,step=True)
    tsol.append(solver.t)
    Tsol.append(solver.y[0])
    NOsol.append(solver.y[idNO+1])
    


