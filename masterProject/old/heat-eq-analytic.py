# -*- coding: utf-8 -*-
"""
Created on Wed Sep  5 17:01:16 2018

Script to plot solutions to a 1 dimensional heat equation 
solved analytically through separation of variables.

Corresponds to conduction through a rod of length L with its ends kept at 
fixed temperatures Th and Tc, initially haviing a uniform temperature of T0

Necessary material properties are k, rho and cp

# Th
#    .
#        .
#            .
#                Tc
# 
# |------ L ------|

@author: Gabriel
"""
import numpy as np, matplotlib.pyplot as plt

import pandas as pd
import os

# -----------------------------------------------------------------------------
# Inputs

Th = 845 + 273.15   # K
Tc = 805 + 273.15   # K
T0 = 825 + 273.15   # K

L = 0.51    # m

k = 20  # W/mK

t = 23.6


# -------------------------------------
# Properties

rho_s = 2610.0  # kg / m3 
rho_g = 1.0     # kg / m3

cp_g = 1.17e3    # J / kgK
cp_s = 1.23e3     # J / kgK

eps = 0.45  # void fraction

rho = eps*rho_g + (1 - eps)*rho_s # kg / m3
cp = eps*cp_g + (1-eps)*cp_s # J / kgK

def S(x,L):
    a = (Tc -Th)/L
    b = Th
    return a*x + b
    
def v(x,t,alpha,L):
    N = int(2000*L)
    S = 0 
    pi = np.pi
    for k in range(N):
        n = k + 1
        psum = L/(n*pi)*((T0 - Th) - (T0 - Tc)*(-1)**n)*np.sin(n*pi/L*x)*np.exp(-n**2*pi**2/L**2*alpha*t)
        S += psum
    return 2/L * S

def u(x,t,k,L=L):
    alpha = k/(rho*cp)    
    return v(x,t,alpha,L) + S(x,L)

N = 200
x = np.linspace(0,L,N)

T = u(x,t,k)

ax=plt.figure(1,[6,7])

plt.xlabel('L [m]')
plt.ylabel('T [$\degree C$]')
plt.title('Solution to the 1D heat equation at t = ' +str(t) + 's')
plt.plot(x,T-273.15)
plt.plot([0,L],[Th-273.15,Tc-273.15],'k--')
plt.grid()

