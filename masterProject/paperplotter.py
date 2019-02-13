#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 15:58:22 2019

Script to generate plots for paper on lateral bed dispersion.

@author: gabgus
"""


import numpy as np, matplotlib.pyplot as plt

# -----------------------------------------------------------------------------

dp=0.95e-3
U0=0.96
hmf=0.45


rhop=2600.0
rho_steam=0.20230481999999997  
mu_air = 44.0e-6
mu_steam= 4.0422497000000007e-05

g=9.82

Ar=rho_steam*(rhop-rho_steam)*g*dp**3/(mu_steam**2)

def ReMf(Ar):
    return np.sqrt(27.2**2 + 0.0408*Ar) - 27.2

U_mf=ReMf(Ar)*mu_steam/(rho_steam*dp)


U_diff=U0-U_mf



nowall= [0.001248780546947, 0.00102684548646, np.nan, 
         0.000549463440657, 0.002143768143151, 0.002197853762626, 
         0.001273352342662, 0.001548487878214, 0.001498536656336, 
         0.001302645152242, 0.002997073312672, 0.002180804779755, 
         0.0035, 0.001948097653237, 0.002075362635415]

udiff_nowall=[0.306, 0.306, 0.306, 0.306, 0.306, 0.473, 0.473, 0.473, 0.473, 
              0.473, 0.639, 0.639, 0.639, 0.639, 0.639]

open30cm=[9.27735E-05, 9.83253E-05, 9.76313E-05, 9.08388E-05, 9.63018E-05, 
9.56205E-05, 0.000199055, 0.000210322, 0.000209326, 0.000131593, 0.000139412, 
0.000138481, 0.000279074, 0.000294833, 0.000293578, 0.000242868, 0.000256644, 
0.000255424]

udiff30cm =[0.306, 0.306, 0.306, 0.306, 0.306, 0.306, 0.473, 0.473, 0.473, 0.473, 
0.473, 0.473, 0.639, 0.639, 0.639, 0.639, 0.639, 0.639]







# Shi and Fan correlation

def shifan(udiff,hmf,dp,rhof,rhop,myf):
    Q = 0.46*(udiff*dp*rhof/myf)**(-0.21)*(hmf/dp)**0.24*((rhop-rhof)/rhof)**(-0.43)
    Dsr = Q*udiff*hmf
    return Dsr

udiff=np.linspace(0.42,0.543)

plt.figure()

#plt.plot(udiff,shifan(udiff,hmf,dp,rho_steam,rhop,mu_steam),label='Shi and Fan correlation')
plt.plot(U_diff,0.009,'o',label='porous plate')

plt.plot([U_diff,U_diff,U_diff,U_diff],[0.001,0.002,0.002,0.003],'o',label= 'nozzles')    
plt.plot([U_diff,U_diff,U_diff,U_diff],[0.0026,0.0024,0.0024,0.0021],'o',label= 'nozzles-qdata')

plt.plot(udiff_nowall,nowall,'o',label='exp-nowall')
plt.xlabel('$u-u_{mf}$ [m/s]')
plt.ylabel('$D_L$ $[m^2/s]$')
plt.legend()
#plt.savefig('./images/fig1.pdf')

plt.figure()
plt.plot(udiff30cm,open30cm,'o',label='exp-30cm')
plt.plot(0.3,0.0004,'o',label='nozzles under wall')

plt.plot([0.3],[0.002],'bo',label= 'nozzles')    
plt.plot([0.3],[0.0024],'bx',label= 'nozzles-qdata')
plt.plot(udiff_nowall,nowall,'o',label='exp-nowall')


plt.xlabel('$u-u_{mf}$ [m/s]')
plt.ylabel('$D_L$ $[m^2/s]$')
plt.legend()
#plt.savefig('./images/fig2.pdf')

plt.figure()
plt.plot(['HP','HP'],[0.001248780546947, 0.00102684548646],'o',)
plt.plot(['LP','LP'],[0.000549463440657, 0.002143768143151],'o',)
plt.plot(['HP','LP'],[0.00225,0.00225],'--')
#plt.savefig('./images/fig3.pdf')

















