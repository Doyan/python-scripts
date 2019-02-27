#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 15:58:22 2019

Script to generate plots for paper on lateral bed dispersion.

@author: gabgus
"""

chalmerspath = '/chalmers/users/gabgus/'
path = chalmerspath

import numpy as np, matplotlib.pyplot as plt
import matplotlib.colors as color

# -----------------------------------------------------------------------------
# Properties
dp=0.95e-3
U0=0.96
hmf=0.45


rhop=2600.0
rho_steam=0.20230481999999997
rho_air= 0.29 
mu_air = 4.553e-5
mu_steam= 4.0422497000000007e-05

g=9.82

# -----------------------------------------------------------------------------
# Simulation values

cases = [1,2,3,4,6,7,8]

dvalues=np.array([0.001,0.002,0.002,0.002,0.01025,0.011,0.024])
qvalues=np.array([1.0,0.2,0.2,0.2,1.0,0.7,1.0])

k_qsignal=np.array([3799.33,3588.61,3495.82,3162.17,10732.36,15862.86,19272.01])
alpha_s =np.array([0.4994 ,0.4983, 0.4985, 0.4986,0.5022, 0.5010, 0.3013]) 

d_qsignal=k_qsignal/alpha_s/rhop/1200






# -----------------------------------------------------------------------------
# Umf calculation

Ar=rho_steam*(rhop-rho_steam)*g*dp**3/(mu_steam**2)

def ReMf(Ar):
    return np.sqrt(27.2**2 + 0.0408*Ar) - 27.2

U_mf=ReMf(Ar)*mu_steam/(rho_steam*dp)


U_diff=U0-U_mf

Uq=U0/U_mf


# Shi and Fan correlation

def shifan(udiff,hmf,dp,rhof,rhop,myf):
    Q = 0.46*(udiff*dp*rhof/myf)**(-0.21)*(hmf/dp)**0.24*((rhop-rhof)/rhof)**(-0.43)
    Dsr = Q*udiff*hmf
    return Dsr


udiff=np.linspace(0.42,0.543)

# -----------------------------------------------------------------------------
# Experimental data - large particles 

u0_exp=np.array([0.499, 0.666, 0.832])
umf_exp=0.193

udiff_exp=u0_exp-umf_exp
uq_exp=u0_exp/umf_exp

# för nowall:  Lp: reps=2, för HP: reps=3
# för 30cm: reps=3 för båda 

nowallHP = [0.001248780546947, 0.00102684548646, np.nan, 
          0.002197853762626, 0.001273352342662, 0.001548487878214, 
          0.002997073312672, 0.002180804779755, 0.0035]

nowallLP=[0.000549463440657, 0.002143768143151, 
          0.001498536656336, 0.001302645152242,
          0.001948097653237, 0.002075362635415]
nowallLP[1]=np.nan


gap30LP = [9.08388E-05, 9.63018E-05, 9.56205E-05,
           0.000131593, 0.000139412, 0.000138481, 
           0.000242868, 0.000256644, 0.000255424]

gap30HP = [9.27735E-05, 9.83253E-05, 9.76313E-05, 
           0.000199055, 0.000210322, 0.000209326, 
           0.000279074, 0.000294833, 0.000293578]

# -----------------------------------------------------------------------------
# Plotting - large particles

fig, ax = plt.subplots(1,2,sharey=True,figsize=(10,5))

fig.suptitle('Porous plate included')
#plt.plot(udiff,shifan(udiff,hmf,dp,rho_steam,rhop,mu_steam),label='Shi and Fan correlation')
ax[0].plot(np.repeat(U_diff,2),dvalues[4:6],'o',label='porous plate')
ax[0].plot(np.repeat(U_diff,4),dvalues[0:4],'o',label= 'nozzles')    
ax[0].plot(np.repeat(U_diff,4),d_qsignal[0:4],'o',label= 'nozzles-qdata')

ax[0].plot(np.repeat(udiff_exp,2),nowallLP,'o',label='exp-LP')

ax[0].set_xlabel('$u-u_{mf}$ [m/s]')
ax[0].set_ylabel('$D_L$ $[m^2/s]$')
ax[0].legend()
ax[0].grid()

ax[1].plot(np.repeat(Uq,2),dvalues[4:6],'o',label='porous plate')
ax[1].plot(np.repeat(Uq,4),dvalues[0:4],'o',label= 'nozzles')    
ax[1].plot(np.repeat(Uq,4),d_qsignal[0:4],'o',label= 'nozzles-qdata')

ax[1].plot(np.repeat(uq_exp,2),nowallLP,'o',label='exp-LP')

ax[1].set_xlabel('$u/u_{mf}$ [-]')
ax[1].set_ylabel('$D_L$ $[m^2/s]$')
ax[1].grid()
fig.savefig( path + 'Pictures/python/fig1.pdf')

# -------------------------------

fig, ax = plt.subplots(1,2,sharey=True,figsize=(10,5))

fig.suptitle('Zoom in on nozzles')

c = [0.3,0.5,0.9]
cd = [0.2,0.2,0.5]

r = [0.9,0.5,0.3]
rd = [0.5,0.2,0.2]

ax[0].plot(np.repeat(udiff_exp,2),nowallLP,'o',color=c,label='exp-nowall-LP')

ax[0].plot(np.repeat(udiff_exp,3),gap30LP,'o',color=r,label='exp-30cm-LP')

ax[0].plot(np.repeat(U_diff,3),dvalues[1:4]*qvalues[1:4],'o',color=rd,label='nozzles under wall')
ax[0].plot(np.repeat(U_diff,3),d_qsignal[1:4]*qvalues[1:4],'x',color=rd,label='nozzles under wall qdata')

ax[0].plot(np.repeat(U_diff,4),dvalues[0:4],'o',color=cd,label= 'nozzles')    
ax[0].plot(np.repeat(U_diff,4),d_qsignal[0:4],'x',color=cd,label= 'nozzles-qdata')


ax[0].set_xlabel('$u-u_{mf}$ [m/s]')
ax[0].set_ylabel('$D_L$ $[m^2/s]$')
#ax[0].legend(loc='upper left')
ax[0].grid()

ax[1].plot(np.repeat(uq_exp,2),nowallLP,
  'o',color=c,label='exp-nowall-LP')

ax[1].plot(np.repeat(uq_exp,3),gap30LP,'o',color=r,label='exp-30cm-LP')

ax[1].plot(np.repeat(Uq,3),dvalues[1:4]*qvalues[1:4],'o',color=rd,label='nozzles under wall')
ax[1].plot(np.repeat(Uq,3),d_qsignal[1:4]*qvalues[1:4],'x',color=rd,label='nozzles under wall qdata')

ax[1].plot(np.repeat(Uq,4),dvalues[0:4],'o',color=cd,label= 'nozzles')    
ax[1].plot(np.repeat(Uq,4),d_qsignal[0:4],'x',color=cd, label= 'nozzles-qdata')


ax[1].set_xlabel('$u/u_{mf}$ [-]')
ax[1].set_ylabel('$D_L$ $[m^2/s]$')
ax[1].grid()

handles, labels = ax[1].get_legend_handles_labels()
fig.legend(handles, labels,loc='upper center' ,ncol=3,bbox_to_anchor=(0.535,0.93),columnspacing=14)
fig.tight_layout(rect=(0,0,1,0.85))

#plt.savefig('./images/fig2.pdf')

# -----------------------------------------------------------------------------

# Experimental data - Sette
dp=0.30e-3

Ar=rho_steam*(rhop-rho_air)*g*dp**3/(mu_air**2)
U_mf=ReMf(Ar)*mu_steam/(rho_steam*dp)

sette_udiff=[0.27,0.39,0.72,0.81,0.97]
h02=[0.00559,0.0100,0.0145,0.0157,0.0201]
h03=[0.0145,0.0168,0.0179,0.0190,0.0280]
h04=[0.019,0.0235,0.0280,0.0335,0.0414]


plt.figure()
plt.plot(sette_udiff,h02,'o',label='sette h0: 0.2m')
plt.plot(sette_udiff,h03,'o',label='sette h0: 0.3m')
plt.plot(sette_udiff,h04,'o',label='sette h0: 0.4m')
plt.plot([0.72],dvalues[6],'x',label='sim')
plt.plot([0.72],d_qsignal[6],'x',label='sim-qsignal')
plt.ylabel('$D_L$ $[m^2/s]$')
plt.xlabel('$u-u_{mf}$ [m/s]')
plt.legend()



