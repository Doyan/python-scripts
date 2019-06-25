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

dvalues=np.array([0.00175,0.002,0.002,0.002,0.01025,0.011,0.024])
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

u0_exp=np.array([0.488, 0.651, 0.815])
umf_exp=0.193

udiff_exp=u0_exp-umf_exp
uq_exp=u0_exp/umf_exp

exp_void_corr = 1.5625

# för nowall:  Lp: reps=2, för HP: reps=3
# för 30cm: reps=3 för båda 

nowallHP = [0.001957344729159, 0.001256181625365,0.003547687321601, 
          0.003303019230456,0.00155774343517,0.004159357549463, 
          0.005627366096332,0.002667866713111,0.005492798646202]

nowallLP=[0.001191593335453,0.002622558297209, 
          0.002749830774123,0.00159357851419,
          0.003605333681627,0.002538875072201]

nowallLP=np.array(nowallLP)* exp_void_corr
#nowallLP[1]=np.nan


gap30LP = [0.000111126746975,0.000117809893996,0.000116976443394,
           0.000160983247542,0.000170548814617,0.000169409714406, 
           0.000297110083034,0.000313962556596,0.000312470917859]

gap30LP=np.array(gap30LP)*exp_void_corr

gap30HP = [0.000113493526582,0.00012028526518,0.000119436297855, 
           0.000243512063063,0.000257294845833,0.000256077285517, 
           0.00034140209181,0.000360680828354,0.000359146016319]



# -----------------------------------------------------------------------------
# Plotting - large particles
#%%
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
#fig.savefig( path + 'Pictures/python/fig1.pdf')
#%%
# -------------------------------


fig, ax = plt.subplots(1,2,sharey=False,figsize=(12,5))

fig.suptitle('Experiment vs Simulation')

c = [0.3,0.5,0.9]
cd = [0.2,0.2,0.5]

r = [0.9,0.5,0.3]
rd = [0.5,0.2,0.2]


A=np.column_stack((np.repeat(uq_exp,2),nowallLP))
A1 = [A[n,:] for n in [0,2,4]]
A2 = [A[n,:] for n in [1,3,5]]



exp_x = [udiff_exp, uq_exp]
sim_x = [U_diff,Uq]

xlabels = ['$u-u_{mf}$ [m/s]','$u/u_{mf}$ [-]']
locs = ['upper left', 'upper right']

for i,axes in enumerate(ax):

    axes.plot(np.repeat(exp_x[i],2),nowallLP,'o',color=c,label='experiment')

    axes.plot(np.repeat(sim_x[i],4),dvalues[0:4],'s',color=cd,label= 'simulation-D-fitted')    
    axes.plot(np.repeat(sim_x[i],4),d_qsignal[0:4],'x',color=cd,label= 'simulation-q-method')

    
    axes.plot(np.repeat(exp_x[i],3),gap30LP,'o',color=r,label='experiment')

    axes.plot(np.repeat(sim_x[i],3),dvalues[1:4]*qvalues[1:4],'s',color=rd,label='simulation-D-fitted')
    axes.plot(np.repeat(sim_x[i],3),d_qsignal[1:4]*qvalues[1:4],'x',color=rd,label='simulation-q-method')

    

    axes.set_xlabel(xlabels[i],fontsize=11)
    axes.set_ylabel('$D_L$ $[m^2/s]$',fontsize=11)
    #axes.legend(loc=locs[i])
    axes.grid()
    
handles, labels = ax[1].get_legend_handles_labels()
fig.legend(handles, labels,loc='upper center' ,ncol=3,bbox_to_anchor=(0.529,0.93),columnspacing=15,fontsize=11)
fig.tight_layout(rect=(0,0,1,0.85))







#%% 

# -----------------------------------------------------------------------------


fig, ax = plt.subplots(1,2,sharey=False,figsize=(12,5))
fig.suptitle('Experiment vs Simulation - Sette')

# Experimental data - Sette
dp=0.30e-3

Ar=rho_steam*(rhop-rho_air)*g*dp**3/(mu_air**2)
U_mf=ReMf(Ar)*mu_steam/(rho_steam*dp)

sette_umf = 0.022
sette_udiff=[0.27,0.39,0.72,0.81,0.97]
sette_uq = (np.array(sette_udiff) + sette_umf) / sette_umf

h02=[0.00559,0.0100,0.0145,0.0157,0.0201]
h03=[0.0145,0.0168,0.0179,0.0190,0.0280]
h04=[0.019,0.0235,0.0280,0.0335,0.0414]

sette_x = [sette_udiff, sette_uq]
sim_x = [0.72, (0.72 + 0.022)/0.022]
xlabels = ['$u-u_{mf}$ [m/s]','$u/u_{mf}$ [-]']

for i,axes in enumerate(ax):
    axes.plot(sette_x[i],h02,'o',label='sette h0: 0.2m')
    axes.plot(sette_x[i],h03,'o',label='sette h0: 0.3m')
    axes.plot(sette_x[i],h04,'o',label='sette h0: 0.4m')
    axes.plot(sim_x[i],dvalues[6],'x',label='sim')
    axes.plot(sim_x[i],d_qsignal[6],'x',label='sim-qsignal')
    axes.set_xlabel(xlabels[i],fontsize=11)
    axes.set_ylabel('$D_L$ $[m^2/s]$',fontsize=11)
    axes.legend()






















