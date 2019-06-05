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

cases = [1,2,3,
         4,6,7,
         8,9,10]

dvalues=np.array([0.00165, 0.002, 0.002,
                  0.002, 0.0098, 0.011,
                  0.0226, 0.00855, 0.0138])

# 1 2,3,4,6 7 8 9 10

kd_values= np.array([0.00075, np.nan, np.nan,
                     np.nan, 0.0023, np.nan, 
                     0.00655, 0.00585, 0.0079])

qvalues=np.array([1.0,0.2,0.2,
                  0.2,1.0,0.7,
                  1.0,1.0,1.0])

d_qsignal=np.array([0.00243839, 0.00230824, 0.00224765, 
                    0.00203272, 0.006174958, 0.0101482, 
                    0.019210092, 0.009364, 0.016194])

DT_method=[np.nan,np.nan,np.nan,
           np.nan,np.nan,np.nan,
           np.nan,0.0058,np.nan]

#------------------------------------------------------------------------------
# Experimental data - large particles 
exp_void_corr = 1.5625

u0_exp=np.array([0.488, 0.651, 0.815])
umf_exp=0.193

udiff_exp=u0_exp-umf_exp
uq_exp=u0_exp/umf_exp

# för nowall:  Lp: reps=2, för HP: reps=3
# för 30cm: reps=3 för båda 

nowallHP = [0.001957344729159, 0.001256181625365,0.003547687321601, 
          0.003303019230456,0.00155774343517,0.004159357549463, 
          0.005627366096332,0.002667866713111,0.005492798646202]

nowallLP=[0.002622558297209,0.001191593335453, 
          0.002749830774123,0.00159357851419,
          0.003605333681627,0.002538875072201]




gap30LP = [0.000111126746975,0.000117809893996,0.000116976443394,
           0.000160983247542,0.000170548814617,0.000169409714406, 
           0.000297110083034,0.000313962556596,0.000312470917859]

gap30HP = [0.000113493526582,0.00012028526518,0.000119436297855, 
           0.000243512063063,0.000257294845833,0.000256077285517, 
           0.00034140209181,0.000360680828354,0.000359146016319]


nowallLP = np.array(nowallLP)*exp_void_corr
gap30LP = np.array(gap30LP)*exp_void_corr

mean_exp = nowallLP.reshape(3,2).mean(axis=1)
err_exp = np.abs(nowallLP.reshape(3,2).T - mean_exp)
# -----------------------------------------------------------------------------
# Umf calculation

Ar=rho_steam*(rhop-rho_steam)*g*dp**3/(mu_steam**2)

def ReMf(Ar):
    return np.sqrt(27.2**2 + 0.0408*Ar) - 27.2

U_mf=ReMf(Ar)*mu_steam/(rho_steam*dp)

U_diff=U0-U_mf

Uq=U0/U_mf

# -----------------------------------------------------------------------------
# colors NOT USED
c = [0.3,0.5,0.9]
cd = [0.2,0.2,0.5]

r = [0.9,0.5,0.3]
rd = [0.5,0.2,0.2]



#%%
# -----------------------------------------------------------------------------
#  sette case


sette_umf = 0.022
sette_udiff=[0.27,0.39,0.72,0.81,0.97]
sette_uq = (np.array(sette_udiff) + sette_umf) / sette_umf

h02=[0.00559,0.0100,0.0145,0.0157,0.0201]
h03=[0.0145,0.0168,0.0179,0.0190,0.0280]
h04=[0.019,0.0235,0.0280,0.0335,0.0414]

sette_x = [sette_udiff, sette_uq]
sim_x = [0.72, (0.72 + 0.022)/0.022]
xlabels = ['$u-u_{mf}$ [m/s]','$u/u_{mf}$ [-]']

plt.figure(1)
plt.plot(sette_uq,h02,'o',label='Experiment')



symbols = ['s','x']
colors = ['g','r']
name=['D-fit','q-signal']

for i,method in enumerate([dvalues, d_qsignal]):
    plt.plot(sim_x[1],method[6],symbols[i],color=colors[i],alpha=0.5,label=f'{name[i]}')
    plt.plot(sim_x[1],method[8],symbols[i],color=colors[i],label=f'{name[i]} large')

plt.legend()
plt.grid(alpha=0.4)

plt.xlabel('$u/u_{mf}$ [-]')
plt.ylabel('$D_L$ $[m^2/s]$')
plt.title('Lighter bed - By method')
plt.savefig('paper/experimental_comparison/sette.pdf',bbox_inches='tight')




#%%
exp_x = [udiff_exp, uq_exp]
sim_x = [U_diff,Uq]

plt.figure(3,(4,4))

plt.plot(np.repeat(exp_x[1],2),nowallLP,'o',label='Experimental')
#plt.errorbar(np.repeat(exp_x[1],1),mean_exp,yerr=err_exp, fmt='o',label='Experimental')

plt.plot(np.repeat(sim_x[1],2),[dvalues[0],d_qsignal[0]],'s',label='Nozzles')
plt.plot(np.repeat(sim_x[1],4),[dvalues[5],dvalues[7],d_qsignal[5],d_qsignal[7]],'x',label='Porous plate')

plt.legend(loc='upper right')
plt.xlabel('$u/u_{mf}$ [-]')
plt.ylabel('$D_L$ $[m^2/s]$')
plt.title('Shallow bed - By boundary condition')
plt.grid(alpha=0.4)
plt.savefig('paper/experimental_comparison/with_porous.pdf',bbox_inches='tight')


#%%

exp_x = [udiff_exp, uq_exp]
sim_x = [U_diff,Uq]



plt.figure(4,(4,4))



plt.plot(np.repeat(exp_x[1],2),nowallLP,'o',label='Experimental')
#plt.errorbar(np.repeat(exp_x[1],1),mean_exp,yerr=err_exp, fmt='o',label='Experimental')

plt.plot(np.repeat(sim_x[1],1),[dvalues[0]],'s',label='Method 2')
plt.plot(np.repeat(sim_x[1],1),[d_qsignal[0]],'x',label='Method 3')
#plt.plot(np.repeat(sim_x[1],4),[dvalues[5],dvalues[7],d_qsignal[5],d_qsignal[7]],'x',label='Porous plate')

plt.legend(loc='upper left')
plt.xlabel('$u/u_{mf}$ [-]')
plt.ylabel('$D_L$ $[m^2/s]$')
plt.title('Shallow bed - Nozzles')
plt.grid(alpha=0.4)
plt.savefig('paper/experimental_comparison/nozzles.pdf',bbox_inches='tight')

#%%
# Experimental vs simulated
plt.figure(5)
exp = nowallLP.reshape(3,2)
limits=[]
for i in [0,1]:
    limits.append((exp[1][i]-exp[0][i])/(uq_exp[1]-uq_exp[0])*(Uq-uq_exp[0]) + exp[0][i])

mean=sum(limits)/2
err = np.abs(np.array(limits) - mean)

plt.errorbar(np.repeat(exp_x[1],1),mean_exp,yerr=err_exp, fmt='o',label='Experimental')
plt.errorbar(Uq,mean,yerr=err[0],fmt='o')

plt.plot(np.repeat(sim_x[1],1),[dvalues[0]],'s',label='D-fit')
plt.plot(np.repeat(sim_x[1],1),[d_qsignal[0]],'x',label='q-signal')



plt.figure(6)

deviation=0.1
x=np.linspace(0,1)

xhigh=x*(1 + deviation)
xlow=x*(1 - deviation)



plt.plot(x,x,'k--',alpha=0.3)
plt.plot(x,xhigh,'k--',alpha=0.3)
plt.plot(x,xlow,'k--',alpha=0.3)


labels= [['Nozzles','_nolabel_'],
         ['Porous','_nolabel_'],
         ['Porous large','_nolabel_'],
         ['Sette','_nolabel_'],
         ['Sette large','_nolabel_']]


for method,label in zip([dvalues,d_qsignal],labels[0]):
    plt.plot(np.repeat(mean,1),method[0],'o',color ='green', label=label)


for method,label in zip([dvalues,d_qsignal],labels[1]):
    plt.plot(np.repeat(mean,1),method[5],'o',color=r, label=label)


for method,label in zip([dvalues,d_qsignal],labels[2]):
    plt.plot(np.repeat(mean,1),method[7],'o',color=rd, label=label)


for method,label in zip([dvalues,d_qsignal],labels[3]):
    plt.plot(np.repeat(h02[2],1),method[6],'o',color=c,label=label)


for method,label in zip([dvalues,d_qsignal],labels[4]):
    plt.plot(np.repeat(h02[2],1),method[8],'o',color=cd,label=label)

#plt.plot([mean,mean,mean,h02[2],h02[2]],
#         [kd_values[0],kd_values[5],kd_values[7],kd_values[6],kd_values[8]],
#         'o',color='red',alpha=0.2,label='k-fit')

#plt.plot([mean],DT_method[7],'ko',label='DT')



plt.legend()
plt.xlim(0.002,0.016)
plt.ylim(0.0001,0.025)
plt.title('Experimental vs Simulated')
plt.xlabel('Experimental  $D_L$ $[m^2/s]$')
plt.ylabel('Simulated  $D_L$ $[m^2/s]$')
plt.savefig('paper/experimental_comparison/diagonal_reduced.pdf')













