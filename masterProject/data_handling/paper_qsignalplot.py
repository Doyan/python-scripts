#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  6 13:11:20 2019

@author: gabgus
"""

import numpy as np
import matplotlib.pyplot as plt
import comparator as C
import scipy.stats as st
from scipy.stats import skewtest, skew

caseNo=10
mfdatapath=C.mfdatapath


qin = np.load('{}c{}_qin.npy'.format(mfdatapath,caseNo))
qut = np.load('{}c{}_qut.npy'.format(mfdatapath,caseNo))
qtime = np.load('{}c{}_qtime.npy'.format(mfdatapath,caseNo))
qmid =(qin+qut)/2

A=[0,0.6,0.6,0.6,0.6,0,0.6,0.6,0.2,0.6,0.2]
dx = [0,0.525,0.525,0.525,0.525,0,0.525,0.525,0.306,1.50,0.72]
dT = 50.0

vof_s=C.vof_s

cfactor=vof_s[caseNo]*1200*2600


start_idx = (qmid!=0).argmax()




sslice=slice(start_idx,int(1.0*len(qtime)))

q = qmid[sslice]
qtime=qtime[sslice]

qx = q/A[caseNo]    
keff = qx * dx[caseNo] / dT
deff= keff/cfactor


q_int = st.t.interval(0.999, len(qx)-1, loc=qx.mean(), scale=st.sem(qx))
k_int = st.t.interval(0.999, len(keff)-1, loc=keff.mean(), scale=st.sem(keff))
d_int=np.array(k_int)/cfactor

N = int(0.1*len(qx))

qroll= C.moving_average(qx,N)

keff_roll = qroll * dx[caseNo] / dT
d_roll=keff_roll/cfactor


plt.plot(qtime,qx/1e6)
plt.plot(qtime[N-1:],qroll/1e6)
plt.grid()
#plt.ylim(0.2,2)
#plt.xlim(10,15)
print('qx = {} MW/m2'.format(qx.mean()/1e6))
plt.show()

plt.plot(qtime,keff,label='raw')
plt.plot(qtime[N-1:],keff_roll,label='rolling mean')
plt.plot([qtime[0],qtime[-1]],[keff.mean(),keff.mean()],'--',label='mean')
plt.xlabel('Time [s]')
plt.ylabel('$k_{eff}$ [W/m/K]')
plt.title('Effective conduction over time')
plt.grid()
plt.xlim(qtime[0],qtime[-1]+0.5)
plt.show()


plt.plot(qtime,deff,label='raw')
plt.plot(qtime[N-1:],d_roll,label='rolling mean')
plt.plot([qtime[0],qtime[-1]],[deff.mean(),deff.mean()],'--',label=f'mean = {deff.mean():.2}')
plt.xlabel('Time [s]')
plt.ylabel('$D_L$ [$m^2/s$]')
plt.title('Lateral dispersion, calculated from heat flux')
plt.grid()
plt.xlim(qtime[0],qtime[-1])
plt.legend(loc='upper left') 
#plt.savefig('paper/qsignal_c9_small.pdf',bbox_inches='tight')
plt.show()








