#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 28 11:15:27 2018

@author: gabgus
"""


import numpy as np, matplotlib.pyplot as plt
import pandas as pd

# -------------------------------------------------------------------------

fnam = 'Qbed_monitor.csv' 

data = pd.read_csv(fnam)

t = data.Time
q = data.Q_bed
#%%
N = len(t)

window = N

Qavg = []
qdiff = []
qdiff.append(0.0)
for i in np.arange(N):
    
    if i <= window:
        norm = i
        qsum = q[0:i].sum()
    else: 
        norm = window
        start = i - window
        qsum = q[start:i].sum()
    
    Qavg.append(qsum / norm)
    

    if i > 0:
        qdiff.append((Qavg[i] - Qavg[i-1]) / (20 * 1800))
    
    
Qsum = []
for i in np.arange(N):
    qsum=q[0:i].sum()
    Qsum.append(qsum) 
    





TH = 860    
TC = 800
    
k_avg = np.array(Qavg) * 1.4 / (0.3*1.2 * (TH - TC))

q100 = 100 * 1800 * 0.3 * 1.2 * (TH-TC) / 1.4


 
def q20(t,q0):
    k = 20 * 1800 * 0.3 * 1.2 * (TH-TC) / 1.4
    return k*t + q0


def q100(t,q0):
    k = 100 * 1800 * 0.3 * 1.2 * (TH-TC) / 1.4
    return k*t + q0




#%%

plt.figure(1)    
plt.plot(data.Time,data.Q_bed)


plt.figure(2)
plt.plot(t,Qavg)
plt.ylim(-1.5e7,0.5e7)
plt.xlim(8, 24)

plt.figure(3)
plt.plot(t,k_avg/1800)
plt.plot([0,19],[20,20],'--')
plt.xlim(8,24)
plt.ylim(-500,200)



plt.figure(4)
plt.plot(t,Qsum)
plt.plot(t,q100(t,0))



plt.figure(5)
plt.plot(t,qdiff)
plt.ylim(-2,2)
