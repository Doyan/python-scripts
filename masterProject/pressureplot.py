# -*- coding: utf-8 -*-
"""
Created on Sun Jun  3 19:12:01 2018

@author: Gabriel
"""

import numpy as np, matplotlib.pyplot as plt
import pandas as pd
# -----------------------------------------------------------------------------
path = './datafiles/'

fnams =['30dp_pressure.csv','50dp_pressure.csv','70dp_pressure.csv']


data = []

for fnam in fnams:
    data.append(pd.read_csv(path+fnam,sep=',',skiprows=1,names=['time','pressure'] ))

#%%

t = []
p = []

pcumm=[]

N = 3500
for frame in data:
    t.append(np.array(frame.time))
    p.append(np.array(frame.pressure))
    cummean = frame.pressure.rolling(N).mean()    
    pcumm.append(cummean.iloc[-1])

start = 6000
stop = len(p[0])

plt.figure(0)
plt.yticks([-14000,-13000, -12000, -11000,-10000,-9000,-8000,-7000,-6000 ])
plt.plot(t[0][start:stop],p[0][start:stop],'-',label='30dp, Mean = {:.3}'.format(pcumm[0]))
plt.plot(t[1][start:stop],p[1][start:stop],'-',label='50dp, Mean = {:.3}'.format(pcumm[1]))
plt.plot(t[2][start:stop],p[2][start:stop],'-',label='70dp, Mean = {:.3}'.format(pcumm[2]))
plt.xlim(24,32)
plt.xlabel('Time [s]')
plt.ylabel('$\Delta P$ [Pa]')
plt.title('Pressure drop comparison')
plt.legend()
plt.savefig('pdrop_zoom.pdf', bbox_inches='tight')


start=0
plt.figure(1)
plt.plot(t[0][start:stop],p[0][start:stop],'-',label='30dp')
plt.plot(t[1][start:stop],p[1][start:stop],'-',label='50dp')
plt.plot(t[2][start:stop],p[2][start:stop],'-',label='70dp')
plt.xlim(0,32)
plt.xlabel('Time [s]')
plt.ylabel('$\Delta P$ [Pa]')
plt.title('Pressure drop comparison')
plt.savefig('pdrop_full.pdf', bbox_inches='tight')
