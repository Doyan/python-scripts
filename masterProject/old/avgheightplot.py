# -*- coding: utf-8 -*-
"""
Created on Sun Jun  3 20:05:10 2018

@author: Gabriel
"""

import numpy as np, matplotlib.pyplot as plt
import pandas as pd
# -----------------------------------------------------------------------------
path = './datafiles/'

fnams =['30dp_avg_bed_height.csv','50dp_avg_bed_height.csv','70dp_avg_bed_height.csv']


data = []

for fnam in fnams:
    data.append(pd.read_csv(path+fnam,sep=',',skiprows=1,names=['time','havg'] ))

#%%

t = []
h = []

hcumm=[]


for frame in data:
    N = frame.shape[0]
    cummean = frame.havg.rolling(N).mean()    
    hcumm.append(cummean.iloc[-1])
    t.append(np.array(frame.time))
    h.append(np.array(frame.havg))




start=0
stop = len(h[0])
plt.figure(0)


plt.plot(t[0][start:stop],h[0][start:stop],'-',label='30dp, Mean = {:.3}'.format(hcumm[0]))
plt.plot([25,32],[0.63,0.63],'--',color='k', label='y = 0.63')

plt.plot(t[1][start:stop],h[1][start:stop],'-',label='50dp, Mean = {:.3}'.format(hcumm[1]))
plt.plot(t[2][start:stop],h[2][start:stop],'-',label='70dp, Mean = {:.3}'.format(hcumm[2]))


plt.xlim(25,32)
plt.xlabel('Time [s]')
plt.ylabel('Average bed height [m]')
plt.title('Height comparison')
plt.legend()
plt.grid()
plt.savefig('havg.pdf', bbox_inches='tight')





