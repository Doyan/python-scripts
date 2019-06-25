#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  8 10:31:55 2019

@author: gabgus
"""

import numpy as np
import matplotlib.pyplot as plt




dpath = 'DTfolder/'
fracs = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4]


planes = []
for elevation in ['high','low']:
    plane = []
    for direction in ['x','y','z']:
        DTmax=[]
        DTmin=[]
        DTmean=[]
        for frac in fracs:    
            DTMat = np.load(dpath+'DTMat_{}_{}_{}_26ish_fft.npy'.format(elevation,direction,frac))
            
            DTmean.append(DTMat.mean(axis=1).mean(axis=0))
            DTmax.append(DTMat.mean(axis=1).max(axis=0))
            DTmin.append(DTMat.mean(axis=1).min(axis=0))
        
        plt.plot(fracs,DTmax,label='max')
        plt.plot(fracs,DTmean,label='mean' )
        plt.plot(fracs,DTmin,label='min')
        
        plt.ylabel('$D_T \;\; [m^2/s]$')
        plt.xlabel('fraction of correlated time series')    
        plt.title('Turbulent dispersion in {}-direction'.format(direction))
        plt.legend()
        plt.grid()
        
        plt.savefig('paper/DT_{}_plane_averaged_{}.pdf'.format(elevation,direction))
        plt.show()
        plane.append(DTmean[2])
    
    planes.append(plane)
    

#%%

planes = np.array(planes)
 
n_groups=3    

fig, ax = plt.subplots()

index = np.arange(n_groups)
bar_width = 0.17

opacity = 0.4
error_config = {'ecolor': '0.3'}


ticklabels=['x','y','z']

ax.zorder = 0
ax.grid(linestyle = '--', alpha=0.5,zorder=0,axis='y')


high= ax.bar(index+ bar_width*0,planes[0] , bar_width,
                alpha=opacity, color='b',
                label='y = 0.3525 m',zorder=40)


low = ax.bar(index + bar_width*1, planes[1], bar_width,
                alpha=opacity, color='g',
                label='y = 0.1125 m ',zorder=41)

ax.set_xticks(index + bar_width/2)
ax.set_xticklabels(ticklabels)
ax.set_yticks([i for i in np.arange(0,0.055,0.005)])
ax.legend()
ax.set_ylabel('$D_T \;\; [m^2/s]$')

fig.suptitle('Turbulent dispersion along the x,y and z directions')
fig.savefig('paper/DT_bar_chart.pdf')





