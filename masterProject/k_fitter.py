# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 16:43:30 2018

k_fitter -

Creates 4 dimensional numpy representations of a series of 
fluent Temperature fields



@author: Gabriel
"""

import numpy as np, matplotlib.pyplot as plt

import pandas as pd
import os
from scipy.interpolate import griddata

# -------------------------------------------------------
datapath = './datafiles/k2000/'

# value to scale integer from filename with
tscaling = 0.05
toffset = 0

filelist = os.listdir(datapath)

times = [ float(file.split('Data-')[1]) for file in filelist]
times = np.array(times) * tscaling + toffset   


filename=filelist[230]
df = pd.read_csv(datapath+filename,names = ['cell','x','y','temp'],skiprows=[0])

# unpack into numpy array
y = np.array(df.y)
x = np.array(df.x)
T = np.array(df.temp)

x0 = 0.0
x1 = 0.51
y0 = 0.0
y1 = 0.4200001
dcell = 0.015

# create uniform grid axis points
gx = np.arange(x0+dcell/2,x1-dcell/2,dcell)
gy = np.arange(y0+dcell/2,y1-dcell/2,dcell)

# plot fluent grid for inspection
plt.plot(x,y,'b.')

plt.plot(gx,np.ones_like(gx)*gy[17],'r.')
plt.plot(np.ones_like(gy)*gx[15],gy,'y.')
plt.plot(np.ones_like(gy)*gx[18],gy,'y.')
plt.show()

# make meshgrid from grid axes
xi, yi = np.meshgrid(gx,gy)

# interpolate fluent grid onto structured python grid

z0 = griddata((x,y),T,(xi,yi),method='nearest')

# create mask to remove values "inside wall"
mask = (xi > 0.24) & (xi < 0.27) & (yi > 0.42 - 0.15)

z0[mask] = np.nan

# plot interpolated grid as surface
plt.contourf(xi,yi,z0,10,cmap='bwr', vmin=1073.15, vmax=1123.15)

plt.colorbar()

plt.contour(xi,yi,z0,36,colors='black', vmin=1073.15, vmax=1123.15)

plt.plot(x,y,'b.')
plt.xlim(0.0,0.51)
plt.ylim(0.0,0.42)

xbar = []
for i in range(len(gx)):
    xsum = 0
    n = 0
    for j in range(len(gy)):
        if not np.isnan(z0[j][i]):
            xsum += z0[j][i]
            n = n + 1
    xbar.append(xsum / n) 
plt.show() 
    
plt.plot(gx,xbar)
plt.plot([0,0.51],[1123.15,1073.15],'k--')
plt.plot([0.24,0.24],[1123.15,1073.15],'--',color=[0.8,0.8,0.8])
plt.plot([0.27,0.27],[1123.15,1073.15],'--',color=[0.8,0.8,0.8])
#plt.grid(linestyle='-.',color=[0.9,0.9,0.9])



