# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 16:43:30 2018

k_fitter -

Compares solutions of a heat transfer problem with a certain k to 
a corresponding temperature field from a multiphase simulation 



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
df = pd.read_csv(datapath+filename,names = ['x','y','x1','y1','temp'],skiprows=[0])

# unpack into numpy array
y = np.array(df.y)
x = np.array(df.x)
T = np.array(df.temp)

# create uniform grid axis points
gx = np.arange(0.0,0.51,0.015)
gy = np.arange(0.0,0.42,0.015)

# plot fluent grid for inspection
plt.plot(x,y,'b.')

plt.plot(gx,np.ones_like(gx)*gy[17],'r.')
plt.plot(np.ones_like(gy)*gx[14],gy,'y.')
plt.plot(np.ones_like(gy)*gx[19],gy,'y.')
plt.show()

# make meshgrid from grid axes
xi, yi = np.meshgrid(gx,gy)

# interpolate fluent grid onto structured python grid

z0 = griddata((x,y),T,(xi,yi),method='nearest')
z1 = griddata((x,y),T,(xi,yi),method='linear')

#nanmask = np.isnan(z1)
#z1[nanmask] = z0[nanmask]    



# create mask to remove values "inside wall"
mask = (xi > 0.24) & (xi < 0.27) & (yi > 0.42 - 0.15)

z1[mask] = np.nan

# plot interpolated grid as surface
plt.contourf(xi,yi,z0,20,cmap='bwr', vmin=1073.15, vmax=1123.15)

plt.colorbar()
#
#plt.contour(xi,yi,z1,20,colors='black', vmin=1073.15, vmax=1123.15)
#
#plt.plot(x,y,'b.')


#xbar = []
#for i in range(len(gx)):
#    xsum = 0
#    n = 0
#    for j in range(len(gy)):
#        if not np.isnan(zi[i][j]):
#            xsum += zi[gx[i]][gy[j]]
#            n = n + 1
#    xbar.append(xsum / n) 
    







