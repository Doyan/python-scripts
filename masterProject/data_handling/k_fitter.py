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
#datapath = './datafiles/k2000/'
runpath = '/scratch/gabgus/fluent/kgenerator/runs/' 

# parse through and sort contents of folder according to timestep 
# for easy access by later functions
def sortKfolder(runpath,knumber,toffset=0):
    datapath = runpath + 'k' + str(knumber) + '/'

    filelist = os.listdir(datapath)
    
    tscaling = float(filelist[0].split('_')[1])
    times = [ round(float(file.split('ts-')[1]) * tscaling + toffset,4) for file in filelist]
    files = [datapath + file for file in filelist]
    

    #times = np.array(times) * tscaling + toffset   

    fileindex=list(sorted(zip(times,files)))
    tfreq = float(fileindex[1][1].split('ts-')[1]) - float(fileindex[0][1].split('ts-')[1])
    return fileindex, tscaling, tfreq

# load fluent data from a file given timestep
def loadTimeStep(fileindex,stepnum):
    # read corresponding datafile into dataframe
    filename=fileindex[stepnum][1]
    df = pd.read_csv(filename,names = ['cell','x','y','temp'],skiprows=[0])

    # unpack into numpy array
    y = np.array(df.y)
    x = np.array(df.x)
    T = np.array(df.temp)
    return x,y,T


# interpolate unordered fluent data onto structured callable grid
# inserts nan where grid extends onto area not covered by fluent
def getOnGrid(x,y,T):
    # Spec boundary of grid
    x0 = 0.0
    x1 = 0.51
    y0 = 0.0
    y1 = 0.4200001

    # Spec wall extent for masking 
    wallxmin = 0.24
    wallxmax = 0.27
    wally = 0.15

    # Spec cell width
    dcell = 0.015

    # create uniform grid axis points
    gx = np.arange(x0+dcell/2,x1-dcell/2,dcell)
    gy = np.arange(y0+dcell/2,y1-dcell/2,dcell)


    # make meshgrid from grid axes
    xi, yi = np.meshgrid(gx,gy)# plot fluent grid for inspection

    # interpolate fluent grid onto structured python grid

    zi = griddata((x,y),T,(xi,yi),method='nearest')

    # create mask to remove values "inside wall"
    mask = (xi > wallxmin) & (xi < wallxmax) & (yi > y1 - wally)

    zi[mask] = np.nan
    return (xi,yi,zi), (gx, gy)

# Make average along y-coordinate
def yAverage(z0):
   Nj = z0.shape[0]
   Ni = z0.shape[1]
    
   xbar = []
   for i in range(Ni):
       xsum = 0
       n = 0
       for j in range(Nj):
           if not np.isnan(z0[j][i]):
               xsum += z0[j][i]
               n = n + 1
               
       xbar.append(xsum / n)
   return xbar

# Interpolate all fluent grids in folder onto ordered grid and stack into 
# array with order Temp[time][x][y] 
def packKfolder(runpath,knumber):    
    fileindex,tscaling,tfreq = sortKfolder(runpath,knumber)
    
    zlist = []
    zbar = []
    for sampleNo in range(len(fileindex)):
        x,y,T = loadTimeStep(fileindex,sampleNo)
        
        M, l = getOnGrid(x,y,T)
        (xi,yi,z0) = M
        
        xbar = yAverage(z0)
        zlist.append(z0)
        zbar.append(xbar)

    zM = np.dstack(zlist).T
    barM = np.squeeze(np.dstack(zbar).T)
    return zM, barM, (xi,yi), l, (tscaling,tfreq)

# -------------- Main loop ---------------------------------------------------
#%%
knumber = 1000
savepath= ./datafiles/kfiles


zM, barM, grid, l, tdata = packKfolder(runpath,knumber)

(xi,yi)=grid
(gx,gy)=l     
(tscaling, tfreq) = tdata

def saveData(savepath, knumber):
    return


time = 12.5
sample = int(time / tscaling / tfreq)   

z0=zM[sample].T
xbar=barM[sample]
# ------------ Plotting -------------------------------------

def furbish():
    plt.xlabel('x-coordinate [m]')
    plt.title('Temperature gradient at t = ' + str(time) + 's')
    plt.xlim(0.0,0.51)
    return

# plot interpolated grid as surface
plt.contourf(xi,yi,z0,10,cmap='bwr', vmin=1073.15, vmax=1123.15)
plt.colorbar()

plt.contour(xi,yi,z0,33,colors='black', vmin=1073.15, vmax=1123.15)

plt.ylabel('y-coordinate [m]')
furbish()
plt.ylim(0.0,0.42)
plt.show() 
    
plt.plot(gx,xbar)
plt.plot([0,0.51],[1123.15,1073.15],'k--')
plt.plot([0.24,0.24],[1123.15,1073.15],'--',color=[0.8,0.8,0.8])
plt.plot([0.27,0.27],[1123.15,1073.15],'--',color=[0.8,0.8,0.8])

plt.ylabel('Temperature [K]')
furbish()
plt.show()
#plt.grid(linestyle='-.',color=[0.9,0.9,0.9])


# plot fluent grid for inspection
#plt.plot(x,y,'b.')
#
#plt.plot(gx,np.ones_like(gx)*gy[17],'r.')
#plt.plot(np.ones_like(gy)*gx[15],gy,'y.')
#plt.plot(np.ones_like(gy)*gx[18],gy,'y.')
#
#plt.xlim(0.0,0.51)
#plt.ylim(0.0,0.42)
#plt.show()


