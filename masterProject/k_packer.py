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
import subprocess
# -------------------------------------------------------
# Paths

# path to save structured collection of data
savepath= './datafiles/kfiles/'

# path to folder with saved fluent runs
runpath = '/scratch/gabgus/fluent/kgenerator/runs/' 

# path to fluent runscript
scriptname = '/scratch/gabgus/fluent/kgenerator/krun.sh'

# path to multiphase simulation data
mfdatapath = './datafiles/mfdata/'

# -------------------------------------------------------------------------
# Functions


# convenience function to generate foldername string
def genfonam(knumber):
    return runpath + 'k' + str(knumber) + '/'

# convenience function to generate knumber from foldername string
def genknumber(folderpath):
    return os.path.basename(os.path.dirname(folderpath)).split('k')[1]

# parse through and sort contents of folder according to timestep 
# for easy access by later functions
def sortKfolder(folderpath,toffset=0):
    
    filelist = os.listdir(folderpath)
    
    tscaling = float(filelist[0].split('_')[1])
    times = [ round(float(file.split('ts-')[1]) * tscaling + toffset,4) for file in filelist]
    files = [folderpath + file for file in filelist]
    

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
def packKfolder(folderpath):    
    fileindex,tscaling,tfreq = sortKfolder(folderpath)
    
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


# shoehorn all relevant data into numpy-arrays and save using .npy format
# files are prefixed with k and time info.    
def saveData(savepath,folderpath,include_static=False):
    zM, barM, grid, l, tdata = packKfolder(folderpath)
    (tscaling, tfreq) = tdata
    (gx,gy) = l
    grid = np.dstack(grid).T
    
    knumber=os.path.basename(os.path.dirname(folderpath)).split('k')[1]
    
    
    prefix = savepath + '{}_{:.2e}_{}_'.format(knumber,tscaling,int(tfreq))
    
    np.save(prefix + '2D-Temp',zM)
    np.save(prefix + '1D-Temp',barM)
    
    if include_static:
        np.save(prefix + '2D-grid',grid)
        np.save(prefix + 'x-coords',gx)
        np.save(prefix + 'y-coords',gy)
    
    return

# request for new runfolder to be made by fluent 
def requestK(knumber,scriptname):
    current= os.getcwd()
    scriptdir=os.path.dirname(scriptname)
    
    os.chdir(scriptdir)
    
    c=subprocess.run([scriptname, str(knumber)])
    
    os.chdir(current)
    
    return c.returncode, c.args

# check what k's we have already parsed through, 
# returns set which can be checked against
def checkDone(savepath):
    savelist = os.listdir(savepath)
    
    done = set()
    for saved in savelist:
        k = saved.split('_')[0]
        if not k in done:
            done.add(k)
    return done

# Function to make sure we have parsed all runs so far
def saveAll(runpath,savepath):
    runlist = os.listdir(runpath)

    done = checkDone(savepath)

    for runfolder in runlist:
        k = runfolder.split('k')[1]
        if not k in done:
            print('processed k:' + k) 
            folderpath = runpath+runfolder + '/'
            saveData(savepath,folderpath)
    return checkDone(savepath)

# Request k to be added to parsed collection. Runs fluent if needed.
def addK(knumber):
    done = checkDone(savepath)
    if str(knumber) in done:
        print('k = ' + str(knumber) + ' already exists in savepath')
        return
    else: 
        print('Making fluent run for k = ' + str(knumber) + ' W/mK \n')
        print('please wait...')
        ecode,args = requestK(knumber,scriptname)
        folderpath= genfonam(knumber)
        saveData(savepath,folderpath)
        if ecode == 0:
            print('\nSuccesfully added k = ' + str(knumber) + ' W/mK')
            return
        else:
            print('\nk-generation failed')
            return 1
        
        
# -------------- Main loop ---------------------------------------------------
#%%

caseNo=2

# value to scale integer from filename with
tscaling = 1e-6

# value to add to captured time to offset time before averaging process
toffsets=[0,-3.6,-3.3,-3.3,-3.5]

toffset = toffsets[caseNo]





# read and preprocess text file to match 0 - L scale 
def getDf(sampleNo,datalist):
    filepath = datapath + datalist[sampleNo][1]

    df = pd.read_csv(filepath,names=['x','Ti','Tavg'])
    dcell = df.x.diff().shift(-1)[0]
    x0 = df.x[0]
    df.x = df.x + dcell/2 - x0
    L = df.x[33] + dcell/2 - (df.x[0] - dcell/2)
    return df, L, x0, dcell





#zM, barM, grid, l, tdata = packKfolder(folderpath)

#(xi,yi)=grid
#(gx,gy)=l     
#(tscaling, tfreq) = tdata



#
#time = 12.5
#sample = int(time / tscaling / tfreq)   
#
#z0=zM[sample].T
#xbar=barM[sample]
## ------------ Plotting -------------------------------------
#
#def furbish():
#    plt.xlabel('x-coordinate [m]')
#    plt.title('Temperature gradient at t = ' + str(time) + 's')
#    plt.xlim(0.0,0.51)
#    return
#
## plot interpolated grid as surface
#plt.contourf(xi,yi,z0,10,cmap='bwr', vmin=1073.15, vmax=1123.15)
#plt.colorbar()
#
#plt.contour(xi,yi,z0,33,colors='black', vmin=1073.15, vmax=1123.15)
#
#plt.ylabel('y-coordinate [m]')
#furbish()
#plt.ylim(0.0,0.42)
#plt.show() 
#    
#plt.plot(gx,xbar)
#plt.plot([0,0.51],[1123.15,1073.15],'k--')
#plt.plot([0.24,0.24],[1123.15,1073.15],'--',color=[0.8,0.8,0.8])
#plt.plot([0.27,0.27],[1123.15,1073.15],'--',color=[0.8,0.8,0.8])
#
#plt.ylabel('Temperature [K]')
#furbish()
#plt.show()
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


