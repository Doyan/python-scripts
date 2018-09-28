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
    
    tscaling = round(float(filelist[0].split('_')[1]),6)
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
    x1 = 0.54
    y0 = 0.0
    y1 = 0.4200001

    # Spec wall extent for masking 
    wallxmin = 0.255
    wallxmax = 0.285
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
    
    
    prefix = savepath + '{}_'.format(knumber)
    
    np.save(prefix + '2D-Temp',zM)
    np.save(prefix + '1D-Temp',barM)
    
    if include_static:
        np.save(savepath + '2D-grid',grid)
        np.save(savepath + 'x-coords',gx)
        np.save(savepath + 'y-coords',gy)
    
    # append tdata to its own registry-file
    tdatapath=savepath + 'tdata.csv'
    
    try: 
        tdf = pd.read_csv(tdatapath,index_col='k')
    except FileNotFoundError:
        tdf = pd.DataFrame(columns=['k','ts','freq'])
        tdf=tdf.set_index('k')
       
    tdf.loc[knumber] = {'freq': tfreq, 'ts': tscaling}
    tdf.to_csv(tdatapath)
    
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

# function to retrieve k-data in 1 or 2-D together with its time related data.
def loadK(knumber,dim='1d'):
    
    # choice of dimension
    if (dim == '1d') or (dim == '1D') or (dim < 2):
        filename = savepath + '{}_1D-Temp.npy'.format(knumber)
        gridname = savepath + 'x-coords.npy'
    elif (dim == '2d') or (dim == '2D') or (dim < 3):
        filename = savepath + '{}_2D-Temp.npy'.format(knumber)
        gridname = savepath + 'grid.npy'
    else:
        print('\ninvalid dimension')
        
    # load data    
    kdata = np.load(filename)
    grid = np.load(gridname)
    
    # retrieve time info
    tdata = pd.read_csv(savepath + 'tdata.csv',index_col='k')
    tfreq = tdata.loc[knumber].freq
    tscaling = tdata.loc[knumber].ts
    

    
    return kdata, tfreq, tscaling, grid



#----- Functions for handling MF Data ---------------------------------------        

# list multiphase folder and return sorted array of times and filenames
def sortMFolder(folderpath,tscaling,toffset):
    filelist=os.listdir(folderpath)

    tsteps = []
    for file in filelist:
        tsteps.append(file.split('_')[2].split('.')[0])

    times = [round(float(tstep)*tscaling + toffset,4) for tstep in tsteps]

    fileindex = list(sorted(zip(times,filelist)))
    return fileindex

# read a certain timestep file in a certain casefolder and return lists 
def fetchMdata(caseNo,sampleNo):

    folderpath = mfdatapath + 'case0' + str(caseNo) + '/'


    # value to scale integer from filename with
    tscaling = 1e-6

    # value to add to captured time to offset time before averaging process
    toffsets=[0,-3.6,-3.3,-3.3,-3.5]
    toffset = toffsets[caseNo]


    fileindex = sortMFolder(folderpath,tscaling,toffset)

    filename = fileindex[sampleNo][1]
    time = fileindex[sampleNo][0]
    df = pd.read_csv(folderpath+filename,names=['x','Ti','Tavg'])

    # unpack into arrays
    x = np.array(df.x)
    Ti = np.array(df.Ti)
    Tavg = np.array(df.Tavg)
    
    dcell = np.diff(x)[0]
    x0 = x[0]
    x = x - x0 + dcell/2 + dcell
    
    
    return x, Ti, Tavg, time
        
# -------------- Main loop ---------------------------------------------------
#%%



fileindex,_,_ = sortKfolder(runpath + 'k1000/')

x,y,T=loadTimeStep(fileindex,20)

M,l=getOnGrid(x,y,T)




caseNo=2
sampleNo = 229


knumber=2000


kdata,tfreq,tscaling,grid = loadK(knumber)


T0 = []
T1 = []
Sno = []

kT0 = []
kT1 = []

for sno in range(230):
    mx,mT,_,time = fetchMdata(caseNo,sno)
    T0.append(mT[0])
    T1.append(mT[-1])
    Sno.append(sno)
    ts = sno +1
    kT0.append(kdata[ts][1])
    kT1.append(kdata[ts][-2])

Sno=np.array(Sno)
kT0=np.array(kT0)
kT1=np.array(kT1)

knumber=2000


kdata,tfreq,tscaling,grid = loadK(knumber)

x0=grid[0]

mx, mT,_, time = fetchMdata(caseNo,sampleNo)

ts = sampleNo + 1




## ------------ Plotting -------------------------------------
#
def furbish():
    plt.plot([0.255,0.255],[1123.15,1073.15],'--',color=[0.8,0.8,0.8])
    plt.plot([0.285,0.285],[1123.15,1073.15],'--',color=[0.8,0.8,0.8])
    
    plt.plot([0.18,0.18],[1123.15,1073.15],'--',color=[0.8,0.8,0.8])
    plt.plot([0.36,0.36],[1123.15,1073.15],'--',color=[0.8,0.8,0.8])
    
    
    plt.xlabel('x-coordinate [m]')
    plt.title('Temperature gradient at t = ' + str(time) + 's')
    plt.xlim(0.0,0.54)
    return

kT=kdata[ts]

plt.plot(grid,kT)
plt.plot(mx,mT)
furbish()

plt.figure()
plt.plot(Sno,T0,Sno,T1)


dT=3
plt.plot(Sno+1,kT0-dT,Sno+1,kT1+dT)

plt.plot([Sno[0],Sno[-1]],[1123.15,1123.15],'r--')

plt.plot([Sno[0],Sno[-1]],[1073.15,1073.15],'b--')
plt.xlim(0,200)

#Low
#plt.ylim(1070,1090)

#high
#plt.ylim(1110,1124)


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


