#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 21:16:05 2018

@author: gabgus
"""

import numpy as np, pandas as pd
from scipy.interpolate import griddata

import os


# -------------------------------------------------------------------------
# Case info / metadata

savepath='/scratch/gabgus/mfdata_all/meshed/'


readSpecs={'mfolder': '/scratch/gabgus/mfdata_all/',
           'delims': [0,1,1,1,1,0,1,1,1,1,1], 
           'ending': ['error','.csv','.csv','.csv','.csv','missing','.csv','.csv','.csv','.csv','.csv']}

readSpecs['foldername'] = ['error', 
         readSpecs['mfolder'] + 'dysor/case01/',
         readSpecs['mfolder'] + 'dysor/case02/',
         readSpecs['mfolder'] + 'dysor/case03/', 
         readSpecs['mfolder'] + 'dysor/case04/',
         'missing',
         readSpecs['mfolder'] + 'pp_real/noWall/', 
         readSpecs['mfolder'] + 'pp_real/wWall/',
         readSpecs['mfolder'] + 'sette/',
         readSpecs['mfolder'] + 'larger_domains/noWall/',
         readSpecs['mfolder'] + 'larger_domains/sette/']
           
Wlimits =  {'w0': ['error' , np.nan, 0.4875, 0.4125, 0.4575, 0.4875, np.nan, 0.4875, np.nan,np.nan,np.nan],
           'w1': ['error' , np.nan, 0.5175, 0.4425, 0.4875, 0.5175, np.nan, 0.5175, np.nan,np.nan,np.nan],
           'wy': 0.3,
           'hasWall': [False, False, True, True, True, True, False, True, False, False, False]}

tInfo = {'toffset': ['error', 3.6, 3.3 , 3.3, 3.5, 'missing', 5.0, 5.0, 4.5, 5.45, 2.75],
         'Tstart': ['error',0,0,0,0,'missing',50,50,29,108,55],
         'noSamples' : ['error',193,191,191,188,'missing',787,756,196,520,207],
         'sstep': ['error',1,1,1,1,'missing',4,4,1,1,1]}


# -------------------------------------------------------------------------
# Functions

def sortCsvByTime(caseNo,delimiter='_'):
    datafolder = readSpecs['foldername'][caseNo]
    ndelimiters = readSpecs['delims'][caseNo]
    ending = readSpecs['ending'][caseNo]
    
    filelist = [datafolder + file for file in os.listdir(datafolder) if file.endswith(ending)]
    times = [float(filename.split(delimiter)[ndelimiters].split(ending)[0]) for filename in os.listdir(datafolder)]
    return list(sorted(zip(times,filelist)))



def loadTimestep(caseNo,filelist,sampleNo,quantity='temp'):
    
    filename = filelist[sampleNo][1]
    time = filelist[sampleNo][0]
    time = round(time - tInfo['toffset'][caseNo],4) 
    
    data=pd.read_csv(filename)
    data.columns=[name.strip() for name in data.columns]

    x=np.array(data['x-coordinate'])
    y=np.array(data['y-coordinate'])
    z=np.array(data['z-coordinate'])
    
    sdict={}
    
    sdict['temp'] = np.array(data['phase-2-temperature'])
    sdict['vof']  = np.array(data['phase-2-vof'])
    sdict['uds']  = np.array(data['phase-2-uds-0-scalar'])
    sdict['xvel'] = np.array(data['phase-2-x-velocity'])
    sdict['yvel'] = np.array(data['phase-2-y-velocity'])
    sdict['zvel'] = np.array(data['phase-2-z-velocity'])
    
    return x,y,z,sdict[quantity],time

def makeMesh(x,y,z):
    
    xpoints = np.unique(np.round(x,5))
    #ypoints = np.unique(np.round(y,3))
    zpoints = np.unique(np.round(z,5))

    dx =np.diff(xpoints)
    #dy =np.diff(ypoints)
    dz =np.diff(zpoints)

    # only consider x and z because fluent makes funny stuff with the ys
    dcell = np.concatenate((dx,dz)).mean()


    x0=round(np.min(x) - dcell/2,4)
    x1=round(np.min(x) + np.max(x) + dcell,4)

    y0=round(np.min(y) - dcell/2,4)
    y1=round(np.min(y) + np.max(y) + dcell,4)

    z0=round(np.min(z) - dcell/2,4)
    z1=round(np.min(z) + np.max(z) + dcell,4)

    gx = np.arange(x0+dcell/2,x1-dcell/2,dcell)
    gy = np.arange(y0+dcell/2,y1-dcell/2,dcell)
    gz = np.arange(z0+dcell/2,z1-dcell/2,dcell)

    xi, yi, zi = np.meshgrid(gx,gy,gz)
    return (xi,yi,zi),(gx,gy,gz), dcell
    


def packScalar(caseNo,scalar):
    filelist = sortCsvByTime(caseNo)

    x,y,z,s,time = loadTimestep(caseNo,filelist,0,quantity=scalar)

    M,g,dcell = makeMesh(x,y,z)

    Smat = []
    t = []

    Tstart=tInfo['Tstart'][caseNo]
    noSamples = tInfo['noSamples'][caseNo]
    sstep=tInfo['sstep'][caseNo]

    for i in range(Tstart,noSamples,sstep):
        print('{:.3f}'.format((float(i) - Tstart) / (noSamples-Tstart)))
        x,y,z,s,time = loadTimestep(caseNo,filelist,i,quantity=scalar)

        Si =  griddata((x,y,z),s,M,method='nearest')

        xi = M[0]
        yi = M[1]
        
        # mask out the cells "inside the wall"
        if Wlimits['hasWall'][caseNo]:
            mask = (xi > Wlimits['w0'][caseNo]) & (xi < Wlimits['w1'][caseNo]) & (yi > Wlimits['wy'])
            Si[mask] = np.nan
    
        Smat.append(Si)
        t.append(time)

    Smat = np.array(Smat)
    t = np.array(t)
    g = np.array(g)
    M = np.array(M)
    return Smat,M,g,t,dcell


def saveScalar(caseNo,scalar):
    Smat,M,g,t,dcell=packScalar(caseNo,scalar)
    Sname = savepath + 'c{}_{}'.format(caseNo,scalar)
    meshname = savepath + 'c{}_mesh'.format(caseNo)
    gridname = savepath + 'c{}_gridpoints'.format(caseNo)
    timename = savepath + 'c{}_times'.format(caseNo)
    
    np.save(Sname,Smat)
    np.save(meshname,M)
    np.save(gridname,g)
    np.save(timename,t)
    return 

def loadMdata(caseNo,scalar='temp'):
    mfdatapath=savepath
    if scalar == 'grid':
        M = np.load('{}c{}_mesh.npy'.format(mfdatapath,caseNo))
        g = np.load('{}c{}_gridpoints.npy'.format(mfdatapath,caseNo))
        t = np.load('{}c{}_times.npy'.format(mfdatapath,caseNo))
        return M,g,t
        
    Smat = np.load('{}c{}_{}.npy'.format(mfdatapath,caseNo,scalar))
    M = np.load('{}c{}_mesh.npy'.format(mfdatapath,caseNo))
    g = np.load('{}c{}_gridpoints.npy'.format(mfdatapath,caseNo))
    t = np.load('{}c{}_times.npy'.format(mfdatapath,caseNo))
    return Smat,M,g,t
    

def packQsignal(fileindex,dcell,prefix):
    Q1=[]
    Q2=[]
    t=[]
    for i,(time,file) in enumerate(fileindex):
        data=pd.read_csv(file)
        data.columns=[name.strip() for name in data.columns]
        Q1.append(data['udm-1'].sum()*(dcell**3))
        Q2.append(data['udm-2'].sum()*(dcell**3))
        t.append(time)

    np.save('/scratch/gabgus/mfdata_all/meshed/{}_qin'.format(prefix),np.array(Q1))
    np.save('/scratch/gabgus/mfdata_all/meshed/{}_qut'.format(prefix),np.array(Q2))
    np.save('/scratch/gabgus/mfdata_all/meshed/{}_qtime'.format(prefix),np.array(t))
    return


# -------------------------------------------------------------------------
    
#
for caseNo in [9,10]:
    for scalar in ['temp','vof','uds','xvel','yvel']:
        saveScalar(caseNo,scalar)
    
    filelist = sortCsvByTime(caseNo)
    x,y,z,s,time = loadTimestep(caseNo,filelist,0,quantity=scalar)

    M,g,dcell = makeMesh(x,y,z)
    
    prefix = 'c{}'.format(caseNo)
    packQsignal(filelist,dcell,prefix)










# Smat,M,g,t,dcell=packScalar(7,'temp')


#plt.plot([gx[xslice.start],gx[xslice.stop - 1]],[1123.15,1073.15],'--k')
#
##dflist=[]    
##   
##csvlist = [file for file in os.listdir(exportpath) if file.endswith('.csv')]
##
##for csv in csvlist:
##    colnames=['cno','x','y','z','uds','xvel','yvel','zvel','vof','T']
##    df=pd.read_csv(exportpath + csv,skiprows=[0],names=colnames)
##    dflist.append(df)
##    
##
##
##dflist[0].uds
