#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 15:15:21 2018

@author: gabgus
"""

import numpy as np, matplotlib.pyplot as plt

import pandas as pd

from scipy.interpolate import griddata

from datetime import datetime

import os
import shutil
import tarfile

# -------------------------------------------------------
# Paths
os.chdir('/chalmers/users/gabgus/python-scripts/masterProject')

# path to save structured collection of data
savepath= '/scratch/gabgus/fluent/dgenerator/parsed/'

# path to folder with saved fluent runs
runpath = '/scratch/gabgus/fluent/dgenerator/runs/' 

# path to multiphase simulation data
mfdatapath='/scratch/gabgus/mfdata_all/meshed/'

# -------------------------------------------------------------------------
# Global variables

wLimits =  {'w0': ['error' , np.nan, 0.255+0.225, 0.180+0.225, 0.225+0.225, 0.255+0.225, np.nan, 0.255+0.225, np.nan],
            'w1': ['error' , np.nan, 0.285+0.225, 0.210+0.225, 0.255+0.225, 0.285+0.225, np.nan, 0.285+0.225, np.nan],
            'wy': 0.3,
            'hasWall': [False, False, True, True, True, True, False, True, False]}

# ----------------------------------------------------------------------------
# Functions

# load multiphase mesh
def loadMmesh(caseNo):
    M = np.load('{}c{}_mesh.npy'.format(mfdatapath,caseNo))
    g = np.load('{}c{}_gridpoints.npy'.format(mfdatapath,caseNo))
    t = np.load('{}c{}_times.npy'.format(mfdatapath,caseNo))
    return M,g,t

# get slices to apply on multiphase mesh make it correspond to k-domain
def getSlice(caseNo):
    Dlimits = {'x0': ['error' , 0, 0, 0, 0, 'missing', 0, 0, 0], 
           'x1': ['error' , 67, 67, 67, 67, 'missing', 67, 67, 100],
           'y1': ['error' , 38, 38, 38, 38, 'missing', 32, 32, 60],
           'y0': ['error' , 6, 6, 6, 6, 'missing', 0, 0, 0]}
    
    xslice=slice(Dlimits['x0'][caseNo],Dlimits['x1'][caseNo])  #max66
    yslice=slice(Dlimits['y0'][caseNo],Dlimits['y1'][caseNo]) #max67
    zslice=1            #max86
    return xslice, yslice, zslice

# parse through and sort contents of folder according to timestep 
# for easy access by later functions
def sortDfolder(folderpath):
    filelist = [folderpath + file for file in os.listdir(folderpath)]
    times = [float(filename.split('-')[1]) for filename in os.listdir(folderpath)]
    return list(sorted(zip(times,filelist)))


# load fluent data from a file given timestep
def loadTimeStep(fileindex,stepnum):
    # read corresponding datafile into dataframe
    filename=fileindex[stepnum][1]
    df = pd.read_csv(filename,names = ['cell','x','y','conc'],skiprows=[0])

    # unpack into numpy array
    y = np.array(df.y)
    x = np.array(df.x)
    C = np.array(df.conc)
    return x,y,C


def getdMesh(caseNo):
    M,g,t = loadMmesh(caseNo)
    xslice, yslice, zslice = getSlice(caseNo)

    x2d=M[0][yslice,xslice,zslice]
    y2d=M[1][yslice,xslice,zslice]
    
    dcell = np.diff(g[0]).mean()
    x0 = x2d[0][0]
    
    xi = x2d - x0 + dcell/2
    yi = y2d
    
    return (xi,yi), dcell,t

# Interpolate all fluent grids in folder onto ordered grid and stack into 
# array with order Conc[time][x][y] 
def packDfolder(folderpath):    
    fileindex = sortDfolder(folderpath)
    
    caseNo = int(folderpath.split('/')[-4].split('case')[1])
    
    (xi,yi), dcell,t_mf = getdMesh(caseNo)
    
    zlist = []
    t = []
    for sampleNo,post in enumerate(fileindex):
        time = np.round(post[0],4)
        if time in np.round(t_mf,decimals=4):
            
            x,y,C = loadTimeStep(fileindex,sampleNo)
            
            Ci = griddata((x,y),C,(xi,yi),method='nearest')
            
            # create mask to remove values "inside wall"
            if wLimits['hasWall'][caseNo]:
    
                # Spec wall extent for masking 
                wallxmin = wLimits['w0'][caseNo]
                wallxmax = wLimits['w1'][caseNo]
                wally = wLimits['wy']
    
                mask = (xi > wallxmin) & (xi < wallxmax) & (yi > wally)
                Ci[mask] = np.nan
        
            zlist.append(Ci) 
            t.append(time)
            
            
    zM = np.array(zlist)
    t = np.array(t)
    return zM,(xi,yi),t

# shoehorn all relevant data into numpy-arrays and save using .npy format
# files are prefixed with d and stored according to case info.    
def saveData(caseNo,q,dnumber,overwrite=True,include_static=True,dim1=False):
    
    # create path for nested folder-structure
    fspec= 'case{}/q{}/'.format(caseNo,q)
    
    runfolder = runpath + fspec + 'd{}/'.format(dnumber)
    savefolder = savepath + fspec 
   
    
    # pack the folder and do minor rearrangement of the data before saving.
    zM,grid,t = packDfolder(runfolder)

    # prefix created files with their d-number 
    prefix = savefolder + '{}_'.format(dnumber)
   
    # make folder structure
    os.makedirs(savefolder,exist_ok=True)
    
    if not overwrite:
        print('Case already exists, not overwritten')
        return
    
    
    # save concentration fields
    if dim1:
        z1d = np.nanmean(zM,axis=1)
        np.save(prefix + '1D-Conc',z1d)
    else:
        np.save(prefix + '2D-Conc',zM)
        
    return

# check what k's we have already parsed through, 
# returns set which can be checked against
def getParsed():
    exclude=set(['gridfiles','mfdata','tdata.csv','backup'])
    done = set()
    caselist = [case for case in os.listdir(savepath) if not case in exclude]
    for case in caselist:
        casenum=case.split('case')[1]
        qlist = os.listdir(savepath + case)

        for q in qlist:
            klist = [file.split('_')[0] for file in os.listdir(savepath + case + '/' + q) if not file in exclude]
            #remove duplicates
            klist = list(set(klist))
            for k in klist:
                runstring = 'c{}_{}_d{}'.format(casenum,q,k)
                done.add(runstring)
    return done

# Check to see what runs we have stored in the run directory
# returns set which can be checked against
def getRuns():
    done = set()
    caselist = [case for case in os.listdir(runpath)]
    for case in caselist:
        casenum=case.split('case')[1]
        qlist = os.listdir(runpath + case)
        for q in qlist:
            klist = [file for file in os.listdir(runpath + case + '/' + q)]
            for k in klist:
                runstring = 'c{}_{}_{}'.format(casenum,q,k)
                done.add(runstring)
    return done

def rmRun(caseNo,q,dnumber):
    runstring = 'c{}_q{}_d{}'.format(caseNo,q,dnumber)
    runs = getRuns()
    if runstring in runs:
        ans = input('Really clear run {}? (y/n)'.format(runstring))
    
        if ans == 'y':     
            rundir = '{}case{}/q{}/d{}/'.format(runpath,caseNo,q,dnumber)
            shutil.rmtree(rundir)
        else:
            print('run not removed')
        return
    print('Run {} is not in collection, cannot remove.'.format(runstring))
    return

# -------------------------------------------------------------------------
# Functions for working with Fluent 

## Request k to be added to parsed collection. Runs fluent if needed.
def addRun(caseNo,dnumber,q,verbose=True, keepruns=False,dim1=False):
    
    runstring = 'c{}_q{}_d{}'.format(caseNo,q,dnumber)
    parsed = getParsed()
    runs = getRuns()
    ecode = 0
    if verbose:
        print('\n')
    if runstring in parsed:
        if verbose:
            print('requested case: ' + runstring + ' is already in parsed collection')
        return 0
    elif runstring in runs:
        if verbose:
            print('parsing saved run for case: ' + runstring)
        saveData(caseNo,q,dnumber,dim1=dim1)
    
    else: 
        if verbose:
            print('Making fluent run for case: ' + runstring)
            print('please wait...')
        try:
            ecode = requestK(caseNo,dnumber,q)
        except NameError:
            from d_generator import runCase as requestD
            ecode = requestD(caseNo,dnumber,q)
                
        saveData(caseNo,q,dnumber,dim1=dim1)
    
    if verbose:
        if ecode == 0:
            print('\nSuccesfully added case ' + runstring + ' to collection')
        else: 
            print('\n case generation failed')
    
    if not keepruns:
        rundir = '{}case{}/q{}/d{}/'.format(runpath,caseNo,q,dnumber)
        shutil.rmtree(rundir)
    return ecode

# add cases in a grid defined by the values in klist and qlist
def addCaseGrid(case,dlist,qlist,verbose=False,dim1=False):
    N=len(qlist)*len(dlist)
    tstart=datetime.now().time()
    print('{} - Processing {} cases'.format(tstart,N))
    i=1
    for q in qlist:
        for d in dlist:
            print('doing case {} / {}'.format(i,N))
            addRun(case,np.round(d,7),np.round(q,5),verbose,dim1=dim1)
            
            i +=1
        print ('current time:')
        print(datetime.now().time())
        print('started:')
        print(tstart)
        print('fraction done')
        print(float(i)/float(N))
       
    return

# function to remove all runs from fluent runfolder
def rmAllRuns():
    excepts=set('backup')
    caselist = [case for case in os.listdir(runpath) if case not in excepts]
    
    ans = input('Really clear all D runs? (y/n)')
    
    if ans == 'y':
        for case in caselist:
            shutil.rmtree(runpath + case)
            os.mkdir(runpath + case)
        print('All runs cleared')
    else:
        print('Runs not cleared')
    return

def rmLogs(force=False):
    logpath=runpath.split('runs/')[0] + 'logs/'
    
    if force:
        ans='y'
    else:    
        ans = input('Really clear all D logs? (y/n)')
    
    if ans == 'y':
        shutil.rmtree(logpath)
        os.mkdir(logpath)
        print('Logfolder cleared')
    else:
        print('Logs not cleared')
    return
    

# function to remove all parsed runs from save folder
def rmAllParsed():
    excepts=set(['gridfiles','mfdata','backup'])
    caselist = [case for case in os.listdir(savepath) if case not in excepts]
    
    ans = input('Really clear D collection? (y/n)')
    
    if ans == 'y':
        for case in caselist:
            shutil.rmtree(savepath + case)
            os.mkdir(savepath + case)
        print('Collection cleared')
    else:
        print('Collection not cleared')
    return

# make a backup tarfile of selected folders in parsed collection
def makeBackup(backupname,dirlist):
    TAR_FILENAME = savepath +'backup/' + backupname + '.tar.gz'

    current= os.getcwd()
    os.chdir(savepath)
    
    with tarfile.open(TAR_FILENAME,'w:gz', format=tarfile.PAX_FORMAT) as tar_file:
        print('Creating archive... Please wait..\n' )
        print('Time at start: {}'.format(datetime.now().time()))
        for folder in dirlist:
            tar_file.add(folder)
            print('{} - Done adding: "{}"'.format(datetime.now().time(),folder))

    os.chdir(current)
    return

# restore a certain backup of selected folders
def restoreBackup(backupname):
    print('\nExtracting tarball... Please wait... \n')
    current= os.getcwd()
    
    tar_filename=savepath + 'backup/' + backupname + '.tar.gz'
    tar = tarfile.open(tar_filename)
    
    os.chdir(savepath)
    tar.extractall()
    tar.close()
    os.chdir(current)
    print('Backup restored!')
    return

# -----------------------------------------------------------------------------
# just for testing

# function to retrieve Temp-data in 2-D for a given k-case.
def loadD(caseNo,dnumber,q): 
    runstring = 'c{}_q{}_d{}'.format(caseNo,q,dnumber)
    if not runstring in getParsed():
        print('\ncase: ' + runstring + ' not in parsed collection \n Cannot load.')
        ans = input('\nMake fluent run for it? (y/n) ')
        if ans == 'y':
            ecode=addRun(caseNo,dnumber,q)
            if ecode == 1:
                print('Run failed for case: {} \n Case not added\n Cannot load.'.format(runstring))
                return
        else:
            print('Not running case, Cannot load')
            return
            
    
    fspec= 'case{}/q{}/'.format(caseNo,q)
    savefolder = savepath + fspec
    filename = savefolder + '{}_2D-Conc.npy'.format(dnumber)
    
    return np.load(filename)


# -----------------------------------------------------------------------------




















