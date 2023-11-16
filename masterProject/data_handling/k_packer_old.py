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

from scipy.interpolate import griddata


import os
import shutil
import tarfile




# -------------------------------------------------------
# Paths
os.chdir('/chalmers/users/gabgus/python-scripts/masterProject')

# path to save structured collection of data
savepath= '/scratch/gabgus/fluent/kgenerator/parsed/'

# path to folder with saved fluent runs
runpath = '/scratch/gabgus/fluent/kgenerator/runs/' 

# path to fluent runscript and associated journal files
scriptname = '/scratch/gabgus/fluent/kgenerator/krun.sh'
stencilpath = '/scratch/gabgus/fluent/kgenerator/journal_stencil.jou'
journalpath = '/scratch/gabgus/fluent/kgenerator/active_journal.jou'

# path to multiphase simulation data
mfdatapath='/scratch/gabgus/mfdata_all/meshed/'

# -------------------------------------------------------------------------
# Global variables

wLimits =  {'w0': ['error' , np.nan, 0.255, 0.180, 0.225, 0.255, np.nan, 0.255, np.nan],
            'w1': ['error' , np.nan, 0.285, 0.210, 0.255, 0.285, np.nan, 0.285, np.nan],
            'wy': 0.3,
            'hasWall': [False, False, True, True, True, True, False, True, False]}

Glimits = {'x0': ['error' , 0.0, 0.0, 0.0, 0.0, 'missing', 0.0, 0.0, 0.0], 
           'x1': ['error' , 0.54, 0.54, 0.54, 0.54, 'missing', 0.54, 0.54, 0.452],
           'y1': ['error' , 34, 34, 34, 34, 'missing', 28, 28, 80],
           'y0': ['error' , 6, 6, 6, 6, 'missing', 0, 0, 0]}

# -------------------------------------------------------------------------
# Functions

# load multiphase mesh
def loadMmesh(caseNo):
    M = np.load('{}c{}_mesh.npy'.format(mfdatapath,caseNo))
    g = np.load('{}c{}_gridpoints.npy'.format(mfdatapath,caseNo))
    t = np.load('{}c{}_times.npy'.format(mfdatapath,caseNo))
    return M,g,t

# get slices to apply on multiphase mesh make it correspond to k-domain
def getSlice(caseNo):
    Mlimits = {'x0': ['error' , 15, 15, 15, 15, 'missing', 15, 15, 24], 
           'x1': ['error' , 51, 51, 51, 51, 'missing', 51, 51, 76],
           'y1': ['error' , 34, 34, 34, 34, 'missing', 28, 28, 80],
           'y0': ['error' , 6, 6, 6, 6, 'missing', 0, 0, 0]}
    
    xslice=slice(Mlimits['x0'][caseNo],Mlimits['x1'][caseNo])  #max66
    yslice=slice(Mlimits['y0'][caseNo],Mlimits['y1'][caseNo]) #max67
    zslice=1            #max86
    return xslice, yslice, zslice

# parse through and sort contents of folder according to timestep 
# for easy access by later functions
def sortKfolder(folderpath):
    filelist = [folderpath + file for file in os.listdir(folderpath)]
    times = [float(filename.split('-')[1]) for filename in os.listdir(folderpath)]
    return list(sorted(zip(times,filelist)))

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

def getkMesh(caseNo):
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
# array with order Temp[time][x][y] 
def packKfolder(folderpath):    
    fileindex = sortKfolder(folderpath)
    
    caseNo = int(folderpath.split('/')[-4].split('case')[1])
    
    (xi,yi), dcell,t_mf = getkMesh(caseNo)
    
    zlist = []
    t = []
    for sampleNo,post in enumerate(fileindex):
        time = np.round(post[0],4)
        if time in t_mf:
            
            x,y,T = loadTimeStep(fileindex,sampleNo)
            
            Ti = griddata((x,y),T,(xi,yi),method='nearest')
            
            # create mask to remove values "inside wall"
            if wLimits['hasWall'][caseNo]:
    
                # Spec wall extent for masking 
                wallxmin = wLimits['w0'][caseNo]
                wallxmax = wLimits['w1'][caseNo]
                wally = wLimits['wy']
    
                mask = (xi > wallxmin) & (xi < wallxmax) & (yi > wally)
                Ti[mask] = np.nan
        
            zlist.append(Ti)
            t.append(time)
    zM = np.array(zlist)
    t = np.array(t)
    return zM,(xi,yi),t


# shoehorn all relevant data into numpy-arrays and save using .npy format
# files are prefixed with k and stored according to case info.    
def saveData(caseNo,kfrac,knumber,overwrite=True,include_static=True):
    
    # create path for nested folder-structure
    fspec= 'case{}/q{}/'.format(caseNo,kfrac)
    
    runfolder = runpath + fspec + 'k{}/'.format(knumber)
    savefolder = savepath + fspec 
   
    
    # pack the folder and do minor rearrangement of the data before saving.
    zM,grid,t = packKfolder(runfolder)

    # prefix created files with their k-number 
    prefix = savefolder + '{}_'.format(knumber)
   
    # make folder structure
    os.makedirs(savefolder,exist_ok=True)
    
    if not overwrite:
        print('Case already exists, not overwritten')
        return
    
    # save temperature fields
    np.save(prefix + '2D-Temp',zM)
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
                runstring = 'c{}_{}_k{}'.format(casenum,q,k)
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
    

## Function to make sure we have parsed all runs so far
def saveAll():
    parsed=getParsed()
    runs=getRuns()
    notParsed=runs.difference(parsed)
    if len(notParsed) == 0:
        print('All runs parsed, no action needed')
        return
    print('saved case:')
    for runstring in notParsed:
        caseNo=runstring.split('_')[0].split('c')[1]
        kfrac=runstring.split('_')[1].split('q')[1]
        knumber=runstring.split('_')[2].split('k')[1]
        
        saveData(caseNo,kfrac,knumber)
        print(runstring)
    return 


# -------------------------------------------------------------------------
# Functions for working with Fluent 

## Request k to be added to parsed collection. Runs fluent if needed.
def addRun(caseNo,knumber,kfrac,verbose=True):
    
    runstring = 'c{}_q{}_k{}'.format(caseNo,kfrac,knumber)
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
        saveData(caseNo,kfrac,knumber)
    
    else: 
        if verbose:
            print('Making fluent run for case: ' + runstring)
            print('please wait...')
            try:
                ecode = requestK(caseNo,knumber,kfrac)
            except NameError:
                from k_generator import runCase as requestK
                ecode = requestK(caseNo,knumber,kfrac)
                
            saveData(caseNo,kfrac,knumber)
    
    if verbose:
        if ecode == 0:
            print('\nSuccesfully added case ' + runstring + ' to collection')
        else: 
            print('\n case generation failed')
    return ecode

# run list of cases each case in the format [caseNo,knumber,kfrac]
def addList(caselist):
    for case in caselist:
        caseNo=case[0]
        knumber=case[1]
        kfrac=case[2]
        addRun(caseNo,knumber,kfrac)
    return

# add cases in a grid defined by the values in klist and qlist
def addCaseGrid(case,klist,qlist,verbose=True):
    N=len(qlist)*len(klist)
    print('Processing {} cases'.format(N))
    i=1
    for q in qlist:
        for k in klist:
            print('doing case {} / {}'.format(i,N))
            addRun(case,k,q,verbose)
            i +=1
    return

# function to remove all runs from fluent runfolder
def rmAllRuns():
    excepts=set('backup')
    caselist = [case for case in os.listdir(runpath) if case not in excepts]
    
    ans = input('Really clear all runs? (y/n)')
    
    if ans == 'y':
        for case in caselist:
            shutil.rmtree(runpath + case)
            os.mkdir(runpath + case)
        print('All runs cleared')
    else:
        print('Runs not cleared')
    return

# function to remove all parsed runs from save folder
def rmAllParsed():
    excepts=set(['gridfiles','mfdata','backup'])
    caselist = [case for case in os.listdir(runpath) if case not in excepts]
    
    ans = input('Really clear collection? (y/n)')
    
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
        for folder in dirlist:
            tar_file.add(folder)
            print('Done adding: "{}"'.format(folder))

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

#-------------------------------------------------------------------------- 
# Retrieval functions
                    
# function to retrieve Temp-data 2-D for a given k-case.
def loadK(caseNo,knumber,kfrac):
    
    runstring = 'c{}_q{}_k{}'.format(caseNo,kfrac,knumber)
    if not runstring in getParsed():
        print('\ncase: ' + runstring + ' not in parsed collection \n Cannot load.')
        ans = input('\nMake fluent run for it? (y/n) ')
        if ans == 'y':
            ecode=addRun(caseNo,knumber,kfrac)
            if ecode == 1:
                print('Run failed for case: {} \n Case not added\n Cannot load.'.format(runstring))
                return
        else:
            print('Not running case, Cannot load')
            return
            
    
    fspec= 'case{}/q{}/'.format(caseNo,kfrac)
    
    savefolder = savepath + fspec 
   
    filename = savefolder + '{}_2D-Temp.npy'.format(knumber)
    
    return np.load(filename)


# --------- Functions for examining boundary temps ----------------------------------------------- 

# Fetch lists with boundary temperatures for each case over range of samples
def getBoundaryTemp(caseNo,s0,s1,knumber=10000,kfrac=0.5):
    kdata,tfreq,tscaling,grid = loadK(caseNo,knumber,kfrac)
    T0 = []
    T1 = []
    kT0 = []
    kT1 = []
    for sno in range(s0,s1):
        mx,mT,_,time = fetchMdata(caseNo,sno)
        T0.append(mT[0])
        T1.append(mT[-1])
        ts = sno +1
        kT0.append(kdata[ts][1])
        kT1.append(kdata[ts][-2])
    
    return np.array(T0), np.array(T1), np.array(kT0), np.array(kT1)

def calcAverageTempDiff(T0,T1,kT0,kT1):
    
    T0bar = sum(T0[100:-1])/len(T0[100:-1])
    T1bar = sum(T1[100:-1])/len(T1[100:-1])


    kT0bar = sum(kT0[100:-1])/len(kT0[100:-1])
    kT1bar = sum(kT1[100:-1])/len(kT1[100:-1])

    dt0 = kT0bar-T0bar
    dt1 = kT1bar - T1bar
    return dt0, dt1


def getCorrection(w,knumber=10000,kfrac=0.5):
    dt0=[0,0,0,0,0]
    dt1=[0,0,0,0,0]
    dtavg=[0,0,0,0,0]

    for c in [1,2,3,4]:
        T0,T1,kT0,kT1=getBoundaryTemp(c,0,230,knumber,kfrac)
        dt0[c],dt1[c] = calcAverageTempDiff(T0,T1,kT0,kT1) 
        w = 0.5
        dtavg[c] = (w*dt0[c]- (1-w)*dt1[c])  

    return dt0, dt1, dtavg

# ------------------ Optimization -----------------------------------------

# Function to compute the error between a generated k-profile and 
# its corresponding multiphase case.
def getDeviation(caseNo,knumber,kfrac,sI=[0,229],xI=[0.02,0.52]):
    
    mx = loadMdata(caseNo,'x')
    mT = loadMdata(caseNo,'T')
    
    kx = loadK(caseNo,knumber,kfrac,'x')
    kT= loadK(caseNo,knumber,kfrac,'T')
    
    mxmask = (mx >= xI[0]) & (mx <= xI[1])
    kxmask = (kx >= xI[0]) & (kx <= xI[1])
    
    ns= sI[1] -sI[0]
    
    tCumErr2=np.zeros((ns,))
    
    n = 0
    
    for sno in range(sI[0],sI[1]):
        xErr2 = (kT[sno][kxmask] - mT[sno][mxmask])**2
        tCumErr2[n]=xErr2.sum()
        n += 1
        
    caseErr2=tCumErr2.sum()
    return caseErr2, len(kT[0][kxmask]), len(range(sI[0],sI[1]))

#-------------------- Plotting --------------------------------------------

# Plots wall and low k zone overlays together with some other notation etc.
def furbish(caseNo,time):
    wxmin, wxmax, lkmin, lkmax = wLimits(caseNo)
    
    if not caseNo == 1:
        plt.plot([wxmin,wxmin],[1123.15,1073.15],'--',color=[0.8,0.8,0.8])
        plt.plot([wxmax,wxmax],[1123.15,1073.15],'--',color=[0.8,0.8,0.8])
    
        plt.plot([lkmin,lkmin],[1123.15,1073.15],'--',color=[0.8,0.8,0.8])
        plt.plot([lkmax,lkmax],[1123.15,1073.15],'--',color=[0.8,0.8,0.8])
        
    plt.xlabel('x-coordinate [m]')
    plt.ylabel('Temperature [K]')
    plt.title('Temperature gradient at t = ' + str(time) + 's')
    plt.xlim(0.0,0.54)
    return

# plot a chosen k-solution together with its corresponding mf-case 
# at a certain timestep. Also marks the wall and low-k zone for each case.
    
def plotComparison(caseNo,sampleNo,kfrac,knumber):
    
    kdata,tfreq,tscaling,grid = loadK(caseNo,knumber,kfrac)

    mx, mT,_, time = fetchMdata(caseNo,sampleNo)

    ts = sampleNo + 1    
    
    kT=kdata[ts]

    plt.plot(grid,kT,label='k = {} W/mK'.format(knumber))
    plt.plot(mx,mT,'r.',label='multiphase data')
    plt.legend()
    furbish(caseNo,time)
    return

def plotK(caseNo,sampleNo,kfrac,knumber):
    kdata,tfreq,tscaling,grid = loadK(caseNo,knumber,kfrac)
    ts = sampleNo + 1    
    
    kT=kdata[ts]
    time = ts*tscaling*tfreq
    plt.plot(grid,kT,label='k = {} W/mK'.format(knumber))
    furbish(caseNo,time)
    return

# Plot the response surface for a grid of runs specified by kpoint and qpoints
# Can generate new runs if specified to.
def plotCaseGrid(caseNo,krange,qrange,sI=[0,229],xI=[0.02,0.52],generate=False):
 
    K,Q = np.meshgrid(krange,qrange)
    if generate:
        addCaseGrid(caseNo,krange,qrange,verbose=False)
        
    ERR = np.zeros_like(K)
    
    for j in range(len(krange)):
        for i in range(len(qrange)):
            ERR[i][j],nx,nt = getDeviation(caseNo,krange[j],qrange[i],sI=sI)
            
    
    imin=np.argmin(ERR)
    
    idx=np.unravel_index(imin,ERR.shape)
    plt.figure(1,(6,5))
    plt.contourf(K,Q,ERR,80)
    plt.colorbar()
    plt.contour(K,Q,ERR,20)
    plt.plot(krange[idx[1]],qrange[idx[0]],'ro')
    plt.xlabel('Conduction [W/mK]')
    plt.ylabel('$k_1 / k_2$')
    plt.grid()
    plt.savefig('/scratch/gabgus/images_etc/test.pdf', bbox_inches='tight')
    plt.show()
    
    return K,Q,ERR,idx


def plotBoundaryTemp(caseNo,sI=[0,229],yI=[1070,1090]):
    T0,T1,kT0,kT1 = getBoundaryTemp(1,sI[0],sI[1])

    Sno = np.array([i for i in range(sI[0],sI[1])])

    plt.figure()
    plt.plot(Sno,T0,Sno,T1)
    
    dt0,dt1,dtavg=getCorrection(0.5)


    dT=dtavg[caseNo]
    plt.plot(Sno+1,kT0-dT,Sno+1,kT1+dT)

    plt.plot([Sno[0],Sno[-1]],[1123.15,1123.15],'r--')

    plt.plot([Sno[0],Sno[-1]],[1073.15,1073.15],'b--')

    plt.xlim(sI[0],sI[1])
    plt.ylim=(yI[0],yI[1])
    return

def returnBestOnGrid(caseNo,nk=20,nq=25,kstep=250,qstep=0.05,generate=False):
    # some convenience
    if type(nk) == int:
        kI = [0,nk]
    else:
        kI = nk
    
    if type(nq) == int:
        qI=[0,nq]
    else:
        qI=nq

    # define gridpoints
    krange=[int(kstep + kstep*i) for i in range(kI[0],kI[1])]
    qrange=[round(qstep + qstep*i,4) for i in range(qI[0],qI[1])]
    
    # plot grid and return minimum
    K,Q,ERR,idx=plotCaseGrid(caseNo,krange,qrange,sI=[0,229],xI=[0.02,0.52],generate=generate)
    
    kfit = krange[idx[1]]
    qfit = qrange[idx[0]]
    err = ERR[idx[0]][idx[1]]
        
    return kfit, qfit, err, k_uni(kfit,qfit)


# -------------- Main loop ---------------------------------------------------
#%%

testfolder = '/scratch/gabgus/fluent/kgenerator/runs/case2/q1.0/k2000/'

TM,M,t = packKfolder(testfolder)

Ti=TM[150]

g=(M[0][0],M[1].T[0])

ybar = np.nanmean(Ti,axis=0)

def plot2d(x2d,y2d,zbar,time,quantity='temperature'):
    plt.contourf(x2d,y2d,zbar,30)
    plt.colorbar()
    plt.title(time)

    plt.contour(x2d,y2d,zbar,40,color='k')
    plt.xlabel('x-coordinate [m]')
    plt.ylabel('y-coordinate [m]')
    plt.title('Z-averaged {} at {}s'.format(quantity,time))
    plt.show()
    return

plot2d(M[0],M[1],Ti,t[150])

def plot1d(gx,ybar,time,quantity='temperature'):
    plt.plot(gx,ybar)
    plt.xlim(gx[0],gx[-1])
    plt.title('Average {} along x at {}s'.format(quantity,time))
    plt.xlabel('x-coordinate [m]')
    plt.ylabel('{}'.format(quantity))
    plt.show()
    return

plot1d(g[0],ybar,t[150])

#
#def lossFunc(caseNo,theta,kI=kI,qI=qI):
#    # Data ranges in time and space
#    sI=[0,229],
#    xI=[0.02,0.52]
#    
#    # Dimensionalising
#    knumber=int((1-theta[0])*kI[0] + theta[0]*kI[1])
#    kfrac=round((1-theta[1])*qI[0] + theta[1]*qI[1],4)
#    
#    addRun(caseNo,knumber,kfrac)
#    
#    return getDeviation(caseNo,knumber,kfrac)*1e-5
#
#nmax  = 20
##def fdsaSearch(caseNo,knumber,kfrac,nmax=20):   
## Non-dimensionalising
#kstar= (knumber - kI[0])/(kI[1] - kI[0])
#qstar= (kfrac - qI[0])/(qI[1] - qI[0])
#
#theta = np.array([kstar,qstar])
#
#thetalist = []
#thetalist.append(theta)
#
#print('\nRunning search')
#for k in range(nmax):
#    print('k = ' + str(k))
#    alpha = 0.602
#    gamma = 0.101
#
#    A=1.0
#    a=0.1
#    c=0.01
#
#    ak=a/((k+1+A)**alpha)
#    ck=c/((k+1)**gamma)
#        
#    # Defining evaluation points, 
#    # u=upper, e=east, l=lower, w=west 
#    ue = theta + ck*np.array([1,0])
#    uw = theta - ck*np.array([1,0])
#    le = theta + ck*np.array([0,1])
#    lw = theta - ck*np.array([0,1])
#    
#    
#    # Finite difference gradient estimation
#    grad=np.array([(lossFunc(caseNo,ue) - lossFunc(caseNo,uw))/(2*ck),
#                   (lossFunc(caseNo,le) - lossFunc(caseNo,lw))/(2*ck)])
#    
#    #  stepping to next value
#    theta_new = theta - ak*grad
#    
#    theta=theta_new
#    
#    thetalist.append(theta)
#
## Dimensionalising
#knumber=(1-theta[0])*kI[0] + theta[0]*kI[1]
#kfrac=(1-theta[1])*qI[0] + theta[1]*qI[1]
    
#return knumber, kfrac, thetalist















## function to fetch limits for the wall for each case
#def wLimits(caseNo):
#    if caseNo not in [2,3,4]:
#        return 0,0,0,0
#    elif caseNo == 1:
#        lkmin = 0.015*(18 - 5)
#        lkmax = 0.015*(18 + 5)
#        return 0,0,lkmin,lkmax
#    
#    wallpos = [0,0,18,13,16]
#    
#    dcell = 0.015
#    
#    wxmin = dcell*(wallpos[caseNo] - 1)
#    wxmax = dcell*(wallpos[caseNo] + 1)
#    
#    lkmin = dcell*(wallpos[caseNo] - 5)
#    lkmax = dcell*(wallpos[caseNo] + 5)
#    return wxmin,wxmax,lkmin,lkmax




def k_uni(k,q):
    dcell = 0.015
    
    L=34*dcell
    L15=(16-4)*2*dcell
    L24=8*dcell
    L3=2*dcell
    
    A1 = 1.3*0.42
    A2 = 0.42-0.15
    alpha = A2/A1
    
    kuni=k*q*alpha/(q*alpha*L15/L + alpha*L24/L + L3/L)
    return kuni



#----- Functions for handling MF Data ---------------------------------------        

# list multiphase folder and return sorted array of times and filenames
#def sortMFolder(folderpath,tscaling,toffset):
#    filelist=os.listdir(folderpath)
#
#    tsteps = []
#    for file in filelist:
#        tsteps.append(file.split('_')[2].split('.')[0])
#
#    times = [round(float(tstep)*tscaling + toffset,4) for tstep in tsteps]
#
#    fileindex = list(sorted(zip(times,filelist)))
#    return fileindex
#
## read a certain timestep file in a certain casefolder and return lists 
#def fetchMdata(caseNo,sampleNo):
#
#    folderpath = mfdatapath + 'case0' + str(caseNo) + '/'
#
#
#    # value to scale integer from filename with
#    tscaling = 1e-6
#
#    # value to add to captured time to offset time before averaging process
#    toffsets=[0,-3.6,-3.3,-3.3,-3.5]
#    toffset = toffsets[caseNo]
#
#
#    fileindex = sortMFolder(folderpath,tscaling,toffset)
#
#    filename = fileindex[sampleNo][1]
#    time = fileindex[sampleNo][0]
#    df = pd.read_csv(folderpath+filename,names=['x','Ti','Tavg'])
#
#    # unpack into arrays
#    x = np.array(df.x)
#    Ti = np.array(df.Ti)
#    Tavg = np.array(df.Tavg)
#    
#    dcell = np.diff(x)[0]
#    x0 = x[0]
#    x = x - x0 + dcell/2 + dcell
#    
#    
#    return x, Ti, Tavg, time
#
#
#def packMdata(caseNo):
#    T = []
#    t = []
#    for sno in range(229):
#        mx,mT,mTavg,time = fetchMdata(caseNo,sno)
#        T.append(mT)
#        t.append(time)
#        
#    prefix = savepath + 'mfdata/c{}_'.format(caseNo)
#    np.save(prefix + 'T',np.array(T))
#    np.save(prefix + 't',np.array(t))
#    np.save(prefix + 'x',mx)
#    return
#
#



    
#%%
## ------------ Plotting -------------------------------------
#

    

#plt.plot(k,ck)


## plot interpolated grid as surface
#plt.contourf(xi,yi,z0,10,cmap='bwr', vmin=1073.15, vmax=1123.15)
#plt.colorbar()
#
#plt.contour(xi,yi,z0,33,colors='black', vmin=1073.15, vmax=1123.15)caseNo=2casecaseNo=2

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


