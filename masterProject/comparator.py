#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 18:10:56 2018

Script to compare data from multiphase simulations with analogous data 
from conduction or dispersion driven cases and compute the difference
across several different timesteps

@author: gabgus
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os

#from k_packer import getParsed, addRun
# -------------------------------------------------------------------------

mfdatapath='/scratch/gabgus/mfdata_all/meshed/'
kfolder = '/scratch/gabgus/fluent/kgenerator/parsed/'
dfolder = '/scratch/gabgus/fluent/dgenerator/parsed/'

dispnames= {'xvel' : 'x velocity',      'yvel' : 'y velocity', 'zvel' :'z velocity',
                 'uds' : 'passive scalar',  'temp' : 'temperature', 'vof' : 'volume fraction' }

Klimits = {'x0': ['error' , 15, 15, 15, 15, 'missing', 15, 15, 24], 
           'x1': ['error' , 51, 51, 51, 51, 'missing', 51, 51, 76],
           'y1': ['error' , 38, 38, 38, 38, 'missing', 32, 32, 60],
           'y0': ['error' , 6, 6, 6, 6, 'missing', 0, 0, 0]}

Dlimits = {'x0': ['error' , 0, 0, 0, 0, 'missing', 0, 0, 0], 
           'x1': ['error' , 67, 67, 67, 67, 'missing', 67, 67, 100],
           'y1': ['error' , 38, 38, 38, 38, 'missing', 32, 32, 60],
           'y0': ['error' , 6, 6, 6, 6, 'missing', 0, 0, 0]}

mtemp=[1098.15,1098.15,1098.15,1098.15,1098.15,1098.15,1098.15,1098.15,1173.15]

# -------------------------------------------------------------------------

def loadMdata(caseNo,scalar='temp'):
    if scalar == 'grid':
        M = np.load('{}c{}_mesh.npy'.format(mfdatapath,caseNo))
        t = np.load('{}c{}_times.npy'.format(mfdatapath,caseNo))
        return M,t
        
    Smat = np.load('{}c{}_{}.npy'.format(mfdatapath,caseNo,scalar))
    M = np.load('{}c{}_mesh.npy'.format(mfdatapath,caseNo))
    t = np.load('{}c{}_times.npy'.format(mfdatapath,caseNo))
    return Smat,M,t
    

def getSlice(caseNo,dk='k'):
    if dk == 'k':
        xslice=slice(Klimits['x0'][caseNo],Klimits['x1'][caseNo])  #max66
        yslice=slice(Klimits['y0'][caseNo],Klimits['y1'][caseNo]) #max67
        zslice=1 #max86
    else:
        xslice=slice(Dlimits['x0'][caseNo],Dlimits['x1'][caseNo])  #max66
        yslice=slice(Dlimits['y0'][caseNo],Dlimits['y1'][caseNo]) #max67
        zslice=1 #max86
           
    return xslice, yslice, zslice

# -------------------------------------------------------------------------
    # K functions
    
def getkParsed():
    from k_packer import getParsed
    done = getParsed()
    return done

def resetKlib():
    from k_packer import rmAllParsed, rmAllRuns, rmLogs
    rmAllParsed()
    rmAllRuns()
    rmLogs()
    return




def addkgrid(caseNo,klist,qlist,verbose=False):
    from k_packer import addCaseGrid
    addCaseGrid(caseNo,klist,qlist,verbose)
    return





# function to retrieve Temp-data in 2-D for a given k-case.
def loadK(caseNo,knumber,kfrac): 
    runstring = 'c{}_q{}_k{}'.format(caseNo,kfrac,knumber)
    if not runstring in getkParsed():
        print('\ncase: ' + runstring + ' not in parsed collection \n Cannot load.')
        ans = input('\nMake fluent run for it? (y/n) ')
        if ans == 'y':
            from k_packer import addRun as addKRun
            ecode=addKRun(caseNo,knumber,kfrac)
            if ecode == 1:
                print('Run failed for case: {} \n Case not added\n Cannot load.'.format(runstring))
                return
        else:
            print('Not running case, Cannot load')
            return
            
    
    fspec= 'case{}/q{}/'.format(caseNo,kfrac)
    
    savefolder = kfolder + fspec 
   
    filename = savefolder + '{}_2D-Temp.npy'.format(knumber)
    
    return np.load(filename)



# -------------------------------------------------------------------------
    # D functions

def getdParsed():
    from d_packer import getParsed
    done = getParsed()
    return done


def resetDlib():
    from d_packer import rmAllParsed, rmAllRuns, rmLogs
    rmAllParsed()
    rmAllRuns()
    rmLogs()
    return

def adddgrid(caseNo,dlist,qlist,verbose=False):
    from d_packer import addCaseGrid
    addCaseGrid(caseNo,dlist,qlist,verbose)
    return



# function to retrieve Conc-data in 2-D for a given d-case.
def loadD(caseNo,dnumber,dfrac): 
    runstring = 'c{}_q{}_d{}'.format(caseNo,np.round(dfrac,4),np.round(dnumber,7))
    if not runstring in getdParsed():
        print('\ncase: ' + runstring + ' not in parsed collection \n Cannot load.')
        ans = input('\nMake fluent run for it? (y/n) ')
        if ans == 'y':
            from d_packer import addRun as addDRun
            ecode=addDRun(caseNo,dnumber,dfrac)
            if ecode == 1:
                print('Run failed for case: {} \n Case not added\n Cannot load.'.format(runstring))
                return
        else:
            print('Not running case, Cannot load')
            return
            
    
    fspec= 'case{}/q{}/'.format(caseNo,dfrac)
    
    savefolder = dfolder + fspec 
   
    filename = savefolder + '{}_2D-Conc.npy'.format(dnumber)
    
    return np.load(filename)



# -------------------------------------------------------------------------
# plotting

class MidpointNormalize(matplotlib.colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        matplotlib.colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))


def plot2d(x2d,y2d,zbar,time,scalar='temp', cmap = 'RdBu_r'):
    midpoints= {'xvel' : 0.0,   'yvel' : 0.0,    'zvel' : 0.0,
                 'uds' : 0.5,   'temp' : 1098.15, 'vof' : 0.63/2 }
    dispnames= {'xvel' : 'x velocity', 'yvel' : 'y velocity', 'zvel' :'z velocity',
                 'uds' : 'passive scalar',   'temp' : 'temperature', 'vof' : 'volume fraction' }
    
    if np.isscalar(time):
        timestring = 'at {}s'.format(time)
    elif time[0] == time[-1]:
        timestring = 'at {}s'.format(time[0])
    else:
        timestring = 'from {}s to {}s'.format(time[0],time[-1])
    
    if scalar == 'vof':
        vmin = 0.0
        vmax = 0.63
    elif scalar == 'uds':
        vmin = 0.0
        vmax =1.0
    else:
        
        valmin=np.abs(np.nanmin(zbar))
        valmax=np.abs(np.nanmax(zbar))
       
        if scalar == 'temp':
            if valmin > midpoints[scalar]:
                midpoints[scalar]=1173.15
       
        
        maxval = max(valmin,valmax)
        maxdiff = maxval - midpoints[scalar]
                
        
        
        vmin = midpoints[scalar] - maxdiff
        vmax = midpoints[scalar] + maxdiff
        
    
    
    norm=MidpointNormalize(vmin,vmax, midpoint=midpoints[scalar]) 
    
    plt.contourf(x2d,y2d,zbar,20,cmap=cmap,norm=norm)
    plt.colorbar()
    plt.contour(x2d,y2d,zbar,60,cmap=cmap,norm=norm)
    plt.xlabel('x-coordinate [m]')
    plt.ylabel('y-coordinate [m]')
    plt.title('Z-averaged {} {}'.format(dispnames[scalar],timestring))

    return
    


def plot1d(gx,ybar,time,scalar='temp',label='1d-data',style='o-'):
    dispnames= {'xvel' : 'x velocity',      'yvel' : 'y velocity', 'zvel' :'z velocity',
                 'uds' : 'passive scalar',  'temp' : 'temperature', 'vof' : 'volume fraction' }
    
    if np.isscalar(time):
        timestring = 'at {}s'.format(time)
    elif time[0] == time[-1]:
        timestring = 'at {}s'.format(time[0])
    else:
        timestring = 'from {}s to {}s'.format(time[0],time[-1])
    
    plt.plot(gx,ybar,style,label=label)
    plt.title('Average {} along x {}'.format(dispnames[scalar],timestring))
    plt.xlabel('x-coordinate [m]')
    plt.ylabel('{}'.format(dispnames[scalar]))
  
    return

def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def plotQsignal(caseNo,plot=True):
    import scipy.stats as st

    qin = np.load('{}c{}_qin.npy'.format(mfdatapath,caseNo))
    qut = np.load('{}c{}_qut.npy'.format(mfdatapath,caseNo))
    qtime = np.load('{}c{}_qtime.npy'.format(mfdatapath,caseNo))
    qmid =(qin+qut)/2

    A=[0,0.6,0.6,0.6,0.6,0,0.6,0.6,0.2]
    dx = [0,0.525,0.525,0.525,0.525,0,0.525,0.525,0.306]
    dT = 50.0

    sslice=slice(int(0.8*len(qtime)),len(qtime))

    q = qmid[sslice]
    qtime=qtime[sslice]
    
    qx = q/A[caseNo]    
    keff = qx * dx[caseNo] / dT

    q_int = st.t.interval(0.999, len(qx)-1, loc=qx.mean(), scale=st.sem(qx))
    k_int = st.t.interval(0.999, len(keff)-1, loc=keff.mean(), scale=st.sem(keff))



    N = int(0.2*len(qx))
    qroll= moving_average(qx,N)
    
    if plot:
        plt.plot(qtime,qx)
        plt.plot(qtime[N-1:],qroll)
        print('qx = {} MW/m2'.format(qx.mean()/1e6))
        plt.show()
        plt.plot(qtime,keff)
    return qx.mean(), q_int, keff, keff.mean(),k_int


def makeFit(caseNo,Smat,M,Sno,
            xslice = 'k', yslice = 'k', zslice = 1):
    
    if xslice == 'k':
        xslice=slice(Klimits['x0'][caseNo],Klimits['x1'][caseNo])  #max66
    if yslice == 'k':    
        yslice=slice(Klimits['y0'][caseNo],Klimits['y1'][caseNo]) #max67
    
    if xslice == 'd':
        xslice=slice(Dlimits['x0'][caseNo],Dlimits['x1'][caseNo])  #max66
    if yslice == 'd':    
        yslice=slice(Dlimits['y0'][caseNo],Dlimits['y1'][caseNo]) #max67
    
    
    x2d=M[0][yslice,xslice,zslice]
    y2d=M[1][yslice,xslice,zslice]
    
    Si = Smat[[Sno]].mean(axis = 0)
    zavg = Si[yslice,xslice,:].mean(axis = 2)
    return x2d,y2d,zavg

def getdeviation2d(caseNo,knumber,q,scalar='temp'):
    Tmat,M,t = loadMdata(caseNo,scalar)
    if scalar == 'temp':
        kdata = loadK(caseNo,knumber,q)
        xslice,yslice, zslice = getSlice(caseNo)
    elif scalar == 'uds':
        kdata = loadD(caseNo,knumber,q)
        xslice,yslice, zslice = getSlice(caseNo,dk='d')
    else:
        print('no stored scalar field for scalar'.format(scalar))
        return
    
    zavg = np.mean(Tmat[:,yslice,xslice,:],axis=3)
    
    errall= np.abs(zavg - kdata)

    terr = np.mean(errall,axis=0)
    xerr = np.nanmean(terr,axis=0)
    errtot = xerr.mean(axis=0)
    return errall,xerr,terr,errtot 

def get1ddata(caseNo,scalar='temp'):    
    Vmat,M,t = loadMdata(caseNo,'vof')
    
    Smat,M,t = loadMdata(caseNo,scalar)
    
    if caseNo < 6:
        yslice=slice(6,250)
        zslice=slice(0,60)
    else:
        yslice=slice(0,250)
        zslice=slice(0,250)
    
    
    vavg = np.mean(Vmat[:,yslice,:,zslice],axis=3)
    zavg = np.mean(Smat[:,yslice,:,zslice],axis=3)    
    
    zdot=zavg*vavg
    
    vofsum=np.nansum(vavg,axis=1)
    zdotsum=np.nansum(zdot,axis=1)
    ybar = zdotsum/vofsum
    
    
    gx = M[0,0,:,0]
    return ybar,gx,t

def getdeviation1d(caseNo,knumber,q,ybar,scalar='temp',tslice=slice(0,500),xmask=slice(0,500)):
    if scalar == 'temp':
        kdata = loadK(caseNo,knumber,q)
    elif scalar == 'uds':
        kdata = loadD(caseNo,knumber,q)
    else:
        print('no stored scalar field for scalar {}'.format(scalar))
   
    kbar = np.nanmean(kdata[tslice], axis = 1)
    ybar = ybar[tslice]

    if scalar == 'temp':
        terr= np.abs(ybar[:][1:len(ybar)-1] - kbar[:][1:len(ybar)-1])
    else:
        terr=np.abs(ybar[:,xmask] - kbar[:,xmask])

    xerr = np.mean(terr,axis=0)
 
    errtot = xerr.mean(axis=0)
    return terr, xerr, errtot



def getResponse(caseNo,klist,qlist,tslice=slice(0,140),scalar='temp',xmask=slice(0,500)):
    err = np.empty([len(klist),len(qlist)])
    
    ybar, gx,t = get1ddata(caseNo,scalar)
    
    if scalar == 'temp':
        xslice=slice(Klimits['x0'][caseNo],Klimits['x1'][caseNo])
    else:
        xslice=slice(Dlimits['x0'][caseNo],Dlimits['x1'][caseNo])
    
    ybar = ybar[:,xslice]
    
    for i,knumber in enumerate(klist):
        for j,q in enumerate(qlist):
            terr, xerr, errtot = getdeviation1d(caseNo,round(knumber,6),round(q,4),ybar,scalar=scalar,tslice=tslice,xmask=xmask)
            err[i,j]=errtot
    
    X,Y = np.meshgrid(klist,qlist)
    plt.contourf(X,Y,err.T,30)
    plt.colorbar()
    plt.contour(X,Y,err.T,40)

    imin=np.argmin(err)
    
    idx=np.unravel_index(imin,err.shape)
    plt.plot(klist[idx[0]],qlist[idx[1]],'ro')
    plt.show()

    return klist[idx[0]], qlist[idx[1]], (X,Y,err.T)
        


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
    


# -------------------------------------------------------------------------

#%%3
# base grid=    d: 0.001 - 0.001 - -> 0.034 
#               q: 0.1 - 0.1 --> 1.2
#

#for cno in [1,2,3,4]:
#    #dlist = np.append(np.array([0.0001]),np.arange(0.0005,0.0045,0.0005))
#    #qlist = np.arange(np.array([0.0005]),0.15,0.45,0.01)
#    dlist=np.append(np.array([0.0005]),np.arange(0.001,0.015,0.001))
#    qlist=np.arange(0.1,1.2,0.1)
#    _,_, _,keff,_ = plotQsignal(cno,False)
#    Deff=keff/1100/2600/0.51
#
#    #adddgrid(cno,dlist,qlist)
#    dnumber,q,Z = getResponse(cno,dlist,qlist,scalar='uds',xmask=slice(0,36))
#    print('D = {} m2/s'.format(dnumber))
#    print('D_eff = {} m2/s'.format(Deff))
    
#for cno in [6,7]:
#    dlist = np.arange(0.007,0.02,0.00025)
#    qlist = np.arange(0.2,1.1,0.05)
#    
#    #adddgrid(cno,dlist,qlist)
#    dnumber,q,Z = getResponse(cno,dlist,qlist,scalar='uds')
#    print('D = {} m2/s'.format(dnumber))
#    print('D_uni = {} m2/s'.format(k_uni(dnumber,q)))
#
#for cno  in [8]:
#    dlist = np.arange(0.015,0.03,0.00025)
#    qlist = np.arange(0.4,1.1,0.05)
#    
#    #adddgrid(cno,dlist,qlist)
#    dnumber,q,Z = getResponse(cno,dlist,qlist,scalar='uds')
#    print('D = {} m2/s'.format(dnumber))
#    print('D_uni = {} m2/s'.format(k_uni(dnumber,q)))
#    
    


# load velocity data for choosen case
caseNo=1
scalar='xvel'

Smat,M,t = loadMdata(caseNo,scalar=scalar)


# Cut away the crazy nozzle. And the nozzle paths. And the space above the bed. 
vinst=Smat[:,6:42,:,:60]


vinst=vinst[:,30,10,20]

vbar = vinst.mean()

vdot = vinst-vbar

def RL(kdot,vinst,vbar,vdot):
    top=0
    ksum=0
    
    if type(kdot) == 'int':
        kmax=kdot
    else:
        kmax=kdot[-1]
    
    for k in range(len(vinst)-kmax):
        top = top + (vinst[k] - vbar)*(vinst[k+kdot] - vbar)
        ksum+=1
    top=top/ksum
    
    return top/np.mean(vdot**2)

karr=np.arange(25)


TL=np.trapz(RL(karr,vinst,vbar,vdot))

DT= np.mean(vdot**2)*TL


plt.plot(t[karr],RL(karr,vinst,vbar,vdot))
plt.title('$R_L$')
plt.xlabel('time [s]')

DTlist=[]
Klist=[]
for K in range(2,40):    
    karr=np.arange(K)
    TL=np.trapz(RL(karr,vinst,vbar,vdot))
    DT=np.mean(vdot**2)*TL
    DTlist.append(DT)
    Klist.append(K)

plt.figure()
plt.plot(Klist,DTlist)
plt.title('$D_T$')
plt.xlabel('time lag [s]')













