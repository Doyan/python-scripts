# -*- coding: utf-8 -*-
"""
Created on Wed Sep  5 17:01:16 2018

Script to plot solutions to a 1 dimensional heat equation 
solved analytically through separation of variables.

Corresponds to conduction through a rod of length L with its ends kept at 
fixed temperatures Th and Tc, initially haviing a uniform temperature of T0

Necessary material properties are k, rho and cp

# Th
#    .
#        .
#            .
#                Tc
# 
# |------ L ------|

@author: Gabriel
"""
import numpy as np, matplotlib.pyplot as plt

import pandas as pd
import os

# -----------------------------------------------------------------------------
# Inputs

Th = 845 + 273.15   # K
Tc = 805 + 273.15   # K
T0 = 825 + 273.15   # K

L = 0.51    # m

k = 20  # W/mK

t = 23.6


# -------------------------------------
# Properties

rho_s = 2610.0  # kg / m3 
rho_g = 1.0     # kg / m3

cp_g = 1.17e3    # J / kgK
cp_s = 1.23e3     # J / kgK

eps = 0.45  # void fraction

rho = eps*rho_g + (1 - eps)*rho_s # kg / m3
cp = eps*cp_g + (1-eps)*cp_s # J / kgK

def S(x,L):
    a = (Tc -Th)/L
    b = Th
    return a*x + b
    
def v(x,t,alpha,L):
    N = int(2000*L)
    S = 0 
    pi = np.pi
    for k in range(N):
        n = k + 1
        psum = L/(n*pi)*((T0 - Th) - (T0 - Tc)*(-1)**n)*np.sin(n*pi/L*x)*np.exp(-n**2*pi**2/L**2*alpha*t)
        S += psum
    return 2/L * S

def u(x,t,k,L=L):
    alpha = k/(rho*cp)    
    return v(x,t,alpha,L) + S(x,L)

N = 200
x = np.linspace(0,L,N)

T = u(x,t,k)

ax=plt.figure(1,[6,7])

plt.xlabel('L [m]')
plt.ylabel('T [$\degree C$]')
plt.title('Solution to the 1D heat equation at t = ' +str(t) + 's')
plt.plot(x,T-273.15)
plt.plot([0,L],[Th-273.15,Tc-273.15],'k--')
plt.grid()

#%%
# Calculation of accumulated deviation from simulated solution
# for array of k-values

caseNo=2

#datapath = './datafiles/case0' + str(caseNo) + '/'
datapath = ''
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

# # read folder and get lists with dataframes and length-parameters sorted after time 
def getDflist(datapath=datapath,tscaling=tscaling,toffset=toffset):
    filelist = os.listdir(datapath)
    
    csvlist =  [f for f in filelist if '.txt' in f]
    times = [ csvfile.split('T_')[1] for csvfile in csvlist]
    times = [ round(float(time.split('.')[0]) * tscaling + toffset,4) for time in times]
    
    datalist=list(sorted(zip(times,csvlist)))
    
    noSamples = len(datalist)
    
    dflist = []
    Llist = []
    x0list = []
    dcelllist = []
    times = []
    
    for sampleNo in range(noSamples):
        df, L, x0, dcell = getDf(sampleNo,datalist)
        dflist.append(df)
        Llist.append(L)
        x0list.append(x0)
        dcelllist.append(dcell)
        times.append(datalist[sampleNo][0])
        
    return times,dflist,Llist,x0list,dcelllist

# compute the deviation from the analytical solution at time t along x
def getXdev(df,time,k,L):
    Nx = len(df)

    xDev = []
    for i in range(Nx):
        xDev.append(abs(df.Ti[i] - u(df.x[i],time,k,L)))
    return sum(xDev) / Nx

# compute the accumulated deviation from the analytical solution along time
def getTimedev(dflist,times,L,k):
    Ntimes = len(times)
    tDev = []    
    for i in range(Ntimes):        
        df = dflist[i]
        time= times[i]
        
        tDev.append(getXdev(df,time,k,L))
    return sum(tDev) / Ntimes, tDev

        
# plot profile from csv together with analytic solution for a guessed k
def plotKguess(sampleNo,k):    
    times,dflist,Llist,x0list,dcelllist = getDflist()
    
    df = dflist[sampleNo]
    L = Llist[sampleNo] - dcelllist[sampleNo]
    time = times[sampleNo]
    
    N = 200
    x = np.linspace(0,L,N)

    T = u(x,time,k,L)

    plt.xlabel('L [m]')
    plt.ylabel('T [$\degree C$]')
    plt.title('Solution to the 1D heat equation at t = ' +str(time) + 's')

    plt.plot(df.x,df.Ti-273.15,'r.')
    plt.plot(x,T-273.15)
    plt.plot([0,L],[Th-273.15,Tc-273.15],'k--')
    plt.grid()
    plt.show()



#times,dflist,Llist,x0list,dcelllist = getDflist()
#
#L = Llist[0] -dcelllist[0]   
#
#plotKguess(229,25000)
#print(dcelllist[0])
### array of k-values to try
#krange = [2389, 2390, 2390.5] 
#
##  loop to try each value and settle for the one with least deviation
#ktrue = 0
#minAck = 100000
#acklist = []
#for k in krange:
#    ack, tDev = getTimedev(dflist,times,L,k)
#    acklist.append(ack)
#    if ack < minAck:
#        minAck=ack
#        ktrue=k
#
#print(ktrue)    


#n=0
#df=dflist[n]
#time= times[n]
#
#klist = []
#for i in range(12):
#    Y = (df.Ti[i] - T0)/(Th - T0)
#    klist.append(rho*cp/(12*(1-Y))*df.x[i]**2/time)



