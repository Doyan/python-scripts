#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 14:40:41 2019

@author: gabgus
"""


import numpy as np
import matplotlib.pyplot as plt
import comparator as C

scalar='temp'

for caseNo in [1,6,9,8,10]:
    ybar,gx,t = C.get1ddata(caseNo,scalar)
    
    xmask = slice(2,-3)
    
    
    dlist= np.arange(0.0001,0.026,0.00005)
    
    
    tint = 10
    
    if caseNo == 9:
        tint = 50
    
    tend = int(np.floor(t.shape[0] / tint) * tint)
    
    tmaxs = np.arange(tint*1,tend+tint,tint)
    print(tmaxs[-1],len(t))
    #tmaxs = [tmaxs[-1]]
    errlist = []
    dbestlist = []
    for tmax in tmaxs:
        tslice = slice(1,int(tmax))
        err = np.empty((len(dlist),))
        
        for i,knumber in enumerate(dlist):
            knumber = knumber*1200*2600*C.vof_s[caseNo]
            errtot = C.getHomoDiff(caseNo,ybar,gx,t,knumber,scalar=scalar,
                                 tslice=tslice,
                                 xmask=xmask,
                                 tol=0.0001,
                                 Tcorr=0)
            err[i]=errtot
    
        imin=np.argmin(err)
        idx=np.unravel_index(imin,err.shape)[0]
        errlist.append(err[idx])
        dbestlist.append(dlist[idx])
        
    plt.grid()
    plt.scatter(tmaxs*0.05,dbestlist,c=errlist)
    plt.colorbar()
    plt.ylim(0,max(dbestlist)*1.1)
    plt.show()

    print(dlist[idx])

plt.plot(dlist,err)
   
plt.plot(dlist[idx],err[idx],'ro')