#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 14:00:45 2019
Paste into comparator script
@author: gabgus
"""
import matplotlib.pyplot as plt

import numpy as np
import comparator as C

scalar = 'temp'
print('\nD - fit')

#xslice = slice(x0,x1)

dbestlist=[]
for caseNo in [1,6,9,8,10]:
    ybar,gx,t = C.get1ddata(caseNo,scalar)
    
    if scalar == 'temp':    
        if caseNo < 5:
            xmask=slice(2,-2)
        else:
            xmask = slice(2,-2)
    else:
        xmask= slice(0,500)
        
    
    dlist= np.arange(0.0001,0.026,0.00005)

    tslice = slice(1,int(len(t)))
    err = np.empty((len(dlist),))
    for i,knumber in enumerate(dlist):
        knumber = knumber*1200*2600*C.vof_s[caseNo]
        errtot = C.getHomoDiff(caseNo,ybar,gx,t,knumber,scalar=scalar,
                             tslice=tslice,
                             xmask=xmask,
                             Tcorr=0,
                             tol=0.0005)
        
        err[i]=errtot

    imin=np.argmin(err)
    idx=np.unravel_index(imin,err.shape)[0]
    
    print('Case {}, D: {}'.format(caseNo,dlist[idx]))

    plt.plot(dlist,err,label='case {}'.format(caseNo))
   
    plt.plot(dlist[idx],err[idx],'ro')
    dbestlist.append(dlist[idx])

plt.xlim(0,0.02)
plt.legend()
plt.show()