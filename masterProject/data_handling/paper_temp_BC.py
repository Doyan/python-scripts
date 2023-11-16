#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 17:01:10 2019
Check for broken BC at frozen temperature edge
@author: gabgus
"""

import numpy as np
import matplotlib.pyplot as plt
import comparator as C


caseNo= 8
scalar= 'temp'

Smat,M,t = C.loadMdata(caseNo,scalar)

sNo=2
xcoord = 25

x2d = M[2][:,xcoord,:]
y2d = M[1][:,xcoord,:]
zbar = Smat[sNo,:,xcoord,:]


mtemp=C.mtemp

norm=C.MidpointNormalize(mtemp[caseNo]-25,mtemp[caseNo]+25, midpoint=1173) 


plt.contourf(x2d,y2d,zbar,30,norm=norm,cmap='RdBu_r')
plt.colorbar()
plt.show()

caseNo= 8
scalar= 'xvel'

Smat,M,t = C.loadMdata(caseNo,scalar)
Vmat,M,T = C.loadMdata(caseNo,'vof')
sNo=2
xcoord = 25

x2d = M[2][:,xcoord,:]
y2d = M[1][:,xcoord,:]
zbar = Smat[sNo,:,xcoord,:]*Vmat[sNo,:,xcoord,:]

C.plot2d(caseNo,x2d,y2d,zbar,t[sNo],scalar=scalar)
plt.show()

caseNo= 8
scalar= 'vof'

Smat,M,t = C.loadMdata(caseNo,scalar)

sNo=2
xcoord = 25

x2d = M[2][:,xcoord,:]
y2d = M[1][:,xcoord,:]
zbar = Smat[sNo,:,xcoord,:]

C.plot2d(caseNo,x2d,y2d,zbar,t[sNo],scalar=scalar)
plt.show()