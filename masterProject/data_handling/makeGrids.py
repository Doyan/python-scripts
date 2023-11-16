#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 15:29:01 2019

@author: gabgus
"""

import numpy as np, matplotlib.pyplot as plt
import comparator as C

# -----------------------------------------------------------------------------


qlist = np.arange(0.1,1.2,0.05)


for caseNo in [1,2,3,4,6,7,8]:
    
    if caseNo < 5:
        dlist = np.round(np.arange(0.0001,0.007,0.0001),6)
    elif caseNo < 8:
        dlist = np.round(np.arange(0.008,0.013,0.0001),6)
    else:
        dlist = np.round(np.arange(0.02,0.03,0.0001),6) 
    
    
    klist = np.round(dlist * 1200 * 2600 * C.vof_s[caseNo],1)
    
    print(f'========  Case: {caseNo}  ===========\n')
    
    #C.adddgrid(caseNo,dlist,qlist,dim1=True)
    #C.addkgrid(caseNo,klist,qlist,dim1=True)
    C.getResponse(caseNo,dlist,qlist,tslice=slice(0,240),scalar='uds',xmask=slice(0,500))
    

print('\n Done!! \n')