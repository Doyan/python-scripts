#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 14:49:14 2019

@author: gabgus
"""
#%%
import numpy as np, matplotlib.pyplot as plt
import comparator as c

# -----------------------------------------------------------------------------
#%%

scalar = 'vof'
sNo=51



for caseNo in [9,10]:
    plt.figure()
    Smat, M, t = c.loadMdata(caseNo,scalar)
    
    c.plot2d(caseNo,M[0][:,:,1],M[1][:,:,1],zbar=Smat[sNo][:,:,46],time=t[sNo],scalar=scalar,cmap='RdBu_r')
    plt.title('Instantaneous volume fraction')
    plt.savefig('paper/characterisation/vof_inst_c{}.pdf'.format(caseNo),bbox_inches='tight')


    




for caseNo in [1]:
    plt.figure()
    Smat, M, t = c.loadMdata(caseNo,scalar)
    
    c.plot2d(caseNo,M[0][6:73,:,1],M[1][6:73,:,1],zbar=Smat[sNo][6:73,:,50],time=t[sNo],scalar=scalar,cmap='RdBu_r')
    plt.title('Instantaneous volume fraction')
    plt.savefig('paper/characterisation/vof_inst_c{}.pdf'.format(caseNo),bbox_inches='tight')















