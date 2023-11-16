#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 16:06:40 2018

@author: gabgus
"""
import os
import pandas as pd, numpy as np

# -------------------------------------------------------------------------

caseNo=1

datapath = '/scratch/gabgus/mfdata_all/dysor/case0{}/inst/'.format(caseNo)
iterpath = '/scratch/gabgus/mfdata_all/dysor/case0{}/iterlist.txt'.format(caseNo)

def renamedata():
    iterlist = pd.read_csv(iterpath,delimiter=' ',names=['iteration','time'],)
    iterlist = iterlist.set_index('iteration')

    filelist = os.listdir(datapath)
    iterations = [int(file.split('_')[1].split('.txt')[0]) for file in filelist]

    fileindex = list(sorted(zip(iterations,filelist)))


    current = os.getcwd()
    os.chdir(datapath)
    for iteration,filename in fileindex:
        time = np.array(iterlist.time.loc[iteration])[0]
        os.rename(filename,'dat_{}.csv'.format(time))
    
    os.chdir(current)
    return

def addheader():
    for caseNo in [2,3,4]:
        datapath = '/scratch/gabgus/mfdata_all/dysor/case0{}/inst/'.format(caseNo)
    
        current = os.getcwd()
        filelist = [file for file in os.listdir(datapath) if file.endswith('.csv')]
        N=len(filelist)
        os.chdir(datapath)
        print('\nrunning...')
        for i,filename in enumerate(filelist):
            data=pd.read_csv(filename,delimiter=' ',
                             names =['x-coordinate', 'y-coordinate','z-coordinate', 
                                     'phase-2-temperature','phase-2-vof','phase-2-x-velocity',
                                     'phase-2-y-velocity','phase-2-z-velocity','phase-2-uds-0-scalar'])
            print('{} of {} files'.format(i,N))
        
            data.to_csv(filename,index=False)
        
        os.chdir(current)
        print('done')
    return