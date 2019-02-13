#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 13:35:58 2018

@author: gabgus
"""

import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

import os
import shutil
import subprocess

# -------------------------------------------------------------------------
datapath = "/media/TOSHIBA EXT/BioShare/pp/noWall"
exportpath = '/scratch/gabgus/mfdata_all/pp_real/noWall'

casfile='/scratch/gabgus/mfdata_all/setteValidation-2-copy.cas'

def pullData(datapath,casfile,exportpath,overwrite=False):
    
    try:
        os.mkdir(exportpath)
    except FileExistsError:
        ans=0
        if not overwrite:
            ans = input('export data folder exists, redo it? (y/n): ')
        
        if (ans == 'y') or overwrite:
            shutil.rmtree(exportpath)
            os.mkdir(exportpath)
        else:
            print('No action taken, error')
    exceptions=set(['pp_wWall_getReal.dat','bioshare_pp_wWall.dat'])
    filelist= [file for file in os.listdir(datapath) if file.endswith('.dat') and file not in exceptions ]
    times = [filename.split('-')[2].split('.dat')[0] for filename in filelist]
        
    datalist=list(sorted(zip(times,filelist)))
    
    current=os.getcwd()
    
    os.chdir(datapath)        
    
    with open('myJournal.jou','w') as journal:
        journal.write('rc {}\n'.format(casfile))
        for time,datfile in datalist:
            journal.write('rd {}\n'.format(datfile))
            
            exportname= '{}dat_{}.csv'.format(exportpath,time)
            
            journal.write('/file/export ascii {} , yes\n'.format(exportname))
            #journal.write('udm-1\n')
            #journal.write('udm-2\n')
            journal.write('phase-2-temperature\n')
            journal.write('phase-2-vof\n')
            journal.write('phase-2-x-velocity\n')
            journal.write('phase-2-y-velocity\n')
            journal.write('phase-2-z-velocity\n')
            journal.write('phase-2-uds-0-scalar\n')
            journal.write('q\nyes\n')
            
        journal.write('exit\nyes\n')
    ecode=0
    with open('printing.log','w') as f:
        print('\n running fluent ... \n check logfile "printing.log" for details. \n please wait ...')
        cmd = ['fluent','3ddp','-g', '-i', 'myJournal.jou']
        ecode = subprocess.call(cmd,stdout=f)
        if ecode == 0:
            print('\n\n Fluent exited successfully! \n Export done.')
        else:
            print('\n\n Error in fluent run. \n check logfile and redo export.')
    os.chdir(current)
    return ecode

dbase="/media/TOSHIBA EXT/BioShare/"
dpaths=[dbase + 'dysorDataSave/case01/',
        dbase + 'dysorDataSave/case02/',
        dbase + 'dysorDataSave/case03/',
        dbase + 'dysorDataSave/case04/',
        dbase + 'pp/wWall/']

ebase="/scratch/gabgus/mfdata_all/"
epaths=[ebase + 'dysor/case01/',
        ebase + 'dysor/case02/',
        ebase + 'dysor/case03/',
        ebase + 'dysor/case04/',
        ebase + 'pp_real/wWall/']

casfiles=['dysor_case1_extract.cas',
          'dysor_case2_extract.cas',
          'dysor_case3_extract.cas',
          'dysor_case4_extract.cas',
          'pp_wWall_getReal_unloaded.cas']


def packQsignal(fileindex,dcell,prefix):
    Q1=[]
    Q2=[]
    t=[]
    for i,(time,file) in enumerate(fileindex):
        data=pd.read_csv(file)
        data.columns=[name.strip() for name in data.columns]
        Q1.append(data['udm-1'].sum()*(dcell**3))
        Q2.append(data['udm-2'].sum()*(dcell**3))
        t.append(time)

    np.save('/scratch/gabgus/mfdata_all/meshed/{}_qin'.format(prefix),np.array(Q1))
    np.save('/scratch/gabgus/mfdata_all/meshed/{}_qut'.format(prefix),np.array(Q2))
    np.save('/scratch/gabgus/mfdata_all/meshed/{}_qtime'.format(prefix),np.array(t))
    return


for dpath,epath,cas in zip(dpaths,epaths,casfiles):
    pullData(dpath,cas,epath,overwrite=True)
    