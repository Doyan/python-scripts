#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 15:51:33 2018

@author: gabgus
"""

import os
import shutil
import subprocess



# -------------------------------------------------------------------------

def runCase(caseNo,dnumber,q):
    
    ddir= '/scratch/gabgus/fluent/dgenerator/'
    
    simspecs = {'sfreq' :[0, 2, 2, 2, 2, 0, 1, 1, 1],
                'imax' : [0,466,466,466,466,0,187,187,176],
                'its' :  [0,0.05,0.05,0.05,0.05,0,0.01,0.01,0.05],
                'ts' :   [0,0.05,0.05,0.05,0.05,0,0.04,0.04,0.05],
                'lxmin': [0,0.42,0.42,0.345,0.39,0,0.42,0.42,0.27],
                'lxmax': [0,0.57,0.57,0.495,0.54,0,0.57,0.57,0.33],
                'pmin':  [0,0.495,0.495,0.495,0.495,0,0.495,0.495,0.3]}
        
    current=os.getcwd()
    
    os.chdir(ddir)
    
    # Make journal file
    # -------------------------------------------------------------------------
    vof_s=[0 ,0.4994 ,0.4983, 0.4985, 0.4986, 0, 0.5022, 0.5010, 0.3013]
    
    
    rho_s = 2613.0
    rho= rho_s*vof_s[caseNo] # kg/m3 bulk density of bed 
    
    cp = 1203.0 # J/kgK bed material cp
    
    yfactor=10 # vertical factor for diffusivity matrix
    yfactor_low = yfactor/q # same vertical diffusivity even though lower lateral
    
    D = dnumber*rho # fluent wants mass diffusivity i.e. kg/m/s
    
    lowd = round(D*q,4) # diffusivity lower around wall
    
    runpath = 'runs/case{}/q{}/d{}/'.format(caseNo,q,dnumber)
    
    # temporal
    samplefreq = simspecs['sfreq'][caseNo]
    timestep = simspecs['ts'][caseNo]
    initialtstep = simspecs['its'][caseNo]
    max_iter = simspecs['imax'][caseNo]
    
    
    # spatial
    pmax = 1.0
    pmin = simspecs['pmin'][caseNo]
    lxmin=simspecs['lxmin'][caseNo]
    lxmax=simspecs['lxmax'][caseNo]
    
    
    
    with open('active_journal.jou','w') as journal:
    # ------------------ setup
        journal.write('rc inputs/c{}.msh \n'.format(caseNo))
        
        journal.write('define/models/unsteady-1st-order? yes \n')
        journal.write('define/user-defined/user-defined-scalars 1 no no no yes "none" , \n')
        
        journal.write('define/materials change-create aluminum bed_material yes \n')
        journal.write(', {} yes , {} , yes , 0 anisotropic , , , {} , {} -1 , \n'.format(rho,cp,yfactor,D))
        
        journal.write('define/materials change-create aluminum bed_material2 yes \n')
        journal.write(', {} yes , {} , yes , 0 anisotropic , , , {} , {} -1 , \n'.format(rho,cp,yfactor_low,lowd))
        
        journal.write('adapt/mark-inout-rectangle , , {} {} , , \n'.format(lxmin,lxmax))
        journal.write('mesh/modify-zones/sep-cell-zone-mark surface_body 0 yes \n')
        journal.write('define/boundary-conditions solid surface_body yes bed_material , , , , , , , , , \n')
        journal.write('define/boundary-conditions solid surface_body:009 yes bed_material2 , , , , , , , , , \n')
        
        journal.write('mesh/modify-zones/slit-interior-between-diff-solids \n')
        journal.write('define/boundary-conditions wall interior-surface_body:004-shadow no no yes \n')
        journal.write('define/boundary-conditions wall interior-surface_body:006-shadow no no yes \n')
        
        journal.write('adapt/mark-inout-rectangle , , {} {} , , \n'.format(pmin,pmax))
        journal.write('solve/initialize/initialize-flow \n')
        journal.write('solve/patch , 1 , , uds-0 1.0 \n')
    
    # --------- solving
        
        journal.write('/file/transient-export/ascii \n')
        journal.write('{}conc , uds-0 q yes yes , {} flow-time ,\n'.format(runpath,samplefreq))
        journal.write('/solve/set time-step {}\n'.format(initialtstep))
        journal.write('/solve/dual-time-iterate 1 50\n')
        journal.write('/solve/set time-step {}\n'.format(timestep))
        journal.write('/solve/dual-time-iterate {} 50\n'.format(max_iter))
    
    # ----------- write and exit
        journal.write('wc Cprofile.cas ok\n')
        journal.write('wd Cprofile.dat ok\n') 
        journal.write('exit yes\n')
         
    # -----------------------------------------------------------------------------
    
    try: 
        os.mkdir(runpath)
    
    except FileNotFoundError:
        qpath ='runs/case{}/q{}/'.format(caseNo,q) 
        os.mkdir(qpath)
        os.mkdir(runpath)
        
    except FileExistsError:
        shutil.rmtree(runpath)
        os.mkdir(runpath)
    
    ecode=0
    logname = 'logs/c{}-q{}-d{}.log'.format(caseNo,q,dnumber)
    with open(logname,'w') as f:
        print('\n running fluent ... \n check logfile {} for details. \n please wait ...'.format(logname))
        cmd = ['fluent','2ddp','-g','-t2', '-i', 'active_journal.jou']
        ecode = subprocess.call(cmd,stdout=f)
        if ecode == 0:
            print('\n\n Fluent exited successfully')
        else:
            print('\n\n Error in fluent run. \n check logfile and redo')
    
    os.chdir(current)
    return ecode









