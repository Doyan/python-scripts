#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 11:10:35 2018

@author: gabgus
"""

import os
import shutil
import subprocess

# -------------------------------------------------------------------------

def runCase(caseNo,knumber,q):
    
    kdir = '/scratch/gabgus/fluent/kgenerator/'

    btemps = {'high' : [0.0, 1123.15, 1123.15, 1123.15, 1123.15,
                        1123.15, 1123.15, 1123.15, 1198.12],
              'low'  : [0.0, 1073.15, 1073.15, 1073.15, 1073.15,
                        1073.15, 1073.15, 1073.15, 1148.17],
              'start': [0,1098.15,1098.15,1098.15,1098.15,1098.15,
                        1098.15,1098.15,1173.15],
              'dT'   : [0.0, 5.1537, 4.6695, 3.8562, 3.8035,
                        0.0, 6.0, 5.6, 7.36]}


    simspecs = {'sfreq' : [0, 2, 2, 2, 2, 0, 1, 1, 1],
                'imax' : [0,466,466,466,466,0,187,187,176],
                'its' : [0,0.05,0.05,0.05,0.05,0,0.01,0.01,0.05],
                'ts' : [0,0.05,0.05,0.05,0.05,0,0.04,0.04,0.05],
                'lxmin': [0,0.195,0.195,0.120,0.165,0,0.195,0.195,0.126],
                'lxmax': [0,0.345,0.345,0.27,0.315,0,0.345,0.345,0.186],
                'htempx' : [0,0.0151,0.0151,0.0151,0.0151,0,0.0151,0.0151,0.0065],
                'ltempx' : [0,0.524,0.524,0.524,0.524,0,0.524,0.524,0.305],
                'lid' : [0,9,9,9,9,9,9,9,9],
                'hid': [0,16,15,15,15,0,16,15,16],
                'cid': [0,21,19,19,19,0,21,19,20]}
    
    
    
    current=os.getcwd()

    os.chdir(kdir)

    # Make journal file
    # -------------------------------------------------------------------------
    runpath = 'runs/case{}/q{}/k{}/'.format(caseNo,q,knumber)
    
    vof_s=[0 ,0.4994 ,0.4983, 0.4985, 0.4986, 0, 0.5022, 0.5010, 0.3013]


    rho_s = 2613.0
    rho= rho_s*vof_s[caseNo] # kg/m3 bulk density of bed !!!
    
    cp = 1203.0 # J/kgK bed material cp

    yfactor=10 # factor for diffusivity matrix
    yfactor_low = yfactor/q # same vertical diffusivity even though lower lateral
        
    lowk = round(knumber*q,4)
    hightemp = round(btemps['high'][caseNo] - btemps['dT'][caseNo],4)
    lowtemp = round(btemps['low'][caseNo] + btemps['dT'][caseNo],4)
    starttemp = btemps['start'][caseNo]
    
    # temporal
    samplefreq = simspecs['sfreq'][caseNo]
    timestep = simspecs['ts'][caseNo]
    initialtstep = simspecs['its'][caseNo]
    max_iter = simspecs['imax'][caseNo]
    
    # spatial
    lxmin=simspecs['lxmin'][caseNo]
    lxmax=simspecs['lxmax'][caseNo]
    htempx=simspecs['htempx'][caseNo]
    ltempx=simspecs['ltempx'][caseNo]
    lid = simspecs['lid'][caseNo]
    hid = simspecs['hid'][caseNo]
    cid = simspecs['cid'][caseNo] 
    
    with open('active_journal.jou','w') as journal:
     # --------------- setup 
        journal.write('rc inputs/c{}.msh \n'.format(caseNo))
        
        journal.write('define/models/unsteady-1st-order? yes\n')
        journal.write('define/models/energy? yes no no no no\n')
        
        journal.write('define/materials change-create aluminum bed_material yes \n')
        journal.write(', {} yes , {} yes anisotropic 1 0 0 {} , {} ,\n'.format(rho,cp,yfactor,knumber))
        journal.write('define/materials change-create aluminum bed_material2 yes \n')
        journal.write(', {} yes , {} yes anisotropic 1 0 0 {} , {} ,\n'.format(rho,cp,yfactor_low,lowk))
        
        journal.write('adapt/mark-inout-rectangle , , {} {} , ,\n'.format(lxmin,lxmax))
        journal.write('mesh/modify-zones/sep-cell-zone-mark surface_body 0 yes \n')
        journal.write('define/boundary-conditions solid surface_body yes bed_material , , , , , , , , , \n')
        journal.write('define/boundary-conditions solid {} yes bed_material2 , , , , , , , , , \n'.format(lid))
        
        journal.write('mesh/modify-zones/slit-interior-between-diff-solids \n')
        journal.write('adapt/mark-inout-rectangle , , 0 {} , , \n'.format(htempx))
        journal.write('mesh/modify-zones/sep-cell-zone-mark surface_body 1 yes \n')
        journal.write('adapt/mark-inout-rectangle , , {} 1 , , \n'.format(ltempx))
        journal.write('mesh/modify-zones/sep-cell-zone-mark surface_body 2 yes \n')
        
        journal.write('/define/boundary-conditions solid {} no no yes yes {} , , , , , , , \n'.format(hid,hightemp))
        journal.write('/define/boundary-conditions solid {} no no yes yes {} , , , , , , , \n'.format(cid,lowtemp))
        journal.write('/solve/initialize/set-defaults temperature {} \n'.format(starttemp))
        
        
        
        # ------------ solving 
        journal.write('/file/transient-export/ascii \n')
        journal.write('{}temp , temperature q yes yes , {} flow-time ,\n'.format(runpath,samplefreq))
        journal.write('/solve/set time-step {}\n'.format(initialtstep))
        journal.write('/solve/initialize/initialize-flow \n')
        journal.write('/solve/dual-time-iterate 1 50\n')
        journal.write('/solve/set time-step {}\n'.format(timestep))
        journal.write('/solve/dual-time-iterate {} 50\n'.format(max_iter))
        
        # write and exit
        journal.write('wc Tprofile.cas ok\n')
        journal.write('wd Tprofile.dat ok\n') 
        journal.write('exit yes\n')
    # -------------------------------------------------------------------------    
    
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
    logname = 'logs/c{}-q{}-k{}.log'.format(caseNo,q,knumber)
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
    
    