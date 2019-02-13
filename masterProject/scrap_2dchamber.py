#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 18 10:18:26 2018

Script for doing quick and dirty plots on the 2dchamber case

@author: gabgus
"""

import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import os
import shutil
import subprocess

# ----------------------------------------------------------------------------
datapath = '/scratch/gabgus/fluent/kluster/case1/datafiles/'
casepath = datapath.split('datafiles/')[0]
casefile = datapath + 'convective3d-6.cas'

def pullImages(datapath,casepath,casefile,overwrite=False):
    filelist= [file for file in os.listdir(datapath) if file.endswith('0.dat')]
    times = [filename.split('-')[2].split('.dat')[0] for filename in filelist]
    
    datalist=list(sorted(zip(times,filelist)))
          
    current = os.getcwd()
    os.chdir(casepath)
    try:
        os.mkdir('images')
    except FileExistsError:
        if not overwrite:
            ans = input('Image folder exists, redo it? (y/n)')
        
        if (ans == 'y') or overwrite:
            shutil.rmtree(casepath + 'images')
            os.mkdir('images')
        else:
            print('No action taken, error')
        
    for contour in range(1,5):
        os.mkdir('images/c{}'.format(contour))
        
    with open('myJournal.jou','w') as jou:
        jou.write('rc {}\n'.format(casefile))
        for time, filename in datalist:
            jou.write('rd {}\n'.format(datapath + filename))
            jou.write('/display/views/restore-view view-0\n')
            jou.write('/display/objects/display contour-1\n')
            jou.write('/display/save-picture ./images/c1/vof_at-{}.png\n'.format(time))
            jou.write('/display/objects/display contour-2\n')
            jou.write('/display/save-picture ./images/c2/tgas_at-{}.png\n'.format(time))
            jou.write('/display/objects/display contour-3\n')
            jou.write('/display/save-picture ./images/c3/cgas_at-{}.png\n'.format(time))
            jou.write('/display/objects/display contour-4\n')
            jou.write('/display/save-picture ./images/c4/ggas_at-{}.png\n'.format(time))
        jou.write('exit\n')
        jou.write('yes\n')
    
    with open('printing.log','w') as f:
        cmd = ['fluent','3ddp', '-t2', '-gu', '-i', 'myJournal.jou']
        ecode = subprocess.call(cmd,stdout=f)
    
    os.chdir(current)
    return ecode    
    
    
#
#
#for case in range(8):
#    datapath = '/scratch/gabgus/fluent/kluster/case{}/datafiles/'.format(case)
#    casepath = datapath.split('datafiles/')[0]
#    casefile = datapath + 'convective3d-6.cas'
#    print('starting on case{}: \n please wait...'.format(case))
#    pullImages(datapath,casepath,casefile,overwrite=True)
    


    















#%%

# Globals

dnstart = [50e3, 125e3, 125e3, 125e3, 50e3, 125e3, 125e3, 125e3]
dnstop = [118e3, 295e3, 295e3, 295e3, 118e3, 295e3, 295e3, 295e3]

# -------------------------------------------------------------------------
# Functions


def getmtotal(case):
    total_mass = pd.read_csv(casepath + 'total_mass.out',skiprows=[0,1,2],names=['tstep','time','mass'],delimiter=' ')
    return total_mass


def getmfr(case,interface):        
    nstart=int(dnstart[case])
    nstop=int(dnstop[case])

    casepath = '/scratch/gabgus/fluent/kluster/case{}/'.format(case)
    mfr = pd.read_csv(casepath + 'mfr.out',skiprows=[0,1,2],names=['tstep','time','mfr_periodic','mfr_gap','mfr_tin','mfr_tout'],delimiter=' ')     
    mfr['mfr'] = mfr[interface] 
    mfrslice=mfr.mfr[nstart:nstop]
    mmean=np.array(mfrslice).mean()

    return mfr, mmean 

def getpin(case):
    nstart=int(dnstart[case])
    nstop=int(dnstop[case])
    
    casepath = '/scratch/gabgus/fluent/kluster/case{}/'.format(case)
    pin = pd.read_csv(casepath + 'p_in.out',skiprows=[0,1,2],names=['tstep','time','p'],delimiter=' ')    
    pmean=pin.p[nstart:nstop].sum()/(nstop-nstart)
    return pin, pmean
    


def plotmfr(case, interface,shift=0.0):
    mfr,mmean = getmfr(case,interface)
    
    ax = plt.plot(mfr.time + shift ,mfr.mfr, label= 'mfr over {}'.format(interface))
    return ax

def plotmean(case,interface):
    nstart=int(dnstart[case])
    nstop=int(dnstop[case])
    
    mfr,mmean = getmfr(case,interface)
    
    plt.plot([mfr.time[nstart],mfr.time[nstop]],[mmean,mmean],'r--',label='mfr_mean = {:.4} kg/s'.format(mmean))
    return

def plotcanvasmfr():
    plt.xlabel('Time [s]')
    plt.ylabel('Massflow [kg/s]')
    plt.title('Massflow over selected interface')
    plt.xlim(4.5,7.9)
    plt.ylim(-15.0,0.8)
    plt.legend(loc='best')
    return

def plotpin(case):
    nstart=int(dnstart[case])
    nstop=int(dnstop[case])
    
    pin,pmean = getpin(case)
    plt.plot(pin.time,pin.p,label='pressure-drop signal')    
    plt.plot([pin.time[nstart],pin.time[nstop]],[pmean,pmean],'r--',label='dP_mean = {:.4} Pa'.format(pmean))
    plt.xlim(pin.time[nstart],pin.time[nstop])

    return 

def plotcanvaspin():
    plt.ylim(3.5e3,8e3)
    plt.legend(loc='best')
    plt.title('Pressure drop, measured pressure at bottom of bed')
    plt.xlabel('Time [s]')
    plt.ylabel('Pressure [Pa]')
    return

# -------------------------------------------------------------------------
#  Main
    

case = 2

interface = 'mfr_gap'
plt.figure(1)
plotmfr(case,interface)
plotmfr(case,'mfr_tout',shift=0.0)
plotmean(case,interface)
#plotmean(case,'mfr_tout')
plotcanvasmfr()
plt.savefig('/scratch/gabgus/fluent/kluster/results/images/dynsand.pdf',bbox_inches='tight')
plt.show()    
#plt.plot(total_mass.time,total_mass.mass,label='Total mass of sand in system')
#plt.ylim(28,30.5)
#plt.xlabel('Time [s]')
#plt.ylabel('Total mass [kg]')
#plt.title('Total mass in system')
#plt.legend(loc='best')

#%%





#sint=0.01
#dt=np.array([5e-5,2e-5,2e-5,2e-5])
#sfreq=sint/dt
#
#dt = np.tile(dt,2)
#sfreq = np.tile(sfreq,2)
#
#vel 
#x1 = 
#x2 





def makeJournal(tstep,freq,casfile,datfile):
    with open('myJournal.jou','w') as journal:
        journal.write('rc {}\n'.format(casfile))
        journal.write('rd {}\n'.format(datfile))
        journal.write('/solve/report-files/edit/mfr-rset file-name mfr.out\nq\n')
        journal.write('/solve/report-files/edit/p-in-rset file-name p_in.out\nq\n')
        journal.write('/solve/report-files/edit/total_mass-rset file-name total_mass.out\nq\n')
        journal.write('/file/auto-save root-name ./datafiles/convective3d\n')
        journal.write('/solve/set timestep {}\n'.format(tstep))
        journal.write('/file/auto-save data-frequency {}\n'.format(freq))
        journal.write('/solve dual-time-iterate 300000 40\n')
    return os.getcwd() + '/myJournal.jou'

def editCfile(cfile,outfile,vel,x1,x2):
    sed1='s/#define RISER_VELOCITY 3.0/#define RISER_VELOCITY {}/g'.format(vel)
    sed2='s/#define DIVIDE_MIDDLE 1.290/#define DIVIDE_MIDDLE {}/g'.format(x1)
    sed3='s/#define DIVIDE_HIGH 2.580/#define DIVIDE_HIGH {}/g'.format(x2)
    
    os.popen('sed -e \'{}\' -e \'{}\' -e \'{}\' {} > {}'.format(sed1,sed2,sed3,cfile,outfile))
    return 



def makeKlusterfolder(dt,sfreq,vel=[3.0, 5.0, 7.0, 9.0]*2,x1=[1.290]*4 + [1.065]*4,x2= [2.580]*4 + [2.805]*4):
    prepdir='/scratch/gabgus/fluent/kluster/'
    casfile='convective3d.cas'
    datfile='convective3d.dat'
    
    
    current=os.getcwd()
    os.chdir(prepdir)
    for i in range(8):
    
        casedir = prepdir + 'case{}'.format(i)
    
        os.mkdir(casedir)
        os.popen('cp {} {}'.format(prepdir + casfile, casedir + '/' + casfile))
        os.popen('cp {} {}'.format(prepdir + datfile, casedir + '/' + datfile))
        
        editCfile(prepdir+'convective.c',casedir+'/convective.c',vel[i],x1[i],x2[i])
                    
        os.chdir(casedir)
        os.mkdir('datafiles')
        makeJournal(casfile,datfile,dt[i],int(round(sfreq[i])))
        os.chdir(prepdir)

    os.chdir(current)
    return


#%%











