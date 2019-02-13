#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  2 11:09:10 2018
Script to correctly sort timestep-marked images and make into gif



@author: gabgus
"""
import os
import subprocess
# -------------------------------------------------------------------------

for caseNo in range(8):
    for contNo in [1,2,3,4]:
        namelist = ['error','vof','tgas','cgas','ggas']


        path = '/scratch/gabgus/fluent/kluster/case{}/images/c{}/'.format(caseNo,contNo)
        outpath = '/scratch/gabgus/fluent/kluster/results/c{}_{}.gif'.format(caseNo,namelist[contNo])

        current = os.getcwd()
        os.chdir(path)
        print('calling convert for c{}_{}.gif \n please wait...'.format(caseNo,namelist[contNo]))
        ecode = subprocess.call('convert -delay 10 $(ls -v *.png) {}'.format(outpath) ,shell=True)
        print(ecode)