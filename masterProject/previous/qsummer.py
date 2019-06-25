#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 14:15:52 2018

Script to sum up cell values for a total added heat flux.
 Works on ascii files from fluent


@author: gabgus
"""

import numpy as np, pandas as pd

import os

# _----------------------------------------------------------------------------


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




mfdatapath = '/scratch/gabgus/mfdata_all/'
folders = ['pp_real/noWall/udm/','pp_real/wWall/udm/','sette/udm/']

dcells = [0.015, 0.015, 0.006]
prefixes = ['c6','c7','c8']

for folder,prefix,dcell in zip(folders,prefixes,dcells):
    datafolder = mfdatapath + folder

    delimiter = '_'
    ndelimiters=1
    ending ='.csv'

    filelist = [datafolder + file for file in os.listdir(datafolder) if file.endswith(ending)]
    times = [float(filename.split(delimiter)[ndelimiters].split(ending)[0]) for filename in os.listdir(datafolder)]

    fileindex = list(sorted(zip(times,filelist)))

    packQsignal(fileindex,dcell,prefix)

