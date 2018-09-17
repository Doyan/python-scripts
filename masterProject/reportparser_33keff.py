#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  6 14:46:23 2018

@author: gabgus
"""

import pandas as pd
import os

# -------------------------------------------------------------------------

df_list = []

#reportdir='/scratch/gabgus/geometry/redone/kluster/reports/'
#reportdir='/scratch/gabgus/geometry/redone/results/images_07-13-2018/reports/'
reportdir='C:/Users/Gabriel/simple-model/reports/'

filelist=os.listdir(reportdir)


df0 = pd.read_csv(reportdir + filelist[0],index_col=0,header =0, usecols=[0,1,2])
df0.case=int(filelist[0].split(sep='_')[1].split(sep='report')[0])
#df0.case=filelist[1].split(sep='model-')[1].split(sep='.')[0]
df0=df0.T

values = df0.iloc[0].copy()
units = df0.iloc[1].copy()

df0.iloc[0]= units
df0.iloc[1] = values
df0=df0.drop(labels=' Unit',axis=0)

unitframe=df0

for file in filelist:
    casestring=file.split(sep='_')[1].split(sep='report')[0]
    #casestring=file.split(sep='model-')[1].split(sep='.')[0]
    df = pd.read_csv(reportdir + file,index_col=0,header =0, usecols=[0,1], names=['Report', 'Value'])
    df=df.T
    df.case=casestring
    df_list.append(df.T)
    df=df.T
    df0=df0.append(df.Value,ignore_index=True)

df0.case.iloc[1] = str(int(df0.case.iloc[1]))
#df0.drop(df0.index[0])


df0.case = pd.to_numeric(df0.case,errors='coerce')

df0=df0.sort_values(by='case')

df0=df0.set_index('case')

df0.to_csv('keff33.csv')
#%%
cols = df0.columns.tolist()
cols_ordered = [cols[16]] + cols[6:12] + cols[-3:-1] + cols[:6] + cols[12:16]  + cols[17:-3] +[cols[-1]] 

df0=df0[cols_ordered]
df0=df0.set_index('case')

df0=df0.drop_duplicates()

df0.to_excel('testU.xlsx')
