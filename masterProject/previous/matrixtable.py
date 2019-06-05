#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 29 09:13:36 2018

@author: gabgus
"""


import numpy as np, matplotlib.pyplot as plt
import pandas as pd


fnam ='/scratch/gabgus/geometry/present/images_05-10-2018/result.csv'


D=pd.read_csv(fnam,sep=',',)

dTavg = D.T_C_avg - D.T_G_avg
dTmax = D.T_C_max - D.T_G_min

table = np.column_stack((D.case,dTavg,dTmax))

print(table)