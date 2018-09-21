# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 16:43:30 2018

k_fitter -

Compares solutions of a heat transfer problem with a certain k to 
a corresponding temperature field from a multiphase simulation 



@author: Gabriel
"""

import numpy as np, matplotlib.pyplot as plt

import pandas as pd
import os

# -------------------------------------------------------
datapath = './datafiles/k2000/'

# value to scale integer from filename with
tscaling = 0.05
toffset = 0

filelist = os.listdir(datapath)

times = [ float(file.split('Data-')[1]) for file in filelist]
times = np.array(times) * tscaling + toffset   


filename=filelist[0]
df = pd.read_csv(datapath+filename)
df.columns = df.columns.str.strip().str.lower().str.replace('-', '_').str.replace('.', '')
df=df.sort_values(by='x_coordinate')