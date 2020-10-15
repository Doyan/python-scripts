# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 22:42:37 2020

@author: Gabriel
"""

import numpy as np

import scipy.stats as stats
import matplotlib.pyplot as plt

# =============================================================================

y = np.array([[
    [130,155,74,180],
    [34,40,80,75],
    [20,70,82,58]],
    
    [[150,188,159,126],
    [136,122,106,115],
    [25,70,58,45]],
    
    [[138,110,168,160],
     [174,120,150,139],
     [96,104,82,60]]
    
    ],dtype=float)


def anova(y,red_df=0):
    a,b,n = y.shape
    
    index = np.array(['A','B','AB','Error','Total'])
    df = np.array([a-1,b-1,(a-1)*(b-1),a*b*(n-1)-red_df,a*b*n-1-red_df])
    SS = np.zeros_like(df,dtype=float)
    
    
    yddd = np.sum(y,axis=(0,1,2))
    yidd = np.sum(y,axis=(1,2))
    ydjd = np.sum(y,axis=(0,2))
    yijd = np.sum(y,axis=2)
    
    
    SS[-1] = np.sum(y**2,axis=(0,1,2)) - yddd**2/(a*b*n)
    SS[0] = 1/(b*n)*np.sum(yidd**2) - yddd**2/(a*b*n)
    SS[1] = 1/(a*n)*np.sum(ydjd**2) - yddd**2/(a*b*n)
    SS_ST = 1/n*np.sum(yijd**2,axis=(0,1)) - yddd**2/(a*b*n)
    SS[2] = SS_ST - SS[0] - SS[1]
    SS[3] = SS[-1] - sum(SS[:-1])
    
    MS = SS/df
    MS[-1] = np.nan
    F0 = MS/MS[-2]
    F0[-2] = np.nan
    
    P0 = np.empty_like(F0)
    for i,f in enumerate(F0):
        if not np.isnan(f):
            P0[i] = stats.f.sf(f,df[i],df[-2])
        else:
            P0[i] = np.nan
    return SS,df,MS,F0,P0,index
    

def print_anova(anova_res,width=(10,4,10,10,10),decimals=(2,0,2,4,9)):
    SS,df,MS,F0,P0,index = anova_res
    
    arrays = [SS,df,MS,F0,P0]
    
    labels = ['SS', 'df','MS','F0','P0']
    
    header = f'\n{"Source":<8}\t'
    for j,label in enumerate(labels):
        header += f'{label:<{width[j]}}\t'
    
    print(header)
    for i,name in enumerate(index):
        s = f'{name:<8}\t'
        for j,array in enumerate(arrays):
            s += f'{round(array[i],decimals[j]):<{width[j]}}\t'
        print(s)


res = anova(y)
print_anova(res)

y_missing = np.array([[3,5,10],
              [4,6,11],
              [np.nan,6,2],
              [7,7,5],
              [0,-4,np.nan],
              [1,2,-6]],dtype=float) 


def estimate_all_missing_values(y,tol=1e-6):
    missing = np.argwhere(np.isnan(y))
    
    y_est = y.copy()
    n_missing= len(missing)
    
    yddbar = np.average(np.nanmean(y))
    
    estimations = [[0.,0.]]
    check = 100
    while check > tol:
        for k,pair in enumerate(missing):
            
            estimation = np.zeros(n_missing)
            i,j = (pair[0], pair[1])
            
            y_est = estimate_missing_value(y_est, i, j)
            estimation[k] = y_est[i,j]
        
        estimations.append(estimation)
        check = sum(np.abs(estimations[-1] - estimations[-2])/yddbar)
        print('running')
    
    df_error_reduced = n_missing
    return y_est, df_error_reduced


def estimate_missing_value(y,i,j):
    y[i,j] = 0
    a,b = y.shape
    
    ydd = np.sum(np.nansum(y))
    yid = np.nansum(y,axis=1)
    ydj = np.nansum(y,axis=0)
    
    value = (a*yid[i]+b*ydj[j]-ydd)/((a-1)*(b-1))
    
    y[i,j] = value
    
    return y


y_est,N_error = estimate_all_missing_values(y_missing)

y_fac = np.transpose(y_est.T.reshape((3,3,2)),axes=(1,0,2))

res = anova(y,N_error)
print_anova(res)




