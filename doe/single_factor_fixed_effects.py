# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 16:44:22 2020

@author: Gabriel
"""

import numpy as np

import scipy.stats as stats
import matplotlib.pyplot as plt
# =============================================================================

def anova(y,n_miss=0):
    
    a = y.shape[0]
    ni = y.shape[1] - np.isnan(y).sum(1)
    
    N = ni.sum()
    
    df_treat = a - 1
    df_error = N - a - n_miss
    df_total = N - 1
    
    
    
    ydd = np.nansum(np.nansum(y))
    yid = np.nansum(y,axis=1)
    
    SS_total = np.nansum(np.nansum(y**2)) - ydd**2/N
    SS_treat = np.nansum(yid**2/ni) - ydd**2/N
    SS_error = SS_total - SS_treat
    
    MS_treat = SS_treat/float(df_treat)
    MS_error = SS_error/float(df_error)
    
    F0 = MS_treat/MS_error
    
    P0 = stats.f.sf(F0, df_treat, df_error)
    
    F_arr  = np.array([F0,np.nan,np.nan])
    MS_arr = np.array([MS_treat,MS_error,np.nan])
    SS_arr = np.array([SS_treat,SS_error,SS_total])
    df_arr = np.array([df_treat,df_error,df_total])
    P_arr = np.array([P0,np.nan,np.nan])
    
    index = 'Treat Error Total'.split(' ')
    return SS_arr,df_arr,MS_arr,F_arr,P_arr,index


def mui_ci(y,alpha,i):
    SS,df,MS,F0,P0,index= anova(y)
    
    ni = y.shape[1] - np.isnan(y).sum(1)
    n_repl = ni[i]
    
    t_val = stats.t.isf(alpha/2,df[1])
    radius = t_val*np.sqrt(MS[1]/n_repl)
    
    y_pred = np.nanmean(y[i,:]) 
    
    return (y_pred - radius,y_pred+radius) 


def contrast_evaluation(y,ci_mat):
    SS,df,MS,F0,P0,index = anova(y)
    
    ni = y.shape[1] - np.isnan(y).sum(1)
    
    
    yidbar = np.average(y,axis=1)
    MS_error = MS[1]
    df_error = df[1]
    
    F_values = []
    P_values = []
    for ci in ci_mat: 
        num = np.sum(ci*yidbar)**2
        denum = MS_error*np.sum(ci**2/ni)
        F0 = num/denum
        F_values.append(F0)
        P_values.append(stats.f.sf(F0, 1, df_error))
    return F_values, P_values    
    

def contrasts_Scheffe(y,alpha,ci_mat):
    SS,df,MS,F0,P0,index = anova(y)
    
    ni = y.shape[1] - np.isnan(y).sum(1)
    
    yidbar = np.average(y,axis=1)
    MS_error = MS[1]
    df_error = df[1]
    df_a = df[0]
    
    
    results = []
    conf_intervals = []
    for ci in ci_mat:
        C = np.sum(ci*yidbar)
        S_C = np.sqrt(MS_error*np.sum(ci**2/ni))
        
        F_alpha = stats.f.isf(alpha, df_a, df_error)
        
        S_alpha = S_C*np.sqrt(df_a*F_alpha)
        
        results.append(np.abs(C) > S_alpha)
        conf_intervals.append((C - S_alpha, C + S_alpha))
    
    
    return results,conf_intervals




def plot_residuals_rtime(y,y_order):
    
    
    yiavg = np.nanmean(y,axis=1)
    y_pred = np.tile(yiavg.reshape((y.shape[0],1)),y.shape[1])
    
    residuals = y -  y_pred 
    
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)

    ax.plot(y_order.flatten(),residuals.flatten(),'.')
    
    
    xmin = min(ax.lines[0].get_xdata()) -1
    xmax = max(ax.lines[0].get_xdata()) +1
    
    ax.plot([xmin,xmax],[0,0],'k')
    ax.set_ylabel('Residuals')
    ax.set_xlabel('Run order')
    ax.set_title('Residuals vs Run order')
    ax.grid()
    ax.set_xlim(xmin,xmax)

    
    return fig,ax




def plot_residuals_predicted(y):
    
    
    yiavg = np.nanmean(y,axis=1)
    y_pred = np.tile(yiavg.reshape((y.shape[0],1)),y.shape[1])
    
    residuals = y -  y_pred 
    
    
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)

    ax.plot(y_pred.flatten(),residuals.flatten(),'.')
    
    xmin = min(ax.lines[0].get_xdata())*0.98
    xmax = max(ax.lines[0].get_xdata())*1.02
    
    ax.plot([xmin,xmax],[0,0],'k')
    ax.set_xlim(xmin,xmax)
    ax.set_ylabel('Residuals')
    ax.set_xlabel('Predicted values')
    ax.set_title('Residuals vs Predicted values')
    ax.grid()
    
    return fig,ax


def plot_residuals_probability(y,fig=None,ax=None):
    if len(y.shape) == 1:
        yiavg = np.nanmean(y)
        
        residuals = y -  yiavg 
    
    else:
        
        yiavg = np.nanmean(y,axis=1)
        y_pred = np.tile(yiavg.reshape((y.shape[0],1)),y.shape[1])
    
        residuals = y -  y_pred 
    
    if fig == None:
        fig = plt.figure(figsize=(8,8))
    
        ax = fig.add_subplot(111)

    
    stats.probplot(residuals.flatten(),plot=ax)
    ax.set_ylabel('Ordered residuals')
    ax.grid()
    
    return fig,ax

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


def t_test_unequal_variance(y,i,j,alpha=0.05):
    
    yiavg = np.average(y,axis=(1))
    y_pred  = np.tile(yiavg.reshape((y.shape[0],1)),y.shape[1])
    
    a,n = y.shape
    
    
    
    S2 = np.sum((y-y_pred)**2,axis=1)/(n-1)
    
    MSi = S2[i]/n
    MSj = S2[j]/n
    
    
    t0 = (yiavg[i] - yiavg[j])/np.sqrt(MSi+MSj)
    df_neq = (MSi + MSj)**2/(MSi**2/(n -1) + MSj**2/(n -1))
    
    P0 = stats.t.sf(np.abs(t0),df_neq)
    print(S2[i],S2[j])
    print(np.abs(t0))
    print(df_neq)
    print(P0*2)
    print(np.abs(t0),stats.t.isf(alpha/2,df_neq))
    
    return




# y = np.array([[575,542,530,539,570],
# [565,593,590,579,610],
# [600,651,610,637,629],
# [725,700,715,685,710]],dtype=float)
    
# y_order = np.array([[13,14,8,5,4],
#            [18,9,6,16,17],
#            [7,19,10,20,1],
#            [2,3,15,11,12]
#            ])
    
# y = np.array([[7.4,6.0,6.9,7.7],
#               [3.1,7.5,5.6,3.0],
#               [5.7,5.9,6.5,5.5]])


#fig,ax = plot_residuals_probability(y[0])
#fig,ax = plot_residuals_probability(y[1],fig,ax)
#fig,ax = plot_residuals_probability(y[2],fig,ax)
#plot_residuals_predicted(y)
#fig,ax = plot_residuals_rtime(y,y_order)




y = np.array([[6,5,7,3,5,8,5,5,4,5,4,6,5,3,2,7,6,8,6,7],
             [7,8,5,9,5,8,6,4,6,7,6,9,5,7,4,6,8,5,8,7]])


res = anova(y)
print_anova(res)


SS,df,MS,F0,P0,index = res

diff = np.average(y[0]-y[1])

t95 = stats.t.isf(0.05/2,38)
t99 = stats.t.isf(0.01/2,38)

estvar = np.sqrt((2*MS[-2])/20) 


Ci_95 = (diff - t95*estvar,diff + t95*estvar)

Ci_99 = (diff - t99*estvar,diff + t99*estvar)
