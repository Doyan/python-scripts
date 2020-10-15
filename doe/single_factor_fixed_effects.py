# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 16:44:22 2020

@author: Gabriel
"""

import numpy as np

import scipy.stats as stats
import matplotlib.pyplot as plt
# =============================================================================

def anova(y):
    
    a = y.shape[0]
    ni = y.shape[1] - np.isnan(y).sum(1)
    
    N = ni.sum()
    
    df_treat = a - 1
    df_error = N - a
    df_total = N - 1
    
    
    
    ydd = np.nansum(np.nansum(y))
    yid = np.nansum(y,axis=1)
    
    SS_total = np.nansum(np.nansum(y**2)) - ydd**2/N
    SS_treat = np.nansum(yid**2/ni) - ydd**2/N
    SS_error = SS_total - SS_treat
    
    MS_treat = SS_treat/df_treat
    MS_error = SS_error/df_error
    
    F0 = MS_treat/MS_error
    
    P0 = stats.f.sf(F0, df_treat, df_error)
    
    F_tup  = (F0)
    MS_tup = (MS_treat,MS_error)
    SS_tup = (SS_treat,SS_error,SS_total)
    df_tup = (df_treat,df_error,df_total)
    
    return P0, F_tup, MS_tup, SS_tup, df_tup 


def mui_ci(y,alpha,i):
    P0,F,MS,SS,df= anova(y)
    
    ni = y.shape[1] - np.isnan(y).sum(1)
    n_repl = ni[i]
    
    t_val = stats.t.isf(alpha/2,df[1])
    radius = t_val*np.sqrt(MS[1]/n_repl)
    
    y_pred = np.nanmean(y[i,:]) 
    
    return (y_pred - radius,y_pred+radius) 


def contrast_evaluation(y,ci_mat):
    P0, F_tup, MS_tup, SS_tup, df_tup = anova(y)
    
    ni = y.shape[1] - np.isnan(y).sum(1)
    
    
    yidbar = np.average(y,axis=1)
    MS_error = MS_tup[1]
    df_error = df_tup[1]
    
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
    P0, F_tup, MS_tup, SS_tup, df_tup = anova(y)
    
    ni = y.shape[1] - np.isnan(y).sum(1)
    
    yidbar = np.average(y,axis=1)
    MS_error = MS_tup[1]
    df_error = df_tup[1]
    df_a = df_tup[0]
    
    
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


def plot_residuals_probability(y):

    yiavg = np.nanmean(y,axis=1)
    y_pred = np.tile(yiavg.reshape((y.shape[0],1)),y.shape[1])
    
    residuals = y -  y_pred 
    
    
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)

    
    stats.probplot(residuals.flatten(),plot=ax)
    ax.set_ylabel('Ordered residuals')
    ax.grid()
    
    return fig,ax


y = np.array([[575,542,530,539,570],
[565,593,590,579,610],
[600,651,610,637,629],
[725,700,715,685,710]],dtype=float)
    
y_order = np.array([[13,14,8,5,4],
           [18,9,6,16,17],
           [7,19,10,20,1],
           [2,3,15,11,12]
           ])
    
plot_residuals_probability(y)
plot_residuals_predicted(y)
fig,ax = plot_residuals_rtime(y,y_order)