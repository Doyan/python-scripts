# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 14:20:27 2020

@author: Gabriel
"""



import numpy as np

import scipy.stats as stats
import matplotlib.pyplot as plt

# =============================================================================

def anova(y,n_missing=0):
    a,b = y.shape
    N = a*b
    
    df = np.array([a-1,b-1,(a-1)*(b-1)-n_missing,N-1],dtype=float)
    df_a = a - 1
    df_b = b - 1
    df_error = df_a*df_b - n_missing
    df_total = N - 1
        
    ydd = np.sum(np.sum(y))
    yid = np.sum(y,axis=1)
    ydj = np.sum(y,axis=0)
    
    SS = np.zeros_like(df)
    
    SS[-1] = np.sum(np.sum(y**2)) - ydd**2/N
    SS[0] = 1/b*np.sum(yid**2) - ydd**2/N
    SS[1] = 1/a*np.sum(ydj**2) - ydd**2/N
    SS[2] = SS[-1] - SS[0] - SS[1]
    
    MS=SS/df
    
    MS[-1] = np.nan
    
    F0 = MS/MS[-2]
    
    F0[-2] = np.nan
    
  
    P0 = np.empty_like(F0)
    for i,f in enumerate(F0):
        if not np.isnan(f):
            P0[i] = stats.f.sf(f,df[i],df[-2])
        else:
            P0[i] = np.nan
            
    index = 'A Block Error Total'.split(' ')
            
    return SS,df,MS,F0,P0,index
    

def contrast_evaluation(y,ci_mat):
    SS,df,MS,F0,P0,index = anova(y)
    
    a,b = y.shape
    
    yidbar = np.average(y,axis=1)
    MS_error = MS[-2]
    df_error = df[-2]
    
    F_values = []
    P_values = []
    for ci in ci_mat: 
        num = np.sum(ci*yidbar)**2
        denum = MS_error/b*np.sum(ci**2)
        F0 = num/denum
        F_values.append(F0)
        P_values.append(stats.f.sf(F0, 1, df_error))
    return F_values, P_values    
    

def contrasts_Scheffe(y,alpha,ci_mat):
    SS,df,MS,F0,P0,index = anova(y)
    
    a,b = y.shape
    
    yidbar = np.average(y,axis=1)
    MS_error = MS[-2]
    df_error = df[-2]
    df_a = df[0]
    
    
    results = []
    conf_intervals = []
    for ci in ci_mat:
        C = np.sum(ci*yidbar)
        S_C = np.sqrt(MS_error*np.sum(ci**2/b))
        
        F_alpha = stats.f.isf(alpha, df_a, df_error)
        
        S_alpha = S_C*np.sqrt(df_a*F_alpha)
        
        results.append(np.abs(C) > S_alpha)
        conf_intervals.append((C - S_alpha, C + S_alpha))
    
    
    return results,conf_intervals



def estimate_missing_value(y,i,j):
    y[i,j] = 0
    a,b = y.shape
    
    ydd = np.sum(np.nansum(y))
    yid = np.nansum(y,axis=1)
    ydj = np.nansum(y,axis=0)
    
    value = (a*yid[i]+b*ydj[j]-ydd)/((a-1)*(b-1))
    
    y[i,j] = value
    
    return y

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

# =============================================================================
    


def plot_residuals_predicted(y):
    
    yidbar = np.average(y,axis=1)
    ydjbar = np.average(y,axis=0)
    yddbar = np.average(np.average(y))
    
    yimat = np.broadcast_to(yidbar.reshape(y.shape[0],1),y.shape)
    yjmat = np.broadcast_to(ydjbar,y.shape)
    
    y_pred = yimat + yjmat - yddbar
    
    residuals = y - y_pred 
    
    
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

    yidbar = np.average(y,axis=1)
    ydjbar = np.average(y,axis=0)
    yddbar = np.average(np.average(y))
    
    yimat = np.broadcast_to(yidbar.reshape(y.shape[0],1),y.shape)
    yjmat = np.broadcast_to(ydjbar,y.shape)
    
    y_pred = yimat + yjmat - yddbar
    
    residuals = y - y_pred 
    
    
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)

    
    stats.probplot(residuals.flatten(),plot=ax)
    ax.set_ylabel('Ordered residuals')
    ax.grid()
    
    return fig,ax

def plot_residuals_ab(y,axis=0):

    yidbar = np.average(y,axis=1)
    ydjbar = np.average(y,axis=0)
    yddbar = np.average(np.average(y))
    
    yimat = np.broadcast_to(yidbar.reshape(y.shape[0],1),y.shape)
    yjmat = np.broadcast_to(ydjbar,y.shape)
    
    y_pred = yimat + yjmat - yddbar
    
    residuals = y - y_pred 
      
    
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)
    
    for i in range(y.shape[axis]):
        if axis == 1:
            ydata = residuals[:,i]
            xstring = 'Blocks'
        else:
            ydata = residuals[i,:]
            xstring = 'Treatments'
        
        ax.plot(i*np.ones_like(ydata)+1,ydata,'b.')
        
    
    xmin = min(ax.lines[0].get_xdata()) -1
    xmax = max(ax.lines[-1].get_xdata()) +1
    
    ax.plot([xmin,xmax],[0,0],'k')
    ax.set_xlim(xmin,xmax)
    ax.set_ylabel('Residuals')
    ax.set_xlabel(xstring)
    ax.set_title('Residuals vs Predicted values')
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




# y = np.array([[90.3, 89.2, 98.2, 93.9, 87.4, 97.9],
#               [92.5, 89.5, 90.6, 94.7, 87.0, 95.8],
#               [85.5, 90.8, 89.6, 86.2, 88.0, 93.4],
#               [82.5, 89.5, 85.6, 87.4, 78.9, 90.7]])


# y = np.array([[685,792,838,875],
#               [722,806,893,953],
#               [733,802,880,941],
#               [811,888,952,1005],
#               [828,920,978,1023]],dtype=float)


# y = np.array([[3,5,10],
#               [4,6,11],
#               [np.nan,6,2],
#               [7,7,5],
#               [0,-4,np.nan],
#               [1,2,-6]],dtype=float) 




y = np.array([[6,5,7,3,5,8,5,5,4,5,4,6,5,3,2,7,6,8,6,7],
             [7,8,5,9,5,8,6,4,6,7,6,9,5,7,4,6,8,5,8,7]])


res = anova(y)
print_anova(res)

