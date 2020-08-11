# -*- coding: utf-8 -*-
"""
Created on Sun Aug  9 18:07:19 2020

@author: Gabriel
"""


import numpy as np

import scipy.stats as stats
import matplotlib.pyplot as plt

# =============================================================================


# y = np.array([
#     [[1,-1,0],[-2,-3,-4],[-2,0,1],[1,4,0]],
#     [[1,-2,-3],[0,4,2],[-1,0,-2],[0,3,2]],
#     [[2,4,0],[-2,0,2],[1,-1,2],[3,2,1]]
#     ],dtype=float)



def anova(y,random=(0,1)):
    a,b,n = y.shape
    
    yddd = np.sum(y,axis=(0,1,2))
    yidd = np.sum(y,axis=(1,2))
    yijd = np.sum(y,axis=2)
    
    df_a = a - 1
    df_bwa = a*(b - 1)
    df_error = a*b*(n - 1)
    df_total = a*b*n - 1
    
    SS_a = 1/(b*n)*np.sum(yidd**2) - yddd**2/(a*b*n)
    SS_bwa = 1/n*np.sum(yijd**2,axis=(0,1)) - 1/(b*n)*np.sum(yidd**2)
    SS_error = np.sum(np.sum(y**2,axis=2),axis=(0,1)) - 1/n*np.sum(yijd**2,axis=(0,1))
    SS_total = np.sum(np.sum(y**2,axis=2),axis=(0,1)) - yddd**2/(a*b*n)
    
    MS_a = SS_a/df_a
    MS_bwa = SS_bwa/df_bwa
    MS_error = SS_error/df_error
    
    
    F = stats.f.sf
    
    
    Fdict = {(0,0): [MS_a/MS_error,MS_bwa/MS_error],
             (0,1): [MS_a/MS_bwa,MS_bwa/MS_error],
             (1,1): [MS_a/MS_bwa,MS_bwa/MS_error]}
    
    Pdict = {(0,0): [F(Fdict[(0,0)][0],df_a,df_error),
                      F(Fdict[(0,0)][1],df_bwa,df_error)],
              (0,1): [F(Fdict[(0,1)][0],df_a,df_bwa),
                      F(Fdict[(0,0)][1],df_bwa,df_error)],
              (1,1): [F(Fdict[(0,1)][0],df_a,df_bwa),
                      F(Fdict[(0,0)][1],df_bwa,df_error)]}
    
    
    return Pdict[random],Fdict[random],(MS_a,MS_bwa,MS_error),(SS_a,SS_bwa,SS_error,SS_total),(df_a,df_bwa,df_error,df_total)





def plot_residuals_predicted(y):
    
    yijbar = np.nanmean(y,axis=2)
    y_pred = np.broadcast_to(yijbar.reshape(*yijbar.shape,1),y.shape)
    
    residuals = y - y_pred 
    
    
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)

    ax.plot(y_pred.flatten(),residuals.flatten(),'.')
    
    a = min(ax.lines[0].get_xdata()) 
    b = max(ax.lines[0].get_xdata())
    
    xmin = a - (b-a)*0.05
    xmax = b + (b-a)*0.05
    
    ax.plot([xmin,xmax],[0,0],'k')
    ax.set_xlim(xmin,xmax)
    ax.set_ylabel('Residuals')
    ax.set_xlabel('Predicted values')
    ax.set_title('Residuals vs Predicted values')
    ax.grid()
    
    return fig,ax


def plot_residuals_probability(y):

    yijbar = np.nanmean(y,axis=2)
    y_pred = np.broadcast_to(yijbar.reshape(*yijbar.shape,1),y.shape)
    
    residuals = y - y_pred 
        
    
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)

    
    stats.probplot(residuals.flatten(),plot=ax)
    ax.set_ylabel('Ordered residuals')
    ax.grid()
    
    return fig,ax

def plot_residuals_ab(y,axis=0):

    yijbar = np.nanmean(y,axis=2)
    y_pred = np.broadcast_to(yijbar.reshape(*yijbar.shape,1),y.shape)
    
    residuals = y - y_pred 
      
    
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)
    
    for i in range(y.shape[axis]):
        if axis == 1:
            ydata = residuals[:,i,:]
            xstring = 'Blocks'
        else:
            ydata = residuals[i,:,:]
            xstring = 'Treatments'
        
        ax.plot((i*np.ones_like(ydata)+1).flatten(),ydata.flatten(),'b.')
        
    
    xmin = min(ax.lines[0].get_xdata()) -1
    xmax = max(ax.lines[-1].get_xdata()) +1
    
    ax.plot([xmin,xmax],[0,0],'k')
    ax.set_xlim(xmin,xmax)
    ax.set_ylabel('Residuals')
    ax.set_xlabel(xstring)
    ax.set_title('Residuals vs Treatments/blocks')
    ax.grid()
    
    return fig,ax





y = np.array([[[26,28,28],[28,28,27]],
             [[25,26,29],[31,30,29]],
             [[33,29,31],[30,31,33]]],dtype=float)


def pure_error_per_treatment(y):
    yijdbar = np.nanmean(y,axis=2)
    yijdbar_mat = np.broadcast_to(yijdbar.reshape(*yijdbar.shape,1),y.shape)
    
    Si = np.sum((y - yijdbar_mat)**2,axis=(1,2))
    return Si/(y.shape[0] -1) 



