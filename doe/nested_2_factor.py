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



def anova(y,blocks_random=True):
    
    try:
        a,b,n = y.shape
    except ValueError:
        y = y.reshape((*y.shape,1))
        
        a,b,n = y.shape
        
        
        
        yddd = np.sum(y,axis=(0,1,2))
        yidd = np.sum(y,axis=(1,2))
        yijd = np.sum(y,axis=2)
    
    
    df = np.array([a -1,a*(b-1),a*b*(n-1),a*b*n -1],dtype=float)
        
    SS =  np.zeros_like(df)
    
    SS[0] = 1/(b*n)*np.sum(yidd**2) - yddd**2/(a*b*n)
    SS[1] = 1/n*np.sum(yijd**2,axis=(0,1)) - 1/(b*n)*np.sum(yidd**2)
    SS[2] = np.sum(np.sum(y**2,axis=2),axis=(0,1)) - 1/n*np.sum(yijd**2,axis=(0,1))
    SS[3] = np.sum(np.sum(y**2,axis=2),axis=(0,1)) - yddd**2/(a*b*n)
    
    
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
    
    if blocks_random:
        F0[0] = MS[0]/MS[1]
        P0[0] = stats.f.sf(F0[0],df[0],df[1])
    
    
    index = 'A B(A) Error Total'.split(' ')
    
    
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





# y = np.array([[[26,28,28],[28,28,27]],
#              [[25,26,29],[31,30,29]],
#              [[33,29,31],[30,31,33]]],dtype=float)


def pure_error_per_treatment(y):
    yijdbar = np.nanmean(y,axis=2)
    yijdbar_mat = np.broadcast_to(yijdbar.reshape(*yijdbar.shape,1),y.shape)
    
    Si = np.sum((y - yijdbar_mat)**2,axis=(1,2))
    return Si/(y.shape[0] -1) 

y = np.array([[[59,57],[46,48]],
     [[46,48],[45,46]],
     [[49,48],[34,35]],
     [[53,52],[46,44]],
     [[39,40],[38,37]]])




y = np.array([[6,5,7,3,5,8,5,5,4,5,4,6,5,3,2,7,6,8,6,7],
             [7,8,5,9,5,8,6,4,6,7,6,99,5,7,4,6,8,5,8,7]])

res = anova(y)
print_anova(res)