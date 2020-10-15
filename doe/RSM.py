# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 22:16:11 2020

@author: Gabriel
"""

import numpy as np

import scipy.stats as stats
import matplotlib.pyplot as plt


# =============================================================================



def fit_firstorder(y_exp,x0,x_radius,variance=variance,df_variance=df_variance):
    
    
    x_natural = lambda x_coded: x_coded*x_radius + x0
    
    
    ci_mat = np.array([[-1,1,-1,1,0],[-1,-1,1,1,0]])
    
    x_exp = x_natural(ci_mat.T)
    
    
    y_resp = np.array([y_exp(x) for x in x_exp])    
        
    ci_ab = np.array([1,-1,-1,1,0])
    ci_intercept = np.array([1,1,1,1,1])
    
    
    # Append the contrasts for computing intercept, also the interaction effect 
    ci_mat = np.vstack((ci_intercept,ci_mat,ci_ab))
    
    
    xTx = ci_mat@ci_mat.T
    xTy = ci_mat@y_resp
    beta = np.linalg.inv(xTx)@xTy
    
    t_lsq = beta/np.sqrt(variance/xTx.sum(axis=1))
    P_lsq = stats.t.sf(t_lsq,df_variance)
    
    print('')
    for name,t,P in zip(['B0','B1','B2','B12'],t_lsq,P_lsq):
        print(name,round(t,4),round(P,6))
    print('')
    return beta, t_lsq, P_lsq


def step_steepest(y_exp,x0,B,x_radius,key=1):
    
    if key ==1:
        x_Delta_coded = np.array([1,B[2]/B[1]])
    else:
        x_Delta_coded = np.array([B[1]/B[2],1])
    
    
    x_Delta_natural = x_Delta_coded*x_radius
    
    def iterate(results,xpoints,dx=x_Delta_natural):
        xpoints.append(xpoints[-1] + dx)
        results.append(y_exp(xpoints[-1]))
    
    
    results =  [y_exp(x0)]
    xpoints = [x0]
    while len(results) < 3 or not results[-3] > results[-1]:
        iterate(results,xpoints)
        print(results[-1],xpoints[-1])
    
    imax = np.argmax(np.array(results))
    
    ymax = results[imax]
    xmax = xpoints[imax]
    print('')
    print(ymax,xmax)
    return ymax, xmax, results, xpoints



# y_exp = lambda x: 12 + 15*x[0] + 6*x[1] + 2*x[0]*x[1] - 2*x[0]**2 - 3*x[1]**2

# variance = 0.2
# df_variance = 10
# radius = 0.5
# x_radius = np.array([1,1])*radius
    
# x0 = np.array([2,1.5])

ci_mat = np.array([[1,0.5,-0.5,-1,-0.5,0.5,0],[0,np.sqrt(0.75),np.sqrt(0.75),0,-np.sqrt(0.75),-np.sqrt(0.75),0]],)



xTx = ci_mat@ci_mat.T



H = ci_mat.T@np.linalg.inv(ci_mat@ci_mat.T)@ci_mat
# Contrasts = ci_mat@y_resp
# ci2_sums = np.sum(ci_mat**2,axis=1,dtype=float)




# # Eq 3.26
# t_values = Contrasts/np.sqrt(variance/1*ci2_sums)
# P_values = stats.t.sf(t_values,1)


# # B = Contrasts/ci2_sums


# beta,t,P = fit_firstorder(y_exp, x0, x_radius)

# y_reg = lambda x: beta[0] + beta[1]*x[0] + beta[2]*x[1] 


# ymax,xmax,results,xpoints = step_steepest(y_exp, x0, beta, x_radius)


# beta,t,P = fit_firstorder(y_exp, xmax, x_radius)



# ymax,xmax,results,xpoints = step_steepest(y_exp, xmax, beta, x_radius,key=2)


# beta,t,P = fit_firstorder(y_exp, xmax, x_radius/10)


# ymax,xmax,results,xpoints = step_steepest(y_exp, xmax, beta, x_radius/10,key=1)



