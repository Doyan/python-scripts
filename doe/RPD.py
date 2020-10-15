# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 16:53:24 2020

@author: Gabriel
"""
import numpy as np

import scipy.stats as stats
import matplotlib.pyplot as plt


# =============================================================================

B = np.array([11,-3,2.5,-2,0.6,0.9,0.4],dtype=float)

s2=0.8
s2z=0.3


F = lambda x,z: B[0] + B[1]*x[0] + B[2]*x[1] + B[3]*x[0]*x[1] + B[4]*z + B[5]*z*x[0] + B[6]*z*x[1]
Ez = lambda x: B[0] + B[1]*x[0] + B[2]*x[1] + B[3]*x[0]*x[1]

Vz  = lambda x: s2z*(B[4] + B[5]*x[0] + B[6]*x[1])**2 + s2

x_mat = np.array([
    [-1,-1],
    [1,-1],
    [-1,1],
    [1,1],
    [0,0],
    [0.5,0.5],
    [-0.5,0.5],
    [-0.5,-0.5],
    [0.5,-0.5],
],dtype=float)


alpha = 0.05

t_alpha = stats.t.isf(0.05,1000)


for x in x_mat:
    s  = f'{(x[0],x[1])}\t\t'
    s += f'{round(Ez(x),5)}\t'
    s += f'{round(Vz(x),5)}\t'
    s += f'{round(Ez(x) - np.sqrt(Vz(x))*t_alpha,5)}'
    print(s)


xrange = np.linspace(-1.,1.,50)
yrange = np.linspace(-1.,1.,50)


xi,yi = np.meshgrid(xrange,yrange)


Z = Ez([xi,yi])
ZV = Vz([xi,yi])
Zf = F([xi,yi],-1)


plt.figure()
plt.contour(xi,yi,Z,40)
plt.colorbar()


plt.figure()
plt.contourf(xi,yi,ZV,[0,1,2,3])
plt.colorbar()
plt.plot([-1,1],[-1,-1],'--k')
plt.plot([-1,-1],[-1,1],'--k')
plt.plot([1,1],[-1,1],'--k')
plt.plot([-1,1],[1,1],'--k')


plt.figure()
plt.contourf(xi,yi,Zf,[5,15,20,25])
plt.colorbar()
plt.plot([-1,1],[-1,-1],'--k')
plt.plot([-1,-1],[-1,1],'--k')
plt.plot([1,1],[-1,1],'--k')
plt.plot([-1,1],[1,1],'--k')


plt.figure()
plt.contour(xi,yi,Z,40)
plt.colorbar()











