#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 23 13:40:29 2018

@author: gabgus
"""

import numpy as np, matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit

# -------------------------------------------------------------------------

fnam = 'reportdata.csv'
fnam2 = 'result'

D=pd.read_csv(fnam)

dTavg = D.T_C_avg - D.T_G_avg
dTmax = D.T_C_max - D.T_G_min


def func(x,a,b):
    return a*(x)**(-1) + b

def curveChecked(f,xdata,ydata):
    popt,cov = curve_fit(f, xdata,ydata)
    residuals = ydata - f(xdata,*popt)
    ss_res = np.sum(residuals**2)
    
    ss_tot = np.sum((ydata-np.mean(ydata))**2)
    r2 = 1 - (ss_res / ss_tot)
    return popt, r2

pavg,r2_avg = curveChecked(func, D.case,dTavg )
pmax,r2_max = curveChecked(func, D.case,dTmax )


x = np.linspace(0,110,100)

plt.plot(D.case,dTavg,'o',label='$\Delta T_{avg}$')
plt.plot(D.case,dTmax,'o',label='$\Delta T_{max}$')
plt.plot(x,func(x,*pavg),'b--', label='$f = %2.3f \; x^{-1} + %2.3f $' % tuple(pavg))
plt.plot(x,func(x,*pmax),'g--', label='$f = %2.3f \; x^{-1} + %2.3f $' % tuple(pmax))


plt.xlabel('Conduction scaling factor')
plt.ylabel('$\Delta$ T [$\degree C$]' )

plt.ylim(0,400)
plt.xlim(1,110)
plt.title('$\Delta$ T vs scaling factor -- x$_{H2O}$ = 0.4, X$_{char}$ = 0.2')
plt.legend()
plt.savefig('dtplot.pdf')