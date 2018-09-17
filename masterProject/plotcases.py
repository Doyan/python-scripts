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


def func(x,a,b):
    return a*(x)**(-1) + b*x

def curveChecked(f,xdata,ydata):
    popt,cov = curve_fit(f, xdata,ydata)
    residuals = ydata - f(xdata,*popt)
    ss_res = np.sum(residuals**2)
    
    ss_tot = np.sum((ydata-np.mean(ydata))**2)
    r2 = 1 - (ss_res / ss_tot)
    return popt, r2

def genCurveData(f,xdata,ydata,I=[0,110],N=100):
    p,r2 = curveChecked(f,xdata,ydata )
    x = np.linspace(I[0],I[1],N)
    y = f(x,*p)

    return x,y,r2,p


def plotFitted(f,xdata,ydata,label='',I=[0,110],N=100):
    if I[0] == 0:
        I[0] = I[0] + 0.00001*I[-1]
        
    x,y,r2,p = genCurveData(f,xdata,ydata,I,N)
    ax = plt.plot(xdata,ydata,'o')
    ax2 = plt.plot(x,y,color=ax[0].get_c(), label=label)
    return ax2

#%%


fnam = 'reportdata.csv'
fnam2 = 'result'

D1=pd.read_csv(fnam2)
D1.spec = '$x_{H_2O} = 0.6$, $X_{ch} = 0.1$'

D2=pd.read_csv(fnam)
D2.spec = '$x_{H_2O} = 0.4$, $X_{ch} = 0.2$'
    


plt.figure(1)
for D in [D1,D2]:
    dTavg = D.T_C_avg - D.T_G_avg
    dTmax = D.T_C_max - D.T_G_min
    
    #ax= plotFitted(func,D.case,dTavg,'$\Delta T_{avg}$ - ' + D.spec)
    #ax[0].set_linestyle('--')
    #color = ax[0].get_c()
    
    ax= plotFitted(func,D.case,dTmax,'$\Delta T_{max}$ - ' + D.spec)
    #ax[0].set_c(color)
    




plt.xlabel('Conduction scaling factor:   ($k_{eff}/1800$ W/mK )')
plt.ylabel('$\Delta$ T [$\degree C$]' )
plt.yticks([10, 50, 100, 150, 200, 250, 300, 400])
plt.xticks([0, 10, 20, 40, 60, 80, 100])

plt.ylim(0,400)
plt.xlim(1,110)
plt.title('$\Delta$ T vs scaling factor')      
   
plt.legend()
plt.grid()
    
plt.savefig('dtmaxplot.pdf')



#    pavg,r2_avg = curveChecked(func, D.case,dTavg )
#    pmax,r2_max = curveChecked(func, D.case,dTmax )
#
#
#    x = np.linspace(0,110,100)
#
#    plt.plot(D.case,dTavg,'o',label='$\Delta T_{avg}$')
#    plt.plot(D.case,dTmax,'o',label='$\Delta T_{max}$')
#    plt.plot(x,func(x,*pavg),'b--', label='$f = %2.3f \; x^{-1} + %2.3f $' % tuple(pavg))
#    plt.plot(x,func(x,*pmax),'g--', label='$f = %2.3f \; x^{-1} + %2.3f $' % tuple(pmax))

#%%
fnam = 'keff33.csv'  

df=pd.read_csv(fnam)


D1 = df[0:11]
D1.spec = 'xH2O = 0.2'
D2 = df[11:22]
D2.spec = 'xH2O = 0.4'

D3 = df[22:33]
D3.spec = 'xH2O = 0.6'


DF = [D1, D2, D3]


def func(x,a,b,c):
    return a*(x)**(-c) + b

for D in DF:    
    dTmax = pd.to_numeric(D.T_C_max) - pd.to_numeric(D.T_G_min)
    ax= plotFitted(func,pd.to_numeric(D.K_eff),dTmax,'$\Delta T_{max}$ - ' + D.spec,[0.01,320000])


plt.ylim(0,300)
plt.xlim(1,200000)
plt.grid()
plt.title('Maximum $\Delta T$ vs effective conduction')
plt.xlabel('$k_{eff}$    [W/mK]')
plt.ylabel('$\Delta T$    [K]')
plt.legend()
plt.savefig('new-k-vs-dT-max.pdf')
plt.show()


for D in DF:
 
    dTavg = pd.to_numeric(D.T_C_avg) - pd.to_numeric(D.T_G_avg)
    
    ax= plotFitted(func,pd.to_numeric(D.K_eff),dTavg,'$\Delta T_{avg}$ - ' + D.spec,[0.01,320000])


plt.ylim(0,300)
plt.xlim(1,200000)
plt.grid()
plt.title('Average $\Delta T$ vs effective conduction')
plt.xlabel('$k_{eff}$    [W/mK]')
plt.ylabel('$\Delta T$    [K]')
plt.legend()
plt.savefig('new-k-vs-dT-avg.pdf')
plt.show()


