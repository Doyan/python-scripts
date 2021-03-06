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

#%%


fnam = 'keff40.csv'  

df=pd.read_csv(fnam)


D1 = df[0:9]
D1.spec = '$x_{H_2O}:\; 0.2, \;\; k_{wall}:\;1$'

D2 = df[9:18]
D2.spec = '$x_{H_2O}:\; 0.2, \;\; k_{wall}:\;0.5$'

D3 = df[18:27]
D3.spec = '$x_{H_2O}:\; 0.1, \;\; k_{wall}:\;1$'

D4 = df[27:36]
D4.spec = '$x_{H_2O}:\; 0.1, \;\; k_{wall}:\;0.5$'

Bonus=df[36:40]

DF = [D1, D2, D3, D4]

def func(x,a,b,c):
    return a*(x)**(-c) + b


plt.figure(1,(8,6))
for D in DF:    
    dTmax = pd.to_numeric(D.T_C_max) - pd.to_numeric(D.T_G_min)
    ax= plotFitted(func,pd.to_numeric(D.K_eff),dTmax,D.spec,[0.01,80000])

dTmax_bonus =pd.to_numeric(Bonus.T_C_max) - pd.to_numeric(Bonus.T_G_min)
dTmax_D2 = pd.to_numeric(D2.T_C_max) - pd.to_numeric(D2.T_G_min)
dTmax_D4 = pd.to_numeric(D4.T_C_max) - pd.to_numeric(D4.T_G_min)


plt.plot(Bonus.K_eff.iloc[0],dTmax_bonus.iloc[0],'d',color=[0.1,0.3,0.1],label='$x_{H_2O}:\; 0.2, \;\; d_{wall}:\;0.15$')
plt.plot(Bonus.K_eff.iloc[1],dTmax_bonus.iloc[1],'x',color='r',label='$x_{H_2O}:\; 0.2, \;\; w_{wall}:\;1.0$')
plt.plot(Bonus.K_eff.iloc[2],dTmax_bonus.iloc[2],'d',color=[0.5,0.5,0.9],label='$x_{H_2O}:\; 0.1, \;\; d_{wall}:\;0.15$')
plt.plot(Bonus.K_eff.iloc[3],dTmax_bonus.iloc[3],'x',color='b',label='$x_{H_2O}:\; 0.1, \;\; w_{wall}:\;1.0$')

plt.xlim(1,80000)
plt.grid()
plt.title('Maximum $\Delta T$ vs effective conduction')
plt.xlabel('$k_{eff}$    [W/mK]')
plt.ylabel('$\Delta T$    [K]')

plt.ylim(0,300)
plt.legend()
plt.savefig('new-k-vs-dT-max.pdf')
plt.show()

plt.figure(2,(8,6))
plt.plot(Bonus.K_eff.iloc[0],dTmax_bonus.iloc[0],'d',color=[0.1,0.3,0.1],label='$x_{H_2O}:\; 0.2, \;\; d_{wall}:\;0.15$')
plt.plot(Bonus.K_eff.iloc[2],dTmax_bonus.iloc[2],'d',color=[0.5,0.5,0.9],label='$x_{H_2O}:\; 0.1, \;\; d_{wall}:\;0.15$')

plt.plot(D2.K_eff.iloc[4],dTmax_D2.iloc[4],'o',color=[0.2,0.5,0.2],label='$x_{H_2O}:\; 0.2, \;\; w_{thinwall}:\;0.8$')

plt.plot(Bonus.K_eff.iloc[1],dTmax_bonus.iloc[1],'x',color='r',label='$x_{H_2O}:\; 0.2, \;\; w_{thinwall}:\;1.0$')

plt.plot(D4.K_eff.iloc[4],dTmax_D4.iloc[4],'o',color='cyan',label='$x_{H_2O}:\; 0.1, \;\; w_{thinwall}:\;0.8$')
plt.plot(Bonus.K_eff.iloc[3],dTmax_bonus.iloc[3],'x',color='b',label='$x_{H_2O}:\; 0.1, \;\; w_{thinwall}:\;1.0$')

plt.grid()
plt.title('Maximum $\Delta T$ vs effective conduction')
plt.xlabel('$k_{eff}$    [W/mK]')
plt.ylabel('$\Delta T$    [K]')
plt.legend()
plt.savefig('bonus-k-vs-dT-max.pdf')


plt.figure(3,(8,6))
for D in DF:
 
    dTavg = pd.to_numeric(D.T_C_avg) - pd.to_numeric(D.T_G_avg)
    
    ax= plotFitted(func,pd.to_numeric(D.K_eff),dTavg,D.spec,[0.01,80000])

dTavg_bonus =pd.to_numeric(Bonus.T_C_avg) - pd.to_numeric(Bonus.T_G_avg)
dTavg_D2 = pd.to_numeric(D2.T_C_avg) - pd.to_numeric(D2.T_G_avg)
dTavg_D4 = pd.to_numeric(D4.T_C_avg) - pd.to_numeric(D4.T_G_avg)


plt.plot(Bonus.K_eff.iloc[0],dTavg_bonus.iloc[0],'d',color=[0.1,0.3,0.1],label='$x_{H_2O}:\; 0.2, \;\; d_{wall}:\;0.15$')
plt.plot(Bonus.K_eff.iloc[1],dTavg_bonus.iloc[1],'x',color='r',label='$x_{H_2O}:\; 0.2, \;\; w_{wall}:\;1.0$')
plt.plot(Bonus.K_eff.iloc[2],dTavg_bonus.iloc[2],'d',color=[0.5,0.5,0.8],label='$x_{H_2O}:\; 0.1, \;\; d_{wall}:\;0.15$')
plt.plot(Bonus.K_eff.iloc[3],dTavg_bonus.iloc[3],'x',color='b',label='$x_{H_2O}:\; 0.1, \;\; w_{wall}:\;1.0$')

plt.ylim(0,300)
plt.xlim(1,80000)
plt.grid()
plt.title('Average $\Delta T$ vs effective conduction')
plt.xlabel('$k_{eff}$    [W/mK]')
plt.ylabel('$\Delta T$    [K]')
plt.legend()
plt.savefig('new-k-vs-dT-avg.pdf')
plt.show()


plt.figure(4,(8,6))
plt.plot(Bonus.K_eff.iloc[0],dTavg_bonus.iloc[0],'d',color=[0.1,0.3,0.1],label='$x_{H_2O}:\; 0.2, \;\; d_{wall}:\;0.15$')
plt.plot(Bonus.K_eff.iloc[2],dTavg_bonus.iloc[2],'d',color=[0.5,0.5,0.9],label='$x_{H_2O}:\; 0.1, \;\; d_{wall}:\;0.15$')

plt.plot(D2.K_eff.iloc[4],dTavg_D2.iloc[4],'o',color=[0.2,0.5,0.2],label='$x_{H_2O}:\; 0.2, \;\; w_{thinwall}:\;0.8$')

plt.plot(Bonus.K_eff.iloc[1],dTavg_bonus.iloc[1],'x',color='r',label='$x_{H_2O}:\; 0.2, \;\; w_{thinwall}:\;1.0$')

plt.plot(D4.K_eff.iloc[4],dTavg_D4.iloc[4],'o',color='cyan',label='$x_{H_2O}:\; 0.1, \;\; w_{thinwall}:\;0.8$')
plt.plot(Bonus.K_eff.iloc[3],dTavg_bonus.iloc[3],'x',color='b',label='$x_{H_2O}:\; 0.1, \;\; w_{thinwall}:\;1.0$')

plt.grid()
plt.title('Average $\Delta T$ vs effective conduction')
plt.xlabel('$k_{eff}$    [W/mK]')
plt.ylabel('$\Delta T$    [K]')
plt.legend()
plt.savefig('bonus-k-vs-dT-avg.pdf')




