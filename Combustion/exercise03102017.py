# -*- coding: utf-8 -*-
"""
Created on Tue Oct 03 15:16:06 2017

@author: Gabriel
"""
import numpy as np, matplotlib.pyplot as plt, cantera as ct
from scipy.stats import beta

#------------------------------------------------------------------------------
zbar=0.6

zvar=[0.160, 0.080, 0.048, 0.022]

z=np.linspace(0,1)

for var in zvar:
    gamma=zbar*(1-zbar)/var
    a=zbar*gamma
    b=gamma*(1-zbar)
    plt.plot(z,beta.pdf(z,a,b),label='$\sigma^2 = $' + str(var))


plt.xlabel('Mixture fraction')
plt.ylabel('Beta-PDF')
plt.title('Mean mixture fraction = 0.6')
plt.legend()    
plt.show()    

#------------------------------------------------------------------------------
# uppg 2
p=ct.one_atm
fnam='ref.txt'

gas = ct.Solution('gri30.xml', 'gri30_mix')

iN2  = gas.species_index('N2')
iO2  = gas.species_index('O2')
iH2O = gas.species_index('H2O')
iF = gas.species_index('CH4')
iCO2= gas.species_index('CO2')

Ms=gas.molecular_weights[[iF, iO2, iCO2, iH2O, iN2]]
Rs=ct.gas_constant/Ms

data  = np.genfromtxt(fnam) #Input 
z = data[:,0]
t = data[:,1]
yf= data[:,2]
yo2= data[:,3]
yco2= data[:,4]
yh2o= data[:,5]
yn2= data[:,6]

ys=data[:,2:7]

rho=[]
for i in range(0,len(t)):
    gas.TP=[t[i],p]
    gas.Y[[iF, iO2, iCO2, iH2O, iN2]]=ys[i]
    rho.append(gas.density)
# p = rho T sum(Rs*Ys) Rs= Rm/Ms

zbar=0.6
zvar=0.022
gamma=zbar*(1-zbar)/var
a=zbar*gamma
b=gamma*(1-zbar)


fz=beta.pdf(z,a,b)
norm=np.trapz(fz)
fz=fz/norm

ysbar=[]
rhobar= np.trapz(rho*fz)
tbar=1/rhobar*np.trapz(rho*t*fz)
yFbar=1/rhobar*np.trapz(rho*ys[:,0]*fz)



