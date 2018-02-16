# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 09:15:11 2017

Analytical investigation of flame temperature in a jet 
"""

import numpy as np, matplotlib as mpl, cantera as ct
#from scipy.optimize import fsolve
plt = mpl.pyplot
colors = mpl.colors

gas = ct.Solution('gri30.xml')
iN2  = gas.species_index('N2')
iO2  = gas.species_index('O2')
iH2O = gas.species_index('H2O')
iF = gas.species_index('CH4')
iCO2= gas.species_index('CO2')

# -----------------------------------------------------------------------------
# Data

v_e = 0.1       # m/s
T_e = 473.0     # K
D = 10.0e-3     # m
mu = 17.0e-6    # Pa s
T_air = 300.0   # K
P = ct.one_atm  # Pa

R= D / 2        # Exit radius

gas.TPX = [T_e, P,"CH4:1"]
rho_e = gas.density

gas.TPX=[T_air,P,"O2:1.0, N2:3.76"]
Y_O2_2=gas.Y[iO2]
gas.set_equivalence_ratio(1,'CH4','O2:1.0, N2:3.76')

Y_Fst = gas.Y[iF]
Y_O2st=gas.Y[iO2]

H_st=gas.enthalpy_mass

gas.HPX=[H_st,P,'CO2:1, H2O:2, N2:6.72']
T_ad=gas.T

Y_Pst=gas.Y[[iCO2, iH2O]]

#M=gas.molecular_weights[[iF, iO2, iH2O, iCO2]]


# -----------------------------------------------------------------------------
# Equations

Re_j = (rho_e * v_e * R) / mu
J_e = rho_e * v_e**2 * np.pi * R**2

# Correlation for fuel massfraction 
def Y_F(x,r):
    Xi = pow(3*rho_e * J_e /(16 * np.pi), 0.5) * 1 / mu * r / x
    Y_F = 0.375 * Re_j * pow(x/R,-1) * pow(1 + Xi**2 /4, -2)
    return Y_F


# -----------------------------------------------------------------------------
# Calculation of T(x,r) assumes linear Burke-schulmann with Tmax=T_ad
Z_st = Y_Fst

a1=(T_ad - T_air) / Z_st       # For Fuel lean mixture 
b1=T_air                        #

a2=(T_e - T_ad)/(1 - Z_st)     # For Fuel rich mixture
b2=T_e - a2                     #

a3=(T_e - T_air)                # For unburnt mixture
b3=T_air                        # 

def T_Z(Z):
    condlist= [Z < Z_st, Z >= Z_st]
    choicelist= [a1*Z + b1, a2*Z + b2]    
    return np.select(condlist,choicelist)

def T_xr(x,r):
    Z=Y_F(x,r)
    T=T_Z(Z)
    return T

def T_u(Z):
    return a3*Z + b3

# -----------------------------------------------------------------------------
# Additional caclulations for easy plotting

# Correlation for Y_F only valid a certain distance from nozzle 
xmin = 0.375 * Re_j * R

# Approximate flame length for convenience when plotting
#def f(x):
#    return Y_F(x,0) - Z_st
#
#Lf = fsolve(f, 0.4)
#Lf=Lf[0]

Lf=0.4130 # m -result of above calculation


# -----------------------------------------------------------------------------
# Making vectors

delta = 0.0001
x = np.arange(xmin*1.3, Lf*1.3, delta)
r = np.arange(-25*R,25*R , delta)
xgrid, rgrid = np.meshgrid(x, r)

zgrid = Y_F(xgrid,rgrid)
Tgrid = T_Z(zgrid)

pzgrid=zgrid/Z_st


V=np.array([0, 0.1, 0.75, 0.95, 1, 1.05, 1.75, 2.0, 4, 6, 8, 10])


# -----------------------------------------------------------------------------
# Plotting

cax=plt.contourf(xgrid,rgrid,pzgrid, V, cmap='viridis',)
plt.contour(xgrid,rgrid,pzgrid,[1],linewidth=0.7)
plt.colorbar(cax)
plt.xlabel('Distance from nozzle, [m]')
plt.ylabel('Radial distance, [m]')
plt.title('Mixture fraction - $Z/Z_{st}$')
plt.savefig('images/zcontour.pdf')
plt.show()

cax=plt.contourf(xgrid,rgrid,Tgrid,cmap='plasma')
plt.contour(xgrid,rgrid,Tgrid,[2471])
#plt.contour(xgrid,rgrid,zgrid,[Z_st],colors='0.9')
plt.colorbar(cax)
plt.xlabel('Distance from nozzle, [m]')
plt.ylabel('Radial distance, [m]')
plt.title('Temperature profile - Burke-Schumann')

plt.legend()
plt.savefig('images/tcontour.pdf')
plt.show()



#print R**2*np.pi*v_e*rho_e

