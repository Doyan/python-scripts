# -*- coding: utf-8 -*-
"""
Created on Thu Oct 05 15:58:24 2017

@author: Gabriel
"""

import numpy as np
import cantera as ct
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
import flameletLoader as fl
from scipy.stats import beta

gas = ct.Solution('gri30.xml')

iN2  = gas.species_index('N2')
iO2  = gas.species_index('O2')
iH2O = gas.species_index('H2O')
iF = gas.species_index('CH4')
iCO2= gas.species_index('CO2')
# -----------------------------------------------------------------------------
# Load Data from CFD
CFD_Data = np.genfromtxt('CFDdataNew0.csv', delimiter = ',', skip_header = 1)
Z_CFD    = CFD_Data[:,1]
k        = CFD_Data[:,0]
varZ     = CFD_Data[:,2]
epsilon  = CFD_Data[:,3]

YCoords  = CFD_Data[:,7]
XCoords  = CFD_Data[:,8]

# Calculate stoichiometric chi based on CFD data
cx=2
chist=cx*epsilon*varZ/k

# for later plotting functions
xres=50
yres=50   

# Groundwork for plotting
    
  
DNozzle = 7.2e-3

    
YGrid = np.linspace(0,0.15,yres)
XGrid = np.linspace(0,0.05,xres)

# Interpolate CFD data onto Python grid:


ZGrid = griddata(XCoords, YCoords, Z_CFD, XGrid, YGrid, interp='linear')
varZGrid = griddata(XCoords, YCoords, varZ, XGrid, YGrid, interp='linear')
kGrid = griddata(XCoords, YCoords, k, XGrid, YGrid, interp='linear')
epsGrid = griddata(XCoords, YCoords, epsilon, XGrid, YGrid, interp='linear')
chistGrid = griddata(XCoords, YCoords, chist, XGrid, YGrid, interp='linear')


# -----------------------------------------------------------------------------
# Augment flamelet lib with density

flamelets = fl.loadFlamelet('flameletLib')

for i in range(1,flamelets.nFlamelets):
    Z=flamelets.getFlameletScalar(i,'Z')
    T=flamelets.getFlameletScalar(i,'TEMPERATURE')
    
    nSpecies=len(gas.species_names)
    nGridpoints=len(Z)
    P=ct.one_atm
    
    Y = np.zeros((nSpecies, nGridpoints))
    for s in range(nSpecies):
        Y[s, :]=flamelets.getFlameletScalar(i, 'massfraction-' + str(gas.species_names[s].lower()))
    
    rho=[]
    for j in range(nGridpoints):
        gas.TPY =T[j], P, Y[:, j]
        rho.append(gas.density_mass)
    
    flamelets.setFlameletScalar(i,'rho',rho)    
# -----------------------------------------------------------------------------
# Calculate stoichiometric mixture fraction

MW=gas.molecular_weights[[iF, iO2, iCO2, iH2O, iN2]]
y_f1=0.16
y_fst=(1*MW[0]/(1*MW[0]+2*MW[1]+2*3.76*MW[4]))

Z_st=y_fst/y_f1

# -----------------------------------------------------------------------------
# Functions

# function to make favre averages of scalars from flamelet lib
def favrage(phi,zbar,varz,chist):
    z=flamelets.interpolateFlamelet(chist,'Z', nInterpPoints = 50)
    rho=flamelets.interpolateFlamelet(chist,'rho', nInterpPoints = 50)
    phi=flamelets.interpolateFlamelet(chist,phi, nInterpPoints = 50)
    
    gamma=zbar*(1-zbar)/varz
    a=zbar*gamma
    b=gamma*(1-zbar)

    fz=beta.pdf(z,a,b)
    
    if np.isinf(fz[0]):
        fz[0]=fz[1]/a*1.5
    
    if np.isinf(fz[-1]):
        fz[-1]=fz[-2]/b*1.5
    
    norm=np.trapz(fz)
    fz=fz/norm
    
    rhobar=np.trapz(fz*rho)

    phi_favre=1/rhobar*np.trapz(fz*phi*rho)
    return phi_favre 

# Function to conveniently plot against x/D, r/D axes     
def plotrDxD(grid,title,zmin=-1,zmax=-1,cmap='viridis'):
    if zmin == -1:
        zmin=np.min(grid)
        
    if zmax == -1:
        zmax=np.max(grid)
    
    
    fig = plt.figure(figsize=(5, 9))
    plt.title(title, fontsize = 13)

    cax = plt.contourf(XGrid/DNozzle,YGrid/DNozzle, grid, levels = np.linspace(zmin,zmax,500),cmap=cmap)
    plt.xlabel('Radial distance - r/D',fontsize = 13)
    plt.ylabel('Axial distance - x/D', fontsize = 13)
    cb=plt.colorbar(cax)  
    
    for c in cax.collections:
        c.set_edgecolor("face")

    
    return fig, cax, cb
   
# Function to conveniently generate scalar grid
def genGrid(phi):
    grid=np.zeros_like(ZGrid)
    for i in range(xres):
        for j in range(yres):
            grid[i,j]=favrage(phi,ZGrid[i,j],varZGrid[i,j],chistGrid[i,j]) 
    return grid

# -----------------------------------------------------------------------------
# Calculate mean Temperature field for scatter

Tmean=[]
for i in range(len(Z_CFD)):
    Tmean.append(favrage('TEMPERATURE',Z_CFD[i],varZ[i],chist[i]))

# -----------------------------------------------------------------------------
    
h2oGrid=genGrid('massfraction-ch4')   
TGrid=genGrid('TEMPERATURE')
# -----------------------------------------------------------------------------
# Plotting
   
f,cax,cb= plotrDxD(ZGrid,'Mean mixture fraction',0,1)
cb.ax.plot([0,1],[Z_st,Z_st],'k',label='$Z_{st}$')

ticks=np.insert(np.arange(0,1.01,0.1),4,Z_st)
tick_str=ticks.astype('|S5')
cb.set_ticks(ticks)
cb.set_ticklabels(tick_str)
plt.contour(XGrid/DNozzle,YGrid/DNozzle,ZGrid,[Z_st])
plt.savefig('images/zbarxy.pdf')


plotrDxD(chistGrid,'$\chi_{st}$')
plt.savefig('images/chistxy.pdf')


plotrDxD(TGrid,'Temperature field',np.min(TGrid),np.max(TGrid),'plasma')
plt.contour(XGrid/DNozzle,YGrid/DNozzle,ZGrid,[Z_st])
plt.savefig('images/tbarxy.pdf')


plt.figure(figsize=(9,10))
cax=plt.scatter(Z_CFD,Tmean,s=0.8,c=chist,cmap='viridis')
cb=plt.colorbar(cax)
plt.grid()
cb.set_label('$\chi_{st}$')
plt.xlabel('Mixture fraction - Z')
plt.ylabel('Average Temperature [K]')
plt.title('Temperature in mixture fraction space')
plt.savefig('images/tscatt.pdf')


plotrDxD(h2oGrid,'H2O-massfraction field')


# Plotting flamelet example
plt.figure()
for fno in [1, 2, 20, 50, 122]:
    Z=flamelets.getFlameletScalar(fno,'Z')
    T=flamelets.getFlameletScalar(fno,'TEMPERATURE')
    chi=flamelets.getFlameletChi(fno)
    plt.plot(Z,T,label='$\chi =$' + str(chi))

plt.plot([Z_st, Z_st],[0,2300],'--k',label='$Z_{st}$')
plt.ylim(250,2250)
plt.legend()
plt.title('A few flamelets')
plt.xlabel('Mixture fraction')
plt.ylabel('Temperature [K]')  
plt.savefig('images/flamex.pdf')
  