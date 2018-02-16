# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 13:43:33 2017

@author: anygren
"""

import numpy as np
import cantera as ct
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
import flameletLoader as fl

# Load Data from CFD
CFD_Data = np.genfromtxt('CFDdataNew0.csv', delimiter = ',', skip_header = 1)
Z_CFD    = CFD_Data[:,1]
k        = CFD_Data[:,0]
varZ     = CFD_Data[:,2]
epsilon  = CFD_Data[:,3]
YCoords  = CFD_Data[:,7]
XCoords  = CFD_Data[:,8]

# Setup A Grid for plotting CFD Data
YGrid = np.linspace(0,0.15,50)
XGrid = np.linspace(0,0.05,50)

# Interpolate CFD data onto Python grid:
DNozzle = 7.2e-3
plt.figure(figsize=(4, 9))
plt.title('Mean Mixture Fraction')
varZGrid = griddata(XCoords, YCoords, varZ, XGrid, YGrid, interp='linear')

# Make A  Contour Plot of Mean Mixture Fraction
ZPlot = plt.contourf(XGrid/DNozzle,YGrid/DNozzle, varZGrid, levels = np.linspace(0,1,100))
plt.xlabel('r/D',fontsize = 15)
plt.ylabel('x/D', fontsize = 15)
plt.colorbar(ZPlot)

# Import flamelet library.
flamelets = fl.loadFlamelet('flameletLib')

# Example Usage: Get the temperature and  Scalar Dissipation Rate of the first flamelet:
Z = flamelets.getFlameletScalar(1, 'Z')
T = flamelets.getFlameletScalar(1, 'TEMPERATURE')
yH2O = flamelets.getFlameletScalar(1, 'massfraction-h2o')
Chi1 = flamelets.getFlameletChi(1)

plt.figure()
plt.plot(Z, T)
plt.grid()
plt.xlabel('Z')
plt.ylabel('T [K]')
plt.title('Temperature in mixture fraction space at Chi = ' + str(Chi1))
plt.show()
