# -*- coding: utf-8 -*-
"""
Created on Fri Sep 08 15:13:23 2017

@author: Gabriel
"""

import cantera as ct, numpy as np, matplotlib.pyplot as plt
from scipy.integrate import ode

gas=ct.Solution('oneStep_ethane.cti')
ns  = gas.n_species
iN2  = gas.species_index('N2')
iO2  = gas.species_index('O2')
iH2O = gas.species_index('H2O')
iF = gas.species_index('C2H6')
iCO2= gas.species_index('CO2')
# 2C2H6 + 7(O2 + 3.76N2) => 4CO2 + 6H2O + 26.32N2
#------------------------------------------------------------------------------
# initial state of the mixture 
T0 = 300
p0 = ct.one_atm
phi=0.7
CR=10

gas.TP = [T0, p0]

gas.set_equivalence_ratio(phi, 'C2H6', 'O2:1.0, N2:3.76')

gamma = gas.cp_mole/gas.cv_mole

p = p0 * pow(CR,gamma)
T1 = T0 * pow(CR,(gamma -1))
gas.TP=[T1, p]

Y0 = gas.Y
rho0 = gas.density
#------------------------------------------------------------------------------
def dy_rho(t, y, rho, gas):
    dydt = np.zeros(y.size)
    T = y[0]
    Y = y[1:]
    
    gas.TDY=[T, rho, Y]
     
    M=gas.molecular_weights
    Cv=gas.cv_mass
    u_s=gas.standard_int_energies_RT * ct.gas_constant * T / M  # kJ/kg
    wdot_s=gas.net_production_rates   # kmol/m3/s
    
    
    dYdt= M/rho * wdot_s 
    dTdt= - (u_s*dYdt).sum() *1/Cv
    
    dydt[0]=dTdt
    dydt[1:]=dYdt
    return dydt

#------------------------------------------------------------------------------
y0 = np.zeros(ns+1)
y0[0] = T1
y0[1:] = Y0

t0=0.
t1=0.005
nt=400

solver = ode(dy_rho).set_integrator('vode', method='bdf', order=15)
solver.set_initial_value(y0, t0).set_f_params(rho0, gas)

gas()
    
tsol = []
Tsol = []

N2sol = []
Fsol =  []
O2sol = []
H2Osol = []
CO2sol = []

while solver.successful() and np.any(solver.t < t1):
    solver.integrate(t1,step=True)
    tsol.append(solver.t)
    Tsol.append(solver.y[0])
    
    N2sol.append(solver.y[iN2+1])
    Fsol.append(solver.y[iF+1])
    O2sol.append(solver.y[iO2+1])
    H2Osol.append(solver.y[iH2O+1])
    CO2sol.append(solver.y[iCO2+1])
    
tsol=np.array(tsol)
plt.plot(tsol*1000,Tsol,label='Temperature')
plt.xlabel('Time [ms]')
plt.ylabel('Temperature [K]')
plt.title('Temperature, over time fixed mass, constant volume')
plt.legend()
plt.savefig('TempA.pdf')
plt.show()

plt.plot(tsol*1000,Fsol,label='fuel')
plt.plot(tsol*1000,O2sol,label='O2')
plt.plot(tsol*1000,CO2sol,label='CO2')
plt.plot(tsol*1000,H2Osol,label='H2O')
plt.xlabel('Time [ms]')
plt.ylabel('Species massfraction')
plt.title('Massfractions over time, fixed mass constant volume')
plt.legend()
plt.savefig('YA.pdf')
plt.show()


gas()

