# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 14:25:59 2017

@author: Gabriel
"""
import cantera as ct, numpy as np, matplotlib.pyplot as plt
from scipy.integrate import ode

gas=ct.Solution('oneStep_nheptan_atom.cti')
ns  = gas.n_species
iN2  = gas.species_index('N2')
iO2  = gas.species_index('O2')
iH2O = gas.species_index('H2O')
iF = gas.species_index('C7H16')
iCO2= gas.species_index('CO2')
# 2C2H6 + 7(O2 + 3.76N2) => 4CO2 + 6H2O + 26.32N2
#------------------------------------------------------------------------------
# Volume calculations
Vmax=3.68e-4    # m3
rc=15          # compression ratio
d=7.0e-2        # m
rpm=500*2         # 360deg/min

Vmin=Vmax/rc    # m3
omega=np.deg2rad(rpm*360/60); # rad/s

a=2*Vmin/(np.pi*d**2)*(rc -1)   # m
L=a*3   # m

def V_t(t):
    V=Vmin + np.pi*pow(d,2)/4*(L+a-a*np.cos(omega*t) - np.sqrt(pow(L,2) - pow(a*np.sin(omega*t),2)))
    return V
# -----------------------------------------------------------------------------
# composition calculations
phi=1

T0=300  # K
P0=ct.one_atm # Pa

gas.TP=[T0, P0]

gas.set_equivalence_ratio(phi, 'C7H16', 'O2:1.0, N2:3.76')
M=gas.molecular_weights     # kg/kmol
Rm=ct.gas_constant    # kJ/kmol K

Y0=gas.Y
mass=gas.density*Vmax

#------------------------------------------------------------------------------
def dy_Vt(t,y,gas):
    dydt = np.zeros(y.size)
    T = y[0]
    V = y[1]
    Y = y[2:]  

    rho = mass / V
    gas.TDY = [T, rho, Y]
    p=gas.P
    wdot_s=gas.net_production_rates     # kmol/(m3 s)
    cv=gas.cv_mass  # kJ/(kg K)
    u_s=gas.standard_int_energies_RT * Rm * T /M  # kJ/kg
    
    dYdt= M / rho * wdot_s
    dVdt= np.pi * pow(d,2)/4 * ((pow(a,2)*omega*np.sin(omega*t)*np.cos(omega*t))/(np.sqrt(pow(L,2) - pow(a*np.sin(omega*t),2))) + a*omega*np.sin(omega*t))
    dTdt= 1/cv * ( -(u_s * dYdt ).sum() - p/mass *dVdt)
    
    dydt[0]=dTdt
    dydt[1]=dVdt
    dydt[2:]=dYdt
    return dydt

#  ----------------------------------------------------------------------------
# solver code
t0 = -np.pi/omega
t1 =t0 +2*np.pi/omega

y0 = np.zeros(ns+2)
y0[0] = T0
y0[1] = V_t(t0)
y0[2:] = Y0


solver = ode(dy_Vt).set_integrator('vode', method='bdf', order=15)
solver.set_initial_value(y0, t0).set_f_params(gas)

gas()
    
tsol = []
Tsol = []
Vsol= []

N2sol = []
Fsol =  []
O2sol = []
H2Osol = []
CO2sol = []
Psol=[]
while solver.successful() and np.any(solver.t < t1):
    solver.integrate(t1,step=True)
    tsol.append(solver.t)
    
    Tsol.append(solver.y[0])
    Vsol.append(solver.y[1])
    
    N2sol.append(solver.y[iN2+2])
    Fsol.append(solver.y[iF+2])
    O2sol.append(solver.y[iO2+2])
    H2Osol.append(solver.y[iH2O+2])
    CO2sol.append(solver.y[iCO2+2])
    Psol.append(gas.P)
    
    
    
    Ysol=np.column_stack((O2sol,Fsol,CO2sol,H2Osol,N2sol))
    Ysol1=np.column_stack((O2sol,Fsol,CO2sol,H2Osol))
gas()
# -----------------------------------------------------------------------------
# plotting
Psol1=np.zeros(np.size(Tsol))
i=0
for T in Tsol:
    Psol1[i]= mass/Vsol[i] * sum(Ysol[i]*Rm/M)*T
    i=i+1

deg=[]
for t in tsol:
    deg.append(t*omega/(2*np.pi)*360)


#ttest=np.linspace(t0,t1,10)

plt.plot(deg,Vsol,label=r'$V( \theta )$')
plt.plot([-180, 180],[Vmin, Vmin], label='V_TDC')
plt.plot([-180,180],[Vmax,Vmax], label='V_BDC') #,ttest,V_t(ttest),'-')
plt.xlabel('Crank angle degree')
plt.ylabel('Volume [m3]')
plt.title('Volume change over 1 cycle')
plt.legend()
plt.savefig('volume.pdf')
plt.show()

plt.plot(deg,Tsol,label='Temperature')
plt.xlabel('Crank angle degree')
plt.ylabel('Temperature [K]')
plt.title('Temperature change over 1 cycle')
plt.legend()
plt.savefig('tempB.pdf')
plt.show()

plt.plot(deg,np.array(Psol1)/ct.one_atm,label='Pressure')
plt.xlabel('Crank angle degree')
plt.ylabel('Pressure [atm]')
plt.title('Pressure change over 1 cycle')
plt.legend()
plt.savefig('pressure.pdf')


plt.show()
plt.plot(deg,Fsol,label='Fuel')
plt.plot(deg,O2sol,label='O2')
plt.plot(deg,H2Osol,label='H2O')
plt.plot(deg,CO2sol,label='CO2')
plt.legend()
plt.xlabel('Crank angle degree')
plt.ylabel('Species massfractions')
plt.title('Gas composition change over 1 cycle')
plt.savefig('compositionB.pdf')
plt.show()
