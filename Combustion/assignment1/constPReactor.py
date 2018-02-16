from cantera import *
from numpy import *
from scipy import integrate

gas = Solution("h2o2_Mueller_Chen.cti")
ns  = gas.n_species
iN2  = gas.species_index('N2')
iO2  = gas.species_index('O2')
iH2O = gas.species_index('H2O')
iF = gas.species_index('H2')

gas.transport_model = 'Mix'
Le = gas.thermal_conductivity / (gas.density * gas.cp_mass * gas.mix_diff_coeffs_mass)

#------------------------------------------------------------------------------
# initila state of the mixture 
T0 = 700.
p0 = one_atm

X0   = zeros(ns)
X0[iO2] = 21.
X0[iN2] = 79. 
X0[iF]  = 42.     # H2

gas.TPX = [T0, p0, X0]

Y0 = gas.Y


# time of integration
t_start = 0.0;
t_end   = 1.e-5;
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# dy/dt zur Integration eines homogenen Reaktors
# bei konstantem Druck:
# t   - aktuelle Zeit
# y   = [T, Y_s] - aktuelle Zustand der Integrationsvariablen
# p   - Druck
# gas - Cantera Objekt 
#------------------------------------------------------------------------------
def dy_p(t,y,p,gas):
  dydt = zeros(y.size)
  T = y[0]
  Y = y[1:]
  gas.TPY = [T, p, Y]
  rho = gas.density
  wdot_s = gas.net_production_rates
  M_s    = gas.molecular_weights
  dydt[1:] = wdot_s * M_s / rho
  h_s    = gas.standard_enthalpies_RT * gas_constant  * T
  cp     = gas.cp_mass
  dydt[0]  = - (h_s * wdot_s).sum() / rho / cp
  return dydt 


#------------------------------------------------------------------------------


# Initialisierung des Zeitintegrators, hier BDF wie Matlab ode15s
y0     = zeros(ns+1)
y0[0]  = T0
y0[1:] = Y0
ode15s = integrate.ode(dy_p)
ode15s.set_integrator('vode', method='bdf', order=15, nsteps=3000)
ode15s.set_f_params(p0,gas)
ode15s.set_initial_value(y0, t_start)


nt = 400
time = linspace(t_start, t_end, nt)

Y = zeros([nt,ns])
Y[0] = Y0

T=zeros(nt)
T[0] = T0

for t in range(1,nt):
  yn = ode15s.integrate(time[t])
  if (t%10 == 0):
    print t, ": time = ", time[t], "s", " T = ", yn[0]

  T[t] = yn[0]

import pylab as py

#py.figure(1)
py.semilogx(time,T)
py.ylabel('T [K]')
py.xlabel('time')
py.text(time.max()*0.1,T.max()*0.9,'T = ')

py.show()





