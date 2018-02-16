import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt

def ode1(t, y, mu):
    # t - Independet variable
    # y - dependent variable
    # mu - parameter
    y1 = y[0]
    y2 = y[1]
    dy1dt = y2
    dy2dt = mu*(1 - y1*y1)*y2 - y1
    return [dy1dt, dy2dt]



y0 = [2, 0] #Initial Values
t0 = 0 # Start Time
mu = 1 # Parameter
sol = ode(ode1).set_integrator('vode', method='bdf')
sol.set_initial_value(y0, t0).set_f_params(mu)

t1 = 20
dt = 1e-3


y1Sol = []
y2Sol = []
tSol = []
while sol.successful() and sol.t < t1:
    sol.integrate(sol.t +dt)
    tSol.append(sol.t)
    y1Sol.append(sol.y[0])
    y2Sol.append(sol.y[1])

plt.plot(tSol,y1Sol, tSol, y2Sol)
plt.show()
