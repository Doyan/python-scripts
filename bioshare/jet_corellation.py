# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 16:29:19 2019
Compares jet penetration depths produced by correllations from litterature.

@author: Gabriel Gustafsson
"""
import numpy as np
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------

rho_steam= 0.45

rho_sand=2650
g=9.81

dp = 0.9e-3

eps_j = 0.01



def Salatini(dp,vj,do,rhop=rho_sand,rhof=rho_steam):
    Fr  = rhof/(rhop-rhof)*vj**2/(g*do)
    return 61.2*Fr**0.4*(1/(1+dp/do)**5.6)*(rhof/rhop)**0.4






Vj = np.linspace(0,0.5)




Lj=Salatini(0.9e-3,15.0,0.04)

plt.figure(1,(10,6))


do=0.04
for do in [0.01,0.02,0.04,0.08,0.1]:
    vj = Vj/(do**2*np.pi/4)
    Lj = Salatini(dp,vj,do)*do

    plt.plot(Vj,Lj,label='$d_o={}$'.format(do))

plt.grid()
plt.legend()
plt.ylabel('$L_j/d_o$')
plt.xlabel('$Vg$ [m3/s]')






# ================================

def boundary(x):
    L = 0.05 # Major fitting param
    
    rho_g = rho_steam
    rho_s = rho_sand
    
    epss0 = 0.5 # check later
    r0 = 0.04
    u0= 50
    theta=17*np.pi/180 #assumed check later
    CD = 0.44 # assumed check later
    
    tol= 1e-4
    
    
    VL=10
    HL=100
    
    
    #b = 0.7616680
    
      
    b0 = x*np.tan(theta) + r0
    um = 1/0.1335*u0*r0/b0
    eps_s=epss0*(1-np.exp(-x/L))
    ut = (4/(3*CD)*rho_s*g*dp/rho_g)**(1/2)
    ub = 10 * ut *  eps_s**3
    
    # left side of momentum balance
    VL = rho_g*u0**2*r0**2
    
    I = [0,1]
    check = 1
    while check > tol:
        b = (I[-1] - I[0])/2 + I[0]
        
        # right side of momentum balance
        xi=b/b0
        
        uc = um*(1-xi**(1.5))**2
        du = uc -ub
        print(du)
    
        HL1 = rho_g * um**2 * b0**2 * (1/2*xi**2 - 8/7*xi**(7/2) + 6/5*xi**5 - 8/13*xi**(13/2)+ 1/8*xi**8) 
        HL2 = -2*rho_g * um * du * b0**2 * (1/2*xi**2 - 4/7*xi**(7/2) + 1/5*xi**5) + 1/2*rho_g*du**2*b**2
        HL3 = (rho_s*eps_s + rho_g*(1-eps_s)) * ub**2 * b**2 *( 1/4*(np.exp(x/(2*L)) - 1)**2 + 1/3*b*(np.exp(x/(2*L)) - 1 ))
        HL = HL1 + HL2 + HL3

        diff = VL - HL
        check = I[-1] - I[0]
        print(b,diff,I,check)
        if diff < 0:
            I[-1] = b
        else:
            I[0] = b
    
    b1  = b*np.exp(-x/(2*L))
    return b, b1, b0, 

B = np.vectorize(boundary)

plt.figure(2,(12,7))
x = np.linspace(0,0.2)
b, b1, b0 = B(x)
plt.plot(x,b,x,b1,x,b0)
plt.ylim(0,0.1)




