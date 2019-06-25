# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 20:37:42 2017
Falling sphere in quiesent fluid
@author: Gabriel
"""
import numpy as np, matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import fsolve
from scipy.integrate import ode

# --------------------------------------------------
           
rho = {'water': 1000, 'air': 1.2, 'glass': 2560, 'oil': 915}# kg/m3 Densities
my = {'water': 0.001, 'air': 1.8e-5, 'oil': 0.0574,'dirty': 0.0125}         # Pa s Viscosities 

dsphere=0.0005      # m diameter of glass sphere
dbubble=0.00488     # m diameter of air bubble

g=9.81          # m/s2 acceleration due to gravity 
#sigma=71.95e-3  # N/m surface tension water-air 
sigma = 0.073

kappa=my['air']/my['water'] # viscosity quota for hadamard drag law 

# ---------------------------------------------
# Selection of properties

psphere={'rhof': 1000, 'rhop': 2560, 'myf': 0.001, 'dp': 0.0005}

# -------------------------------------------------------
# Some helpful functions....
 
# Particle Reynolds number for standalone calc
def Rep(vp,p=psphere):
    return vp*p['rhof']*p['dp']/p['myf']

# Drag correlations and drag force
def Fd(vp,p,draglaw):
    rhof=p['rhof']
    myf=p['myf']
    dp=p['dp']
        
    Rep=np.abs(vp)*rhof*dp/myf
    Eo=g*np.abs(p['rhof']-p['rhop'])*dp**2/sigma
    
    
    cd_stokes=24/Rep
    
    if draglaw == 'schiller-naumann': # from lecture slides
        cd=cd_stokes*(1+0.15*Rep**0.687)
    elif draglaw == 'tran-cong': # from article suggested in book
        D=3/2   # given ratio, c is assumed 1 - only motion i z-direction so circular projection 
        f=D*(1+0.15*(D*Rep)**0.687) + (0.0175*D**2)/(1 +4.25e4*(D*Rep)**(-1.16))
        cd=cd_stokes*f
    elif draglaw == 'stokes-comp': # correction suggested in assigment text
        cd=1/2 + cd_stokes
    elif draglaw == 'stokes':
        cd=cd_stokes
    elif draglaw == 'hadamard':
        D=3/2   # given ratio, c is assumed 1 - only motion i z-direction so circular projection 
        f=D*(1+0.15*(D*Rep)**0.687) + (0.0175*D**2)/(1 +4.25e4*(D*Rep)**(-1.16))
        cd=cd_stokes*f*( (2/3 + kappa)/(1 + kappa) )
        #cd=cd_stokes*(1+0.15*Rep**0.687)*( (2/3 + kappa)/(1 + kappa) )
    elif draglaw == 'toniyama':
        cd = max([min([16/Rep*(1+0.15*Rep**0.687),cd_stokes*2]),(8/3)*(Eo/(Eo+4))])
    else:
        raise ValueError('Invalid Drag law')
        
    return 1/2 * rhof * dp**2 * np.pi/4 * cd * np.abs(-vp)*-vp

# ---------------------------------------------------------------    
# Coefficient for added mass force
    
def cam(p):
    return 1+p['rhof']/(2*p['rhop'])
  
# --------------------------------------------------------------
# # Newtons second law - Force balance

def dvdt(vp,cam,p=psphere,draglaw='schiller-naumann'):
    rhof=p['rhof']
    rhop=p['rhop']
    dp=p['dp']
    
    mp= np.pi*dp**3/6 * rhop # kg, Mass of sphere
    
    Fg= -mp*g           # Force due to gravity
    Fb= mp/rhop*rhof*g  # Force due to bouyancy
     
    return (Fd(vp,p,draglaw) + Fb + Fg )/(mp*cam)

# ---------------------------------------------------------------
# Numerics
y0=0.0000000000001
# Forward Euler solver
def fwe(dvdt,tint,dt,y0=y0):
    y=[y0]
    t=[tint[0]]
    i=0
    while t[i] < tint[1]:
        y.append(y[i] + dt*dvdt(y[i]))
        t.append(t[i]+dt)
        i+=1
    return np.array(t)*1000, np.abs(np.array(y))

# Backward Euler solver
def bwe(dvdt,tint,dt,tol=1e-6,y0=y0):
    y=[y0]
    t=[tint[0]]
    i=0
    while t[i] < tint[1]:
        yg=y[i]
        y1=y[i]+dt*dvdt(y[i])
        while np.abs(y1-yg) > tol:
            yg=y1
            y1=y[i]+dt*dvdt(yg);
        
        y.append(y1)
        t.append(t[i]+dt)
        i+=1
    return np.array(t)*1000, np.abs(np.array(y))



# Helpers 
def findTau(t,y):
    vterm=y[-1]
    vtau=(np.e-1)/np.e*vterm
    tau_p=np.interp(vtau,y,t)
    return tau_p, vtau

def plotSol(t,y,Label,linestyle='-'):
    tau_p, vtau = findTau(t,y)
    
    
    plt.xlabel('Time [ms]')
    plt.ylabel('Particle Velocity [m/s]')
    p=plt.plot(t,y,linestyle,label=Label,)
    po=plt.plot(tau_p,vtau,'o',color=p[0].get_color())
    
    
    
    return p,po

def forcePlot(t,v,p,title,draglaw='schiller-naumann'):
    rhof=p['rhof']
    rhop=p['rhop']
    dp=p['dp']
    
    mp= np.pi*dp**3/6 * rhop # kg, Mass of sphere
    
    Fg= -mp*g*1e6           # Force due to gravity
    Fb= mp/rhop*rhof*g*1e6  # Force due to bouyancy
    
    
    Cam=cam(p)
    Fam=1/2 * mp * rhof/rhop * -dvdt(v,Cam,p,draglaw)*1e6
    Fddata=Fd(v,p,draglaw)*1e6
    
    plt.plot([t[0], t[-1]],[Fg, Fg],label='Gravity, Fg='+ '{:f}'.format(Fg) + ' nN' )
    plt.plot([t[0], t[-1]],[Fb, Fb],label='Bouyancy, Fb='+ '{:f}'.format(Fb) + ' nN' )
    plt.plot(t,Fddata,label=r'Drag, $Fd_{term}$=' + '{:f}'.format(Fddata[-1]) + ' nN')
    plt.plot(t,Fam,label='Added mass force')
    plt.xlabel('Time [ms]')
    plt.ylabel(r'Force [$\mu$N]')
    plt.title(title)
    return
    
    
    
#%%

psphere={'rhof': 1000, 'rhop': 2560, 'myf': 0.001, 'dp': 0.0005}
    
# Example usage
    
# Vector generation glass sphere
tint=[0.,0.4]
dt=(tint[1]-tint[0])/1000

#cam=cam(psphere)
Cam=1
f = lambda vp: dvdt(vp,Cam)
t,y = fwe(f,tint,dt)

# Plotting

plt.figure(1)
#plotSol(texp,vexp,'Experimental Data','--')

plotSol(t,y,'Schiller-Naumann Drag Model')

f = lambda vp: dvdt(vp,Cam,psphere,'stokes')
t,y = fwe(f,[0,0.4],0.4/1000)

plotSol(t,y,'Stokes regime Drag Model')


#%%

plt.figure(2)

tint=[0.,4.0]
dt=(tint[1]-tint[0])/1000

pcoal = {'rhof': 1.2,'rhop':1470,'myf':my['air'],'dp':130e-6}

Cam =cam(pcoal)
f = lambda vp: dvdt(vp,Cam,p=pcoal)
t,y = fwe(f,tint,dt)

plotSol(t,y,'Schiller-Naumann Drag Model')
taup,vterm=findTau(t,y)
print(taup,vterm)

psand = {'rhof': 0.45,'rhop':2600,'myf':my['air'],'dp':0.9e-3}

Cam =cam(psand)
f = lambda vp: dvdt(vp,Cam,p=psand)
t,y = fwe(f,tint,dt)

plotSol(t,y,'Schiller-Naumann Drag Model')
taup,vterm=findTau(t,y)
print(taup,vterm)



    
    
    