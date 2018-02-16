# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 20:37:42 2017

@author: Gabriel
"""
import numpy as np, matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.integrate import ode

# --------------------------------------------------
           
rho = {'water': 1000, 'air': 1.2, 'glass': 2560, 'oil': 915}# kg/m3 Densities
my = {'water': 0.001, 'air': 1.8e-5, 'oil': 0.0574}         # Pa s Viscosities 

dsphere=0.0005      # m diameter of glass sphere
dbubble=0.15e-2     # m diameter of air bubble

g=9.81          # m/s2 acceleration due to gravity 
sigma=71.95e-3  # N/m surface tension water-air 

kappa=my['air']/my['water'] # viscosity quota for hadamard drag law 

# read experimental data
expData=np.genfromtxt('exptData.dat')

texp=expData[:,0]
vexp=expData[:,1]

# -------------------------------------------------
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
    
    
    

#%%
# ---------------------------------------------------------------------------
# Vector generation glass sphere
tint=[0.,0.4]
dt=(tint[1]-tint[0])/1000

#cam=cam(psphere)
Cam=1
f = lambda vp: dvdt(vp,Cam)
t,y = fwe(f,tint,dt)

# Plotting

plt.figure(1)
plotSol(texp,vexp,'Experimental Data','--')

plotSol(t,y,'Schiller-Naumann Drag Model')

f = lambda vp: dvdt(vp,Cam,psphere,'stokes')
t,y = fwe(f,[0,0.4],0.4/1000)

plotSol(t,y,'Stokes regime Drag Model')


f = lambda vp: dvdt(vp,Cam,psphere,'stokes-comp')
t,y = fwe(f,[0,0.4],0.4/1000)

plotSol(t,y,'$C_D = 1/2 + 24/Re_p $')

plt.xlim(0,0.17*1000)

plt.title(r'Velocity profile for falling glass sphere, $u_f = 0$')
plt.legend()
plt.savefig('images/velo-glass.pdf')



# Figure for density ratio vs response time
rhoratio=np.linspace(0.1,10)
tau_vec=rhoratio*dsphere**2/18/(my['water']/rho['water'])

plt.figure(2)
plt.plot(rhoratio,tau_vec*1000,label=r'$ \tau \left( \frac{\rho_p}{\rho_f} \right) $')
plt.xlabel('Density ratio')
plt.ylabel('Particle response time [ms]')
plt.title('Effect of density ratio on response time')
plt.legend()
plt.savefig('images/tau-rho.pdf')


# Generate numerical demo

tint=[0.,0.3]
dt=(tint[1]-tint[0])/100

def odefunc(t,v,cam,p,draglaw):
    return dvdt(v,cam,p,draglaw)

solver = ode(odefunc).set_integrator('dopri5')
solver.set_initial_value(y0,tint[0]).set_f_params(Cam,psphere,'schiller-naumann')

vsol=[]
tsol=[]

while solver.successful() and np.any(solver.t < tint[1]):
    solver.integrate(solver.t + dt)
    tsol.append(solver.t)
    vsol.append(solver.y[0])

vsol=-np.array(vsol)
tsol=np.array(tsol)*1000


fig=plt.figure(3,figsize=(9,4))
ax1=fig.add_subplot(121)
ax2=fig.add_subplot(122)

ax1.plot(tsol,vsol,'--',label='Runge Kutta 4th order')
ax2.plot(tsol,vsol,'--',label='Runge Kutta 4th order')

f = lambda vp: dvdt(vp,Cam)
t,y = fwe(f,tint,dt)
ax1.plot(t,y,label='Forward Euler')
ax2.plot(t,y,label='Forward Euler')

f = lambda vp: dvdt(vp,Cam)
t,y = bwe(f,tint,dt,1e-10)
ax1.plot(t,y,label='Backward Euler')
ax2.plot(t,y,label='Backward Euler')

plt.suptitle('--- Different numerical schemes, dt=' + str(dt*1000) +'ms ---')

ax1.set_title('Entire profile')
ax2.set_title('Zoomed in')

dt=(tint[1]-tint[0])/400

f = lambda vp: dvdt(vp,Cam)
t,y = fwe(f,tint,dt)

ax1.plot(t,y,label='Forward Euler, dt=' + str(dt*1000) + 'ms')
ax2.plot(t,y,label='Forward Euler, dt=' + str(dt*1000) + 'ms')

ax1.set_xlim(0,0.07*1000)
ax1.set_xlabel('Time [ms]')
ax1.set_ylabel('Velocity [m/s]')
ax1.legend()
ax2.set_xlabel('Time [ms]')
ax2.set_ylabel(' ')
ax2.legend()


ax2.set_xlim(0.005*1000,0.03*1000)
ax2.set_ylim(0.03,0.07)
plt.savefig('images/num.pdf')

#------------------------------------------------------------------------------ 
# Falling in air
tint=[0.,1.5]
dt=(tint[1]-tint[0])/1000

plt.figure(4)

f = lambda vp: dvdt(vp,Cam)
t,y = fwe(f,tint,dt)
taup,vterm=findTau(t,y)

ybar=[]
ybar.append(y[-1])
ybar.append(y[-1])
plt.plot(np.array([0.,1.5])*1000,ybar,'--',label=r'Falling in water, $\tau_p$=' + '{:f}'.format(taup) + 'ms')



pair={'rhof': rho['air'], 'rhop': 2560, 'myf': my['air'], 'dp': 0.0005}

f = lambda vp: dvdt(vp,Cam,pair)
t,y = fwe(f,tint,dt)
taup,vterm=findTau(t,y)
print(taup)

plotSol(t,y,r'Falling in air, $\tau_p =$' + '{:f}'.format(taup) + 'ms')
plt.legend()
plt.title('Falling through different fluids')
plt.savefig('images/velo-air.pdf')

#%%
# -----------------------------------------------------------------------------
# Start bubble calcs
pbubble={'rhof': rho['water'], 'rhop': rho['air'], 'myf': my['water'], 'dp': dbubble}

tint=[0.,0.03]
dt=(tint[1]-tint[0])/1000

#Cam=cam(pbubble)
Cam=cam(pbubble)
f = lambda vp: dvdt(vp,Cam,pbubble,'tran-cong')
t,y = fwe(f,tint,dt)

plt.figure(5)
taup,vtau=findTau(t,y)
p,po=plotSol(t,y,r'T-C drag, $v_{term}=$' + '{:f}'.format(y[-1]) +' m/s')
po[0].set_label(r'$\tau_p = $' + '{:f}'.format(taup) + 'ms')

Cam=cam(pbubble)
f = lambda vp: dvdt(vp,Cam,pbubble,'schiller-naumann')
t,y = fwe(f,tint,dt)
taup,vtau=findTau(t,y)

p,po=plotSol(t,y,r'S-N drag, $v_{term}=$' + '{:f}'.format(y[-1]) +' m/s')
po[0].set_label(r'$\tau_p = $' + '{:f}'.format(taup) + 'ms')

Cam=cam(pbubble)
f = lambda vp: dvdt(vp,Cam,pbubble,'hadamard')
t,y = fwe(f,tint,dt)
taup,vtau=findTau(t,y)

p,po=plotSol(t,y,r'H-D drag, $v_{term}=$' + '{:f}'.format(y[-1]) +' m/s')
po[0].set_label(r'$\tau_p = $' + '{:f}'.format(taup) + 'ms')

plt.title('Rising bubble, Added mass force included')
plt.legend()
plt.savefig('images/velo-bubble.pdf')

#------------------------------------------------------------------------------
# Bubble in oil
pbubble['rhof']=rho['oil']
pbubble['myf']=my['oil']



tint=[0.,0.005]
dt=(tint[1]-tint[0])/1000


#Cam=cam(pbubble)
Cam=cam(pbubble)
f = lambda vp: dvdt(vp,Cam,pbubble,'tran-cong')
t,y = fwe(f,tint,dt)

plt.figure(6)

taup,vtau=findTau(t,y)

p,po=plotSol(t,y,r'T-C drag, $v_{term}=$' + '{:f}'.format(y[-1]) +' m/s')
po[0].set_label(r'$\tau_p = $' + '{:f}'.format(taup*1000) + ' ms')

plt.title('Bubble rising in oil')
plt.legend()
plt.savefig('images/oil.pdf')



#%%

plt.figure(7, figsize=(9,15))
tint=[0.,0.4]
dt=(tint[1]-tint[0])/1000

Cam=cam(psphere)
f = lambda vp: dvdt(vp,Cam)
t,y = fwe(f,tint,dt)

plt.subplot(411)
forcePlot(t,-y,psphere,'Comparison of Forces - falling in water')
plt.xlim(0,170)
plt.legend(loc=(0.6,0.1))
ax=plt.gca()
ax.axes.set_xlabel(' ')

plt.subplot(412)
tint=[0.,1.5]
dt=(tint[1]-tint[0])/1000

Cam=cam(pair)
f = lambda vp: dvdt(vp,Cam,pair)
t,y = fwe(f,tint,dt)

forcePlot(t,-y,pair,'Falling in air')
#plt.xlim(0,14000)
plt.legend(loc=(0.6,0.1))
ax=plt.gca()
ax.axes.set_xlabel(' ')
plt.xlim(0,1450)

plt.subplot(413)
tint=[0.,0.03]
dt=(tint[1]-tint[0])/1000


pbubble={'rhof': rho['water'], 'rhop': rho['air'], 'myf': my['water'], 'dp': dbubble}
Cam=cam(pbubble)
f = lambda vp: dvdt(vp,Cam,pbubble,'tran-cong')
t,y = fwe(f,tint,dt)

forcePlot(t,y,pbubble,' Bubble rising in water','tran-cong')
#plt.xlim(0,14000)
plt.legend(loc=(0.59,0.55))
ax=plt.gca()
ax.axes.set_xlabel(' ')
plt.xlim(0,30)

plt.subplot(414)
tint=[0.,0.005]
dt=(tint[1]-tint[0])/1000


pbubble={'rhof': rho['oil'], 'rhop': rho['air'], 'myf': my['oil'], 'dp': dbubble}
Cam=cam(pbubble)
f = lambda vp: dvdt(vp,Cam,pbubble,'tran-cong')
t,y = fwe(f,tint,dt)

forcePlot(t,y,pbubble,'Bubble  rising in oil','tran-cong')
#plt.xlim(0,14000)
plt.legend(loc=(0.59,0.55))
plt.xlim(0,5)
plt.savefig('images/force.pdf')

# -----------------------------------------------------------------------------
# Weber shit

pbubble={'rhof': rho['water'], 'rhop': rho['air'], 'myf': my['water'], 'dp': dbubble}
Cam=cam(pbubble)
f = lambda vp: dvdt(vp,Cam,pbubble,'tran-cong')
t,y = fwe(f,tint,dt)


We=pbubble['rhof']*y**2*dbubble/sigma

Eo=g*np.abs(rho['water']-rho['air'])*dbubble**2/sigma
Re=Rep(y,pbubble)

Mo=np.log(Eo*We**2/pow(Re,4))

hdquota=(3*kappa+2)/(3*kappa+3)

print('HD quota: ' + str(hdquota))



# ---------------------------------------------------
psand={'rhof': 0.74317, 'rhop': 2650, 'myf': 2.6008e-5, 'dp': 0.95e-3}
plt.figure(8)
tint=[0.,4]
dt=(tint[1]-tint[0])/20000

Cam=cam(psand)
f = lambda vp: dvdt(vp,Cam,psand,'schiller-naumann')
t,y = fwe(f,tint,dt)
taup,vtau=findTau(t,y)

p,po=plotSol(t,y,r'S-N drag, $v_{term}=$' + '{:f}'.format(y[-1]) +' m/s')
po[0].set_label(r'$\tau_p = $' + '{:f}'.format(taup) + 'ms')
plt.legend()

plt.figure(9)

forcePlot(t,-y,psand,'sand falling')
plt.legend(loc=(0.59,0.55))


