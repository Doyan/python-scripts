# -*- coding: utf-8 -*-
"""
Created on Mon 12 February 2018

@author: Gabriel
"""
import numpy as np, matplotlib.pyplot as plt

#------------------------------------------------------------------------------
fnam="data.csv" # data from handbook of chemistry and physics at 1bar

file=open(fnam,'r') 

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        pass
    
    return False 

rows=file.readlines()

data={}
for idx,row in enumerate(rows):
    row=row.translate({ord(c): None for c in '"'})
    data[idx]=row.split(',')
    for n, item in enumerate(data[idx]):
        if is_number(item):
            data[idx][n]=float(item)
file.close()

# -------Luftegenskaper------------------------
# data from handbook of chemistry and physics at 1bar
Trange=[data[0][1],data[1][1]]
murange=[data[0][8],data[1][8]]
rhorange=[data[0][2],data[1][2]]

T_air=200 + 273.15

mu=np.interp(T_air,Trange,murange)*1e-6
rho=np.interp(T_air,Trange,rhorange)

#------------Geometry--------------------
L=7.4 
W=5.9
D=L*W/(L+W) # Hydraulic diameter

Hbed=1.0 # @ bubbling fluidisation
Abed=L*W
Vbed=Abed*Hbed

N_or=28*11*3

#----- Partikelegenskaper ----------------------------------------------------
dp=0.95e-3
rhop=2.65e3
g=9.82
vT=7.425 # from mpflow assignment 1 code
phi=0.9 # sphericity
shift=1.34




epsmf=(14*phi)**(-1/3) # voidage at minimum fluidisation

BDcrit=(rhop-rho)*dp**2
Ar=rho*(rhop-rho)*g*dp**3/(mu**2)

Remf=np.sqrt(27.2**2 + 0.0408*Ar) - 27.2
Rec=1.24*Ar**0.45

# Estimation of suitable velocity for bubbling fluidisation
Ustar_mf=Remf/(Ar**(1/3))
Ustar_c=Rec/(Ar**(1/3))
Ustar_mid=np.sqrt(Ustar_mf*Ustar_c)*shift

Ububbling=Ustar_mid*Ar**(1/3)*mu/(rho*dp)
U_mf=Remf*mu/(rho*dp)
U_c=Rec*mu/(rho*dp)
uc=3.0*(rhop*dp)**(1/2)-0.17
Ums=U_mf+0.07*np.sqrt(g*D)
# Estimation of bed voidage at fluidised state - particulate regime.
ReT=vT*rho*dp/mu
n=4.4*ReT**(-0.1)

ui=vT-10**(dp/D)

eps=(Ububbling/ui)**(1/n)

# Correlation for estimation of overall voidage in bubbling fluidisation
epsbar=Ububbling*epsmf/(1.05*Ububbling*epsmf+(1-epsmf)*U_mf)


# Estimation of total mass of sand
Vsand=(1-epsbar)*Vbed
Hbedmf=Vsand/((1-epsmf)*L*W)
msand=Vsand*rhop

Hres=1.0

# --------------- Plot of Re and Ar^1/3---------------------

def ReC(Ar):
    return 1.24*Ar**(0.45)

def ReMf(Ar):
    return np.sqrt(27.2**2 + 0.0408*Ar) - 27.2

def Umid(Ar):
    return np.sqrt(ReMf(Ar)*ReC(Ar)*Ar**(-2/3))





Arr=np.linspace(0.001,100**(3),100)
Ar13=Arr**(1/3)


plt.figure(1)
plt.loglog(Ar13,ReMf(Arr)/Ar13,label='$U^*_{mf}$')
plt.loglog(Ar13,ReC(Arr)/Ar13,label='$U^*_{c}$')
plt.loglog([Ar**(1/3),Ar**(1/3)],[1e-4,1e3],'--k',label='$Current\; Ar$')
plt.loglog(Ar13,Umid(Arr),'--m',label='$U^*_{mid}$')

plt.legend(loc='best')
plt.xlabel('$Ar^{1/3}$')
plt.ylabel('$Re(Ar)/Ar^{1/3}$')

plt.figure(2)
plt.plot(Ar13,ReMf(Arr)/Ar13,label='$U^*_{mf}$')
plt.plot(Ar13,ReC(Arr)/Ar13,label='$U^*_{c}$')
plt.plot([Ar**(1/3),Ar**(1/3)],[1e-4,1e3],'--k',label='$Current\; Ar$')
plt.plot(Ar13,Umid(Arr),'--m',label='$U^*_{mid}$')




plt.legend(loc='best')
plt.xlabel('$Ar^{1/3}$')
plt.ylabel('$Re(Ar)/Ar^{1/3}$')
plt.xlim(0,40)
plt.ylim(0,5)


#-- Uppskattning av relaxationstid------------------------

Rep=Remf

Re_gas=rho*Ububbling*D/mu

alpha=0.63
alphag=1-alpha
tau_stokes=rhop*dp**2/(18*mu)
A=alphag**4.14
nu_rs=0.5*(A-0.06*Rep+np.sqrt((0.06*Rep)**2+0.12*Rep*(2*0.8*alphag**(1.28)-A)+A**2))
Cd=(0.63+4.8/np.sqrt(Rep/nu_rs))**2

f=Cd*Rep*alphag/(nu_rs**2)

tau=tau_stokes/f

#------------- Tuning syamlal o brien
epsMf=epsmf
epsmf=(1-0.63)

U_Mf=dp**2*(rhop-rho)*g/(150*mu)*epsmf**3/(1-epsmf)

Re_ts=(((4.8**2+2.52*(4*Ar/3)**(1/2))**(1/2)-4.8)/1.26)**2
Re_t=rho*U_Mf*dp/(mu*epsmf)

V_r=Re_t/Re_ts
A=epsmf**4.14

K=0.06*Re_t

B= (((2*V_r-A+K**2)**2-K**2-A**2)/(2*K))/2

C2=B/(epsmf**1.28)

C1=1.28+np.log(C2)/np.log(0.85)


epsmf=epsMf
#-------------Printing---------------------------------
print('--------------Air-data--------------')
print('rho= '+ str(rho) + ' kg/m3')
print('mu= '+ str(mu) + ' kg/m3')
print(' ')


print('-------------Particle properties----------------')
print('rho= {:} kg/m3'.format(rhop))
print('d_p= {:.3} mm\n'.format(dp*1000))

print("Geldart B-D = {:.2E} kg/m -> group D".format(BDcrit))
print("vT/Umf = {:.3} -> somewhat large particle".format(vT/U_mf))
print("Ar = {:.3E}".format(Ar))
print("tau_p = {:.2E}s\n".format(tau))

print("Re_p @ minimum fluidisation = {:.3}".format(Remf))
print("Re_p @ terminal velocity = {:.3}\n".format(ReT))

print('----------------Bed Properties------------')
print('total mass of sand = {:.3} kg\n'.format(msand))
print('air velocity @ minimum fluidisation = {:.3}'.format(U_mf))
print('air velocity @ turbulent fluidisation = {:.3}'.format(U_c))
print('air velocity @ bubbling fluidisation = {:.3}\n'.format(Ububbling))

print('air velocity @ slugging fluidisation = {:.3}'.format(Ums))
print('condition for slugging fluidisation: H > {:.3}\n'.format(3.5*D*(1-1/np.sqrt(N_or))))


print('voidage @ minimum fluidisation = {:.3} -> alpha_mf = {:.3}'.format(epsmf,(1-epsmf)))
print('voidage @ bubbling fluidisation = {:.3} -> alpha = {:.3}\n'.format(eps,(1-eps)))
print('Height @ fluidisation = {:.3} m'.format(Hbed))
print('Height @ minimum fluidisation = {:.3} m\n'.format(Hbedmf))
#print('resulting Height @ fluidisation = {:.3} m'.format(Hres))

print('C1 = {:.3}, C2 = {:.3}'.format(C1,C2))






