# -*- coding: utf-8 -*-
"""
Created on Mon Oct 02 01:02:33 2017

Plotting for assignment 2

@author: Gabriel
"""
import numpy as np, matplotlib as mpl, cantera as ct
#from scipy.optimize import fsolve
plt = mpl.pyplot
colors = mpl.colors

gas = ct.Solution('gri30.xml')
iN2  = gas.species_index('N2')
iO2  = gas.species_index('O2')
iH2O = gas.species_index('H2O')
iF = gas.species_index('CH4')
iCO2= gas.species_index('CO2')

# -----------------------------------------------------------------------------
# Data for analytical solution


T_e = 473.0     # K
T_air = 300.0   # K
P = ct.one_atm  # Pa

gas.TPX=[T_air,P,"O2:1.0, N2:3.76"]
Y_O2_2=gas.Y[iO2]
gas.set_equivalence_ratio(1,'CH4','O2:1.0, N2:3.76')

Y_Fst = gas.Y[iF]

z_st = Y_Fst

H_st=gas.enthalpy_mass

gas.HPX=[H_st,P,'CO2:1, H2O:2, N2:6.72']
T_ad=gas.T

Y_st=gas.Y[[iF, iO2, iCO2, iH2O]]

# -----------------------------------------------------------------------------
# Equations

a1=(T_ad - T_air) / z_st       # For Fuel lean mixture 
b1=T_air                        #

a2=(T_e - T_ad)/(1 - z_st)     # For Fuel rich mixture
b2=T_e - a2                     #

a3=(T_e - T_air)                # For unburnt mixture
b3=T_air     

def tempbs(Z):
    condlist= [Z < z_st, Z >= z_st]
    choicelist= [a1*Z + b1, a2*Z + b2]    
    return np.select(condlist,choicelist)

def ybs(Z):
    condlist= [Z < z_st, Z >= z_st]
    
    prodlist= [ Z/z_st * Y_st[2:4], (1-Z)/(1-z_st) * Y_st[2:4]]   
    flist= [0, 1/(1-z_st)*Z - 1/(1-z_st)*z_st ]
    o2list= [-Y_O2_2 / z_st*Z + Y_O2_2, 0]
    
    YF=np.select(condlist,flist)
    YO2=np.select(condlist,o2list)
    YP=np.select(condlist,prodlist)
    return YF, YO2, YP


def T_u(Z):
    return a3*Z + b3

def readDataSheet(fnam):
    data  = np.genfromtxt(fnam) #Input 
    z = data[:,0]
    t = data[:,1]
    yf= data[:,2]
    yo2= data[:,3]
    yco2= data[:,4]
    yh2o= data[:,5]
    return z, t, yf, yo2, yco2, yh2o
# -----------------------------------------------------------------------------
    
Zdata= np.append(np.arange(0, z_st, 0.001), np.arange(z_st, 1, 0.01))
Zdata=Zdata[:, np.newaxis];


YF, YO2, YP=ybs(Zdata)
Tdata= tempbs(Zdata)    


# -----------------------------------------------------------------------------
w=np.array([0.01, 0.1, 1, 10.0, 20.0, 100.0]) # parameter variation width
for val in w:
    fnam= 'width' + str(val) + '.txt'    
    z,t,yf,yo2,yco2,yh2o = readDataSheet(fnam)
    plt.plot(z,t,label='w x ' + str(val))


plt.xlabel('Mixture fraction - Z')
plt.ylabel('Temperature [K]')
plt.title('Parameter variation - width')

plt.legend()
plt.savefig('images/width.pdf')


plt.show()


# parameter variation air to fuel ratio
afr=np.array([0.01, 0.1, 1, 10, 100])
for val in afr:  
    fnam= 'afr' + str(val) + '.txt'    
    z,t,yf,yo2,yco2,yh2o = readDataSheet(fnam)
    plt.plot(z,t,label='afr=' + str(val))

plt.legend()
plt.xlabel('Mixture fraction - Z')
plt.ylabel('Temperature [K]')
plt.title('Parameter variation - air to fuel ratio')

plt.savefig('images/afr.pdf')
plt.show()

m=np.array([0.0001,0.001,0.01, 5.8,6])
for val in m:  
    fnam= 'm' + str(val) + '.txt'    
    z,t,yf,yo2,yco2,yh2o = readDataSheet(fnam)
    plt.plot(z,t,label='$\dot{m_{fuel}} \;x \;' + str(val) +'$')

plt.legend()
plt.xlabel('Mixture fraction - Z')
plt.ylabel('Temperature [K]')
plt.title('Parameter variation - mass flux of fuel')
plt.savefig('images/mdot.pdf')
plt.show()

#plt.plot(z,t,label='$T_b$')
#plt.show()

val=5.0;
fnam= 'm' + str(val) + '.txt'    
z,t,yf,yo2,yco2,yh2o = readDataSheet(fnam)
    
plt.plot(z,yco2,label='$CO_2$')
plt.plot(z,yh2o,label='$H_2O$')
plt.plot(z,yo2,label='$O_2$')
plt.plot(z,yf,label='$CH_4$')
plt.plot([z_st, z_st],[0, 1],'--k',linewidth=0.3)

plt.ylim((0, 0.22))
plt.xlim((0,0.25))
plt.show()


fnam='ref.txt'
z,t,yf,yo2,yco2,yh2o = readDataSheet(fnam)

plt.plot(z,t,label='counterflow')
plt.plot(Zdata,Tdata,'--k',label='B-S',linewidth=0.2)
plt.plot([0, 1],[T_air,T_e],'--k',linewidth=0.2)
plt.xlabel('Mixture fraction - Z')
plt.ylabel('Temperature [K]')
plt.title('Temperature vs Z - Reference case')
plt.legend()
plt.savefig('images/cftref.pdf')
plt.show()

plt.plot(z,yco2,label='$CO_2$')
plt.plot(z,yh2o,label='$H_2O$')
plt.plot(z,yo2,label='$O_2$')
plt.plot(z,yf,label='$CH_4$')
plt.plot(Zdata,YP[:,0],'--k',linewidth=0.1,label='B-S')
plt.plot(Zdata,YP[:,1],'--k',linewidth=0.1)
plt.plot(Zdata,YO2,'--k',linewidth=0.1)
plt.plot(Zdata,YF,'--k',linewidth=0.1)

plt.plot([z_st, z_st],[0, 1],'--k',linewidth=0.5,label='$Z_{st}$')
plt.xlabel('Mixture fraction - Z')
plt.ylabel('Mass fraction')
plt.title('Mass fraction vs Z - Reference case')

plt.legend()

plt.ylim((0, 0.22))
plt.xlim((0,0.25))
plt.savefig('images/cfyref.pdf')
plt.show()

#plt.plot([z_st, z_st],[T_u(z_st), tempbs(z_st)],'--k',label='$Z_{st}$')
#plt.xlim([0, 1])

#
#plt.legend()
#plt.title('Burke-Schumann solution - Temperature')
#plt.xlabel('Mixture fraction - Z')
#plt.ylabel('Temperature [K]')
#plt.savefig('images/b-s-t.pdf')
#plt.show()
#
#
#plt.plot(Zdata,YP[:,0],label='$CO_2$')
#plt.plot(Zdata,YP[:,1],label='$H_2O$')
#plt.plot(Zdata,YO2,label='$O_2$')
#plt.plot(Zdata,YF,label='$CH_4$')
#plt.plot([z_st, z_st],[0, 0.22],'--k',label='$Z_{st}$')
#plt.ylim((0, 0.22))
#plt.xlim((0,1))
#
#plt.title('Burke-Schumann solution - Composition')
#plt.xlabel('Mixture fraction - Z')
#plt.ylabel('Massfraction - Y')
#plt.legend()
#plt.savefig('images/b-s-y1.pdf')
#plt.show()
#
#
#plt.plot(Zdata,YP[:,0],label='$CO_2$')
#plt.plot(Zdata,YP[:,1],label='$H_2O$')
#plt.plot(Zdata,YO2,label='$O_2$')
#plt.plot(Zdata,YF,label='$CH_4$')
#plt.plot([z_st, z_st],[0, 0.22],'--k',label='$Z_{st}$')
#plt.ylim((0, 0.22))
#plt.xlim((0,0.25))
#
#plt.title('Burke-Schumann solution - Composition')
#plt.xlabel('Mixture fraction - Z')
#plt.ylabel('Massfraction - Y')
#plt.legend()
#plt.savefig('images/b-s-y025.pdf')
#plt.show()










