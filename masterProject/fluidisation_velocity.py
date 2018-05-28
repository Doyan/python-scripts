#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 27 14:08:18 2018

@author: gabgus
"""


import numpy as np, matplotlib.pyplot as plt
import pandas as pd


air = pd.read_csv('air_1bar.csv',index_col=False)
steam = pd.read_csv('steam_1bar.csv', index_col=False)

TK = 273.15
T_air = 850 + TK
T_steam = 850 + TK


mu_air=np.interp(T_air,air['T'],air.mu)*1e-6
rho_air = np.interp(T_air,air['T'],air.rho)

mu_steam = np.interp(T_steam,steam['T'],steam.mu)*1e-6
rho_steam = np.interp(T_steam,steam['T'],steam.rho)


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

phi=0.9 # sphericity
shift=1.34

epsmf=(14*phi)**(-1/3) # voidage at minimum fluidisation


def Umf(rho,mu):
    Ar=rho*(rhop-rho)*g*dp**3/(mu**2)

    Remf=np.sqrt(27.2**2 + 0.0408*Ar) - 27.2
    Rec=1.24*Ar**0.45

    U_mf=Remf*mu/(rho*dp)
    U_c=Rec*mu/(rho*dp)
    return U_mf, U_c 

