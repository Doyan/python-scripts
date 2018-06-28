# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 18:45:19 2018

@author: Gabriel
"""
import numpy as np

import pandas as pd

from thermopy import nasa9polynomials as nasa9
db=nasa9.Database()

datapath = 'datafiles/'
sourcepath = './'


# -------------------------------------------------------------------------
# Table values
# -------------------------------------------------------------------------

# thermodynamic data for compounds through thermopy:
h2o = db.set_compound('h2o')
h2ol = db.set_compound('h2o(l)')
h2 = db.set_compound('h2') 
co = db.set_compound('carbon monoxide')
co2 = db.set_compound('carbon dioxide')
ch4 = db.set_compound('ch4')
c2h2 = db.set_compound('C2H2,acetylene')
c2h4 = db.set_compound('c2h4')
c2h6 = db.set_compound('c2h6')
c3h6 = db.set_compound('c3h6,propylene')
n2 = db.set_compound('n2')
o2 = db.set_compound('o2')
C = db.set_compound('c(gr)')

# component lists for different gas mixtures
comps_gas = [h2, co, co2, ch4, c2h2, c2h4, c2h6, c3h6]
comps_air = [co2, h2o, o2, n2]
comps_flue = [co2, h2o, o2, n2]

# 0 deg C in Kelvin, for conversion purposes
TK= 273.15 # K
Tref = 25 # deg C

# Pressure at one atm
p_atm = 101325 # Pa

# molar Gas constant 
R=8.3144

# molar composition of air 
x_air = np.array([0, 0, 1, 3.76])/4.76

# enthalpy of vaporisation for water.
dh_vap = h2o.enthalpy(373.15) - h2ol.enthalpy(373.15) # J/mol

# Heating value of anthracite from "Data och Diagram"
lhv_ch = 32.7                 # MJ/kg to represent char

# molar heat of combustion for carbon (i.e heat of formation for CO2)
dHr_c = 393.5/1000 # MJ/mol 

# Reaction enthalpy for char conversion
dHr_g = 131.48 # kJ/mol  reaction: C(s) + H2O -> CO + H2

# heating values for gas components
# from NREL answer to query:
# http://www.ieatask33.org/app/webroot/files/file/publications/HeatingValue.pdf

names_gas = ['h2', 'co', 'co2', 'ch4', 'c2h2', 'c2h4', 'c2h6', 'c3h6']
lhv_gas = [241.79, 282.9, 0, 802.71, 1256.9, 1323.2, 1428.83, 1926.1]     # kj/mol


# https://doi.org/10.1063/1.1285884
# units: K, mol/dm3, J/mol, J/mol, J/(mol K), J/(mol K), J/(mol K), m/s 
air_props = pd.read_csv(datapath + 'air_1atm_molar.csv',index_col=False)


# from handbook of chemistry and physics at 1bar
# Units: K, kg/m3, kj/kg, kj/(kgK),kj/(kg K),kj/(kg K),m/s, microPa s,mW/(m K) 
steam_props = pd.read_csv(datapath + 'steam_1bar.csv', index_col=False)

# molar weight of elements
MW_CHO = np.array([12.01, 1.007, 15.998])/1e3 # kg/mol

# add lhv to thermopy objects for convenience
for i,comp in enumerate(comps_gas):
    comp.lhv = lhv_gas[i]

# molar masses through thermopy
mw_gas =[]
for comp in comps_gas:
    mw_gas.append(comp.molecular_weight)

# -------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------

# Molar mass of mixture
# (thermopy compound list, molar fraction list,) -> kg/mol    

def mixtureMW(comps,X):
    MWmid = []
    X= X / sum(X)
    for i, comp in enumerate(comps):
        MWmid.append(X[i] * comp.molecular_weight)
    return sum(MWmid)


# Lower heating value of mixture
# (thermopy compound list, molar fraction list,) -> MJ/kg    

def mixtureLHV(comps,X):
    lhv = []
    X= X / sum(X)
    for i, comp in enumerate(comps):
        lhv.append(X[i] * comp.lhv )
    return sum(lhv) / mixtureMW(comps,X)/1000

# Specific enthalpy 
# (thermopy compound, T in Celsius) -> kj / kg

def enthalpy(comp,T): 
    return comp.enthalpy(T+TK)/comp.molecular_weight/1e3

# Specific enthalpy of mixture 
# (thermopy compound list, molar fraction list, T in Celsius) -> kJ/kg

def enthalpyMix(comps,X,T):
    h = []
    X= X / sum(X)
    for i, comp in enumerate(comps):
        h.append(X[i] * comp.enthalpy(T+TK)/1000)
    return sum(h)/ mixtureMW(comps,X)



# Mean specific heat capacity of thermopy compund over interval 
# (thermopy compound, start temp in deg C, end temp in deg C) -> kj/mol    

def meanCp(comp,T1,T2):
    cp = []
    n = 20
    Tspan = np.linspace(min([T1,T2])+TK,max([T2,T1])+TK,n)
    
    for temp in Tspan:
        cp.append(comp.heat_capacity(temp))
    return sum(cp)/n/1e3


# Mean specific heat capacity of mixture over interval
# (thermopy compound list, molar fraction list, ...
#  ... start temp in deg C, end temp in deg C) -> kj/kgK

def meanCpmix(comps, X, T1, T2, mol=False):
    cp = []
    X = X/X.sum()
    for i, comp in enumerate(comps):
        cp.append(X[i] * meanCp(comp,T1, T2))
    if mol:
        return sum(cp)
    else:
        return sum(cp) / mixtureMW(comps,X)    

# Mean Cp for fuel over interval [Tref - Tfuel] using correllation  
# by  Gupta et al 2003, https://doi.org/10.1016/S0016-2361(02)00398-8
# (fuel temperature in celsius) -> kj/kgK  
def cpFuel(Tfuel):
    Tmean = (Tfuel+Tref+2*TK)/2 # K 
    Cp_f = (5.46*(Tmean)/2 - 524.77)/1000 # kj/kgK 
    return Cp_f

# Mean Cp for char over interval [Tref - Tbed] 
# correllation also from Gupta et al 2003 
# (bed temperature in celsius) -> kj/kgK
def cpChar(Tbed):
    Tmean = (Tbed+Tref+2*TK)/2 # K
    Cp_ch = -0.0038*Tmean**2 + 5.98*Tmean - 795.28 # j/kgK
    return Cp_ch/1000 # kj/kgK

# Composition and yield of flue gas based on fuel CHO and equivalence ratio
#  kg element /kgdaf -> mol/(mol gas), (mol gas/mol fuel)
def fluegasComp(CHO_fuel,ER):
    X_fuel = np.array(CHO_fuel) / MW_CHO
    co2 = X_fuel[0]
    h2o = X_fuel[1]/2
    o2 = co2 + h2o/2 - X_fuel[2]/2
    X = np.array([co2, h2o, o2*(ER-1), o2*ER*3.76])
    Xair = X[2:]
    yAir = Xair.sum()
    yFg = X.sum()
    return X/yFg, yFg, yAir 

# Empirical relationship between total and primary air
# total air in mol/s
def primaryAirFraction(total_air):
    tot=np.array([39.22408417, 69.4281724])/3.6 # Nm3/s
    tot = tot*p_atm/(R*(Tref+TK))
    xpr = [0.724031165, 0.4475699422]
    
    k = (xpr[1] - xpr[0])/(tot[1] - tot[0])
    m = xpr[0] - k*tot[0]
    
    Xprimary = k*total_air + m
    return Xprimary

# -------------------------------------------------------------------------
# Definition of base case
# -------------------------------------------------------------------------
class Case(object):
    """Base Case class to start from"""
    
    # Collection of given gas compositions in mol % 
    gas_order = ['h2', 'co','co2','ch4','c2h2','c2h4','c2h6','c3h6']
    gas_comps = {'Anton_high':[40.58, 16.03, 31.54, 8.31, 0.38, 2.88, 0.13, 0.16], 
              'Anton_mid': [36.32, 20.48, 28.65, 9.94, 0.5, 3.62, 0.21, 0.27],
              'Anton_low': [27.51, 25.92, 28.25, 12.16, 0.73, 4.45, 0.26, 0.73]
    }
    
    # Corresponding gas yields in kg/kgdaf      
    gas_yields = {'Anton_high': 0.88, 'Anton_mid': 0.7, 'Anton_low': 0.53}
    
    # ---------------------------------------------------------------------
    # Inputs 
    gascombo='Anton_mid'    # namestring for chosen gas comp and yield
    
    X_gas = np.array(gas_comps[gascombo])/100    # molar fraction, gas composition
    yield_gas = gas_yields[gascombo]                # kg gas/kgdaf
    
    # Fuel properties, change all at once!
    y_char = 0.2622875      # yield char kg/kgdaf
    y_vol = 0.703075        # yield volatiles kg/kgdaf	
    y_ash = 0.0391125       # yield ash kg/kgdaf 
    
    fuel_HHV = 20.50125     # MJ/kg
    fuel_LHV = 19.25125     # MJ/kg
    
    CHO = np.array([52.675, 5.7375, 37.625]) / 100 # mass %, major elements     
    
    # Temperatures
    TG = 800        # deg C
    TC = 850        # deg C
    Tair = 200      # deg C
    Tsteam = 200    # deg C
    Tfuel = 40      # deg C
    
    # Loads
    P_gas = 10      # MW
    P_heat = 85     # MW
    
    # Boiler efficiency (stack losses)
    nboil = 0.89
   
     # Moisture content, wet basis
    #xh2o_g_wet = 0.55
    #xh2o_c_wet = 0.55
    
    # Equivalence ratio
    ER = 1.2
    
    # Velocities
    u_air = 1.2     # m/s
    u_steam = 1.2   # m/s
    
    # Dimensions
    L = 7.398   # m, Length of furnace
    W = 5.9     # m, Width of furnace
    H = 1       # m, Bed height
    
    L_chamber = 2.6             # m, Length of chamber
    W_chamber = 1               # m, Width of chamber
    H_gap = 0.5                 # m, Height of gap
    wall_thickness = 0.15    # m, Thickness of chamber wall
   
    # degree of char conversion
    Xch = 0
    
    # fraction of combustion chambers heat demand satisfied from bed
    frombed = 0.45
    
    # Values used to calculate k_eff
    D = 0.02            # m2/s, Particle dispersion
    cp_solid = 1000     # J/kgK, Heat capacity of silica sand
    porosity = 0.525    # void fraction in bed
    rho_solid = 2800    # kg/m3 bed marerial density 
    
    #----------------------------------------------------------------------
    

    # Initialisation method to allow passing dict or keyword arguments 
    # to change individual attributes in order to create different cases
    def __init__(self, fnam='Cases_simplified_model.xlsx',index=0):
        df = pd.read_excel(fnam,'Sheet1',header=0,skiprows=[1])
        df = df.fillna(value=df.iloc[0])
        dictionary = df.to_dict('records')[index]
        dictionary['index'] = int(dictionary['case_index'])
        for key in dictionary:
            setattr(self, key, dictionary[key])
# -------------------------------------------------------------------------
# Source term calculation (enclose in loop or make function?)
# -------------------------------------------------------------------------
case = Case() # load case specific indata


# ------------ Preprocess parameters --------------------------

# Moisture made into dry basis
case.xh2o_g = case.xH2O_G_wet / (1 - case.xH2O_G_wet)
case.xh2o_c = case.xH2O_C_wet / (1 - case.xH2O_C_wet)
    
# Calculate outgoing fluegas composition and total air demand for this case
x_flue, yFlue, yAir = fluegasComp(case.CHO,case.ER)   # mol fraction, mol/kg fuel


# ---------- Calculate massflows in gasifier ------------------------------

lhv_cg = mixtureLHV(comps_gas,case.X_gas)   #  MJ/kg, lhv of produced gas

m_cg = case.P_gas / lhv_cg                  # kg/s, needed massflow cold gas    

m_gfuel = m_cg / case.yield_gas         # kg/s, needed massflow gasifier fuel

# Massflow of steam
A_G = case.L_chamber * case.W_chamber # Crossectional area of gasifier

V_steam = case.u_steam * A_G

#interpolate in table for density at temp = Tbed
rho_steam = np.interp(case.TG + TK,steam_props.Temp,steam_props.rho) # kg/m3

# mass flow of steam  
m_steam = V_steam * rho_steam

SB = m_steam / m_gfuel # steam to biomass ratio, kg steam / kg daf

# ---------- calculate energy demands in gasifier -------

# Sensible heat of of char kj /kg daf
hchar = case.y_char * cpChar(case.TG) * (case.TG - Tref) 

# Sensible heat of volatiles kj/kg daf
hvol = case.y_vol * meanCpmix(comps_gas,case.X_gas,Tref,case.TG) * (case.TG - Tref)   

# Sensible heat of entering fuel kj/ kgdaf
hfuel = cpFuel(case.Tfuel) * (case.Tfuel - Tref)

# energy for heating fuel kj / kgdaf
qfuel_g = hchar + hvol - hfuel

# energy for heating and drying water kj/kgdaf
qh2o_g = case.xh2o_g*(enthalpy(h2o,case.TG) - enthalpy(h2ol,case.Tfuel))

# energy for heating fluidising steam kj / kgdaf
qsteam_g = SB * (enthalpy(h2o,case.TG) - enthalpy(h2o,case.Tsteam))

# Total heat demand in gasifier
Q_g = m_gfuel*(qh2o_g + qfuel_g + qsteam_g)/1000 # MW

# --------------- mass flows in combustion chamber ------------ 

m_cfuel = (case.P_heat/case.nboil + Q_g - lhv_ch*(1-case.Xch)* case.y_char*m_gfuel) / case.fuel_LHV # needed massflow of fuel to combustion


# Massflow of air and flue gas
th = case.wall_thickness
A_wall = th*case.L_chamber + 2*(th*case.W_chamber + th*th)
A_C = case.L * case.W - A_G - A_wall

V_fm = case.u_air * A_C 


# Molar flow of fluidising media from ideal gas law
n_fm_tot = V_fm * p_atm / (8.3144 * (case.TC+TK))

# molar flow of air needed for this amount of fuel
n_air = m_cfuel*yAir # mol/s

n_primary = primaryAirFraction(n_air)

n_flue = n_fm_tot - n_primary

n_fm = n_primary*x_air + n_flue*x_flue

x_fm = n_fm/n_fm_tot


#  Initial guess for fraction of fluidising media that is flue gas
Xflue = 0.2

# Iteration to find correct amount of flue gas recirculated
diff = 1
i =0
while diff > 1e-6:
    # molar components in fluidising media
    n_fm = n_fm_tot*Xflue*x_flue + n_fm_tot*(1-Xflue)*x_air
    x_fm = n_fm / n_fm_tot
    
    # -------------------- calculate energy demands in combustor ------
    
    # Sensible heat of of char kj /kg daf
    hchar = case.y_char * cpChar(case.TC) * (case.TC - Tref) 
    
    # Sensible heat of volatiles kj/kg daf
    hvol = case.y_vol * meanCpmix(comps_gas,case.X_gas,Tref,case.TC) * (case.TC - Tref)   
    
    # Sensible heat of entering fuel kj/ kgdaf
    hfuel = cpFuel(case.Tfuel) * (case.Tfuel - Tref)
    
    # energy for heating fuel kj / kgdaf
    qfuel_c = hchar + hvol - hfuel
    
    # energy for heating and drying water kj/kgdaf
    qh2o_c = case.xh2o_c*(enthalpy(h2o,case.TC) - enthalpy(h2ol,case.Tfuel))
    
    # Energy for heating fluidising media (kj/mol fluidising media)
    qfm = meanCpmix(comps_flue,x_fm,case.Tair,case.TC,True) * (case.TC - case.Tair)
    
    Q_c = (m_cfuel*(qfuel_c + qh2o_c)*case.frombed + n_fm_tot*qfm)/1000 
    
    # ------------------ Energy balance tally ----------------
    Qtot = Q_c + Q_g
    
    n_o2 = Qtot/dHr_c
    
    # New value for fraction that is flue gas
    Xflue_new = 1 - n_o2/(n_fm_tot*x_air[2])
    
    # Calculate tolerance, advance counter and set iterated variable 
    diff = np.abs(Xflue-Xflue_new)
    i += 1
    Xflue = Xflue_new

# ------------------------------------------------------------------------- 
# Postprocessing
# -------------------------------------------------------------------------
# Formulation of the source terms

S_c = Qtot - 5  # MW, heat released from combustion

sink_c = m_cfuel*(qfuel_c + qh2o_c)*case.frombed /1e3   # MW, heat needed to dry combusting fuel

S_g = m_gfuel*case.y_char*case.Xch*dHr_g/1e3            # MW, heat needed for char gasification

sink_g = m_gfuel*(qh2o_g + qfuel_g)/1e3                 # MW, heat needed to dry gasifying fuel

S_tot = S_c - S_g - sink_g - sink_c

# Calculation of k_eff from case data
k_eff = case.D*case.rho_solid*case.cp_solid*case.porosity

# Fluid properties at 800 - 850 deg C
mu_air = 44.0e-6 # Pa s

tc_air = 68.0e-3 # W/(m K)

cp_air =  1160 # J/(kg K)

mu_steam = np.interp(case.TG + TK,steam_props.Temp,steam_props.mu) * 1e-6 # Pa s

tc_steam = np.interp(case.TG + TK,steam_props.Temp,steam_props.tc) * 1e-3 # W/(m K) 

cp_steam = np.interp(case.TG + TK,steam_props.Temp,steam_props.Cp) * 1000 # J/(kg K) 



# Preprocess velocities to match T = 200 C, 
# since simulation uses ideal gas law

V_fm_in= n_fm_tot / np.interp(case.Tair + TK,air_props.Temp,air_props.rho) / 1000
v_air = V_fm_in / A_C # m/s

V_steam_in =  m_steam / np.interp(case.Tsteam + TK,steam_props.Temp,steam_props.rho)
v_steam = V_steam_in / A_G # m/s

#--------------------------------------------------------------------------
# Write sourcefile for STAR CCM+ to read

data = [1, case.index, S_c, -sink_c, -S_g, -sink_g, v_air, v_steam, k_eff,
        mu_air, mu_steam,tc_air,tc_steam, cp_air, cp_steam, case.cp_solid, case.rho_solid,
        case.TG+TK,case.TC+TK,case.Tair+TK,case.Tsteam+TK,
        case.porosity, case.L_chamber, case.W_chamber, case.H_gap,case.wall_thickness, case.L,
        case.W, case.H]
header = 'row,index,S_c,sink_c,S_g,sink_g,v_air,v_steam,k_eff,mu_air,mu_steam,tc_air,tc_steam,cp_air,cp_steam,cp_solid,rho_solid,TG,TC,Tair,Tsteam,porosity,L_chamber,W_chamber,H_gap,chamber_thickness,L,W,H'

with open(sourcepath + 'source_' + str(case.index) + '.csv', 'w') as f:
    f.write(header)
    f.write('\n')
    for i,val in enumerate(data):    
        f.write(str(val) + ',')
    f.close()

# S_c, S_g sink_c sink_g v_steam v_air k_eff mu_air mu_steam cp_air cp_steam

# cp_solid, rho_solid, porosity

# L_chamber, W_chamber, H_chamber, thickness, L,W,H,

# case_index

# some kind of appending case_name -- index dict


 













    
    
    