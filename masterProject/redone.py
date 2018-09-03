# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 18:45:19 2018

@author: Gabriel
"""
import numpy as np

import pandas as pd

from openpyxl import load_workbook

from thermopy import nasa9polynomials as nasa9
db=nasa9.Database()

datapath = 'datafiles/'     # path to tables etc
sourcepath = 'sourcefiles/'           # path for resulting sourcefiles to STAR

fnam='Cases_simplified_model.xlsx' # filename for cases to read

# Switch for enabling calculation of just bed with no gasification chamber
justbed=False # for debug purpose

# colnames = ['case_index', 'gas_H2', 'gas_CO', 'gas_CO2', 'gas_CH4', 'gas_C2H2', 'gas_C2H4', 'gas_C2H6', 'gas_C3H6', yield_gas, tar_conc, y_char, y_vol, y_ash, fuel_HHV, fuel_LHV, fuel_C, fuel_H, fuel_O, TG, TC, Tair, Tsteam, Tfuel, P_gas, P_heat, nboil, xH2O_G_wet, xH2O_C_wet, ER, u_air, u_steam, H_bed, L_chamber, W_chamber, 'H_gap', 'wall_thickness', 'D', 'porosity', 'rho_solid']
cases = pd.read_excel(fnam,'Sheet1',header=0,skiprows=[1], usecols='B:AQ') # read cases
cases = cases.fillna(value=cases.iloc[0]) # fill empty values with base case values
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
H = db.set_compound('H')
O = db.set_compound('O')
naphtalene=db.set_compound('c10h8')
naphtalene.elements[0] = ('C',10) # thermopy has parsed the nasa polys wrong for 2 digit subscripts it seems.

# component lists for different gas mixtures
comps_gas = [h2, co, co2, ch4, c2h2, c2h4, c2h6, c3h6]
comps_air = [co2, h2o, o2, n2]
comps_flue = [co2, h2o, o2, n2]


# 0 deg C in Kelvin, for conversion purposes
TK= 273.15 # K
Tref = 25# deg C

# Pressure at one atm
p_atm = 101325 # Pa
dP = 6200 # pressure at bottom of bed
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
CHO_dict={'C': 0, 'H': 1, 'O': 2}


# add lhv to thermopy objects for convenience
for i,comp in enumerate(comps_gas):
    comp.lhv = lhv_gas[i]

# molar masses through thermopy
mw_gas =[]
for comp in comps_gas:
    mw_gas.append(comp.molecular_weight)

# append objects with mass fraction of CHO in each species
comps_reacting = comps_gas + [naphtalene, h2o]
for comp in  comps_reacting:
    CHO = np.array([0.0,0.0,0.0])
    for element in comp.elements:
        weight = MW_CHO[CHO_dict[element[0]]] * element[1]
        CHO[int(CHO_dict[element[0]])] = weight
    comp.CHO = CHO / CHO.sum()
        



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
def cpChar(Tbed,Tref=Tref):
    Tmean = (Tbed+Tref)/2 + TK # K
    Cp_ch = -0.0038*Tmean**2 + 5.98*Tmean - 795.28 # j/kgK
    return Cp_ch/1000 # kj/kgK

# Composition and yield of flue gas based on fuel CHO and equivalence ratio
#  kg element /kgdaf -> mol/(mol gas), (mol flue/kgdaf), (mol air/kgdaf)
def fluegasComp(CHO_fuel,ER):
    X_fuel = np.array(CHO_fuel) / MW_CHO
    co2 = X_fuel[0]
    h2o = X_fuel[1]/2
    o2 = co2 + h2o/2 - X_fuel[2]/2
    X = np.array([co2, h2o, o2*(ER-1), o2*ER*3.76])
    yFg = X.sum()       # mol flue/kgdaf
    
    yAir  =  o2*ER + o2*ER*3.76
    return X/yFg, yFg, yAir 

# Empirical relationship between total and primary air
# total air in mol/s -> (total_air/primary_air)
def primaryAirFraction(total_air):
    tot=np.array([39.22408417, 69.4281724])/3.6 # Nm3/s
    tot = tot*p_atm/(R*(Tref+TK))
    xpr = [0.724031165, 0.4475699422]
    
    k = (xpr[1] - xpr[0])/(tot[1] - tot[0])
    m = xpr[0] - k*tot[0]
    
    Xprimary = k*total_air + m
    return Xprimary


# Cp for silica sand from Ihsan Barin. "Thermochemical data of Pure substances"
#  Temperature in K -> (j/(kg K)    
def cp_solid(T):
    coeff= [1.5660000E-07, -8.3755230E-04, 1.5971295E+00, -4.2836000E+01, -1.5288225E+07];
    cp = coeff[0] * 4*T**3 + coeff[1] * 3*T**2 + coeff[2]*2*T + coeff[3]
    return cp


# Elemental massfraction of mixture
def elemComp(comps, X):
    CHO = np.zeros_like(MW_CHO)
    mixMW=mixtureMW(comps,X)
    for i, comp in enumerate(comps):
        CHO = CHO + X[i]*comp.molecular_weight*comp.CHO / mixMW
    return CHO / CHO.sum()
        

# -------------------------------------------------------------------------
# Definition of case as an object for convenience (less brackets)
# -------------------------------------------------------------------------
class Case(object):
    # Initialisation method to allow passing dict or keyword arguments 
    # to change individual attributes in order to create different cases
    def __init__(self, df=cases,index=0):
        dictionary = df.to_dict('records')[index]
        dictionary['index'] = int(dictionary['case_index'])
        for key in dictionary:
            setattr(self, key, dictionary[key])
            
# -------------------------------------------------------------------------
# Source term calculation (enclose in loop or make function?)
# -------------------------------------------------------------------------

index = 0
case = Case(cases,index) # load case specific indata



def calcSource(case):
    # Collection of given gas composition into numpy array
    
    case.gas_order = ['h2', 'co','co2','ch4','c2h2','c2h4','c2h6','c3h6']
    case.x_gas = np.array([case.gas_H2, case.gas_CO, case.gas_CO2, 
                           case.gas_CH4, case.gas_C2H2, case.gas_C2H4, case.gas_C2H6, case.gas_C3H6])/100
    
    # Collection of given fuel elemental composition into numpy array
    CHO = np.array([case.fuel_C, case.fuel_H, case.fuel_O])/100
    case.CHO = CHO#/CHO.sum()
    
    # make char and volatile yield into daf basis
    pa = np.array([case.y_char, case.y_vol])
    pa = pa/pa.sum()
    
    case.y_char = pa[0]
    case.y_vol = pa[1]
    
    
    # specify extra params not in excel
    case.W_bed = 5.9
    case.L_bed = 7.398
    
    # ------------ Preprocess parameters --------------------------
    
    # Moisture made into dry basis
    case.xh2o_g = case.xH2O_G_wet / (1 - case.xH2O_G_wet)
    case.xh2o_c = case.xH2O_C_wet / (1 - case.xH2O_C_wet)
        
    
    # ---------- Calculate massflows in gasifier ------------------------------
    
    lhv_cg = mixtureLHV(comps_gas,case.x_gas)   #  MJ/kg, lhv of produced gas
    
    m_cg = case.P_gas / lhv_cg                  # kg/s, needed massflow cold gas    
    
    m_gfuel = m_cg / case.yield_gas         # kg/s, needed massflow gasifier fuel
    
    # Massflow of steam
    th = case.wall_thickness
    A_wall = th*case.L_chamber + 2*(th*case.W_chamber + th*th)
    A_G = case.L_chamber * case.W_chamber # Crossectional area of gasifier
    
    V_steam = case.u_steam * (A_G + A_wall)
    
    #interpolate in table for density at temp = Tbed
    rho_steam = np.interp(case.TG + TK,steam_props.Temp,steam_props.rho) # kg/m3
    
    # mass flow of steam  
    m_steam = V_steam * rho_steam
    
    SB = m_steam / m_gfuel # steam to biomass ratio, kg steam / kg daf ### ADD steam to fuel ratio SFR
    
    SFR = (m_steam+ m_gfuel*case.xh2o_g) / m_gfuel
    
    # Calculate tar yield and degree of char conversion
    
    n_gas = m_cg / mixtureMW(comps_gas,case.x_gas) # mol/s, cold gas
    
    V_gas = n_gas*R*(Tref+TK)/p_atm # Nm3/s, cold gas
    
    m_tar = case.tar_conc * V_gas / 1000 # kg/s, tar
    
    # Carbon balance
    gas_CHO=elemComp(comps_gas,case.x_gas)
    
    char_C = case.CHO[0] - case.yield_gas*gas_CHO[0] - m_tar/m_gfuel*naphtalene.CHO[0]
    
    steam_O = m_cg/m_gfuel*gas_CHO[2] + m_tar/m_gfuel*naphtalene.CHO[2] - case.CHO[2] 
    
    X_ch = (case.y_char - char_C)/case.y_char 
    X_steam = (steam_O*m_gfuel)/(m_steam*h2o.CHO[2])
    
    # ---------- calculate energy demands in gasifier -------
    
    # Sensible heat of of char kj /kg daf
    hchar = case.y_char * cpChar(case.TG) * (case.TG - Tref) 
    
    # Sensible heat of volatiles kj/kg daf
    hvol = case.y_vol * meanCpmix(comps_gas,case.x_gas,Tref,case.TG) * (case.TG - Tref)   
    
    # Sensible heat of entering fuel kj/ kgdaf
    hfuel = cpFuel(case.Tfuel) * (case.Tfuel - Tref)
    
    
    cp_w_uni = case.y_char*cpChar(case.T_whole) + case.y_vol* meanCpmix(comps_gas,case.x_gas,Tref,case.T_whole)
    cp_w_TG =  case.y_char*cpChar(case.TG,case.T_whole) + case.y_vol* meanCpmix(comps_gas,case.x_gas,case.T_whole,case.TG)
    
    K_fuel = m_gfuel*cp_w_TG*case.Xbed
    const_fuel = m_gfuel*(cp_w_uni*(case.T_whole-Tref)-hfuel)- K_fuel*(case.T_whole+TK)
    
    K_moist= m_gfuel*case.xh2o_g*case.Xbed*meanCp(h2o,case.TG,case.T_whole)/h2o.molecular_weight
    const_moist = m_gfuel*case.xh2o_g*(enthalpy(h2o,case.T_whole) - enthalpy(h2ol,Tref)) - K_moist*(case.T_whole+TK)
    
    const_g = (const_fuel + const_moist)/1e3 
    K_g = (K_fuel+K_moist)/1e3
    
    sink_gnew=const_g+ K_g*(case.TG+TK)
    
    # energy for heating fuel kj / kgdaf
    qfuel_g = hchar + hvol - hfuel
    
    # energy for heating and drying water kj/kgdaf
    qh2o_g = case.xh2o_g*(enthalpy(h2o,case.TG) - enthalpy(h2ol,case.Tfuel))
    
    # energy for heating fluidising steam kj / kgdaf
    qsteam_g = SB * (enthalpy(h2o,case.TG) - enthalpy(h2o,case.Tsteam))
    
    #  energy needed for char conversion reaction kj/kgdaf
    qchar = case.y_char*X_ch*dHr_g
    
    # Total heat demand in gasifier
    Q_g = sink_gnew + m_gfuel*( + qsteam_g + qchar)/1000 # MW
    
    # --------------- mass flows in combustion chamber ------------ 
    
    m_bonuschar=(1-X_ch)* case.y_char*m_gfuel
    m_cfuel = (case.P_fuel + Q_g - lhv_ch*m_bonuschar) / case.fuel_LHV # needed massflow of fuel to combustion
    
    
    # For generating case without gasification chamber
    if justbed:
        m_bonuschar = 0
        m_cfuel = case.P_heat/case.nboil / case.fuel_LHV 
        
    CHO_mix = m_bonuschar* np.array([1,0,0]) + m_cfuel*case.CHO
    CHO_mix = CHO_mix / CHO_mix.sum()
    
    # Calculate outgoing fluegas composition and total air demand for this case
    x_flue, yFlue, yAir = fluegasComp(CHO_mix,case.ER)   # mol fraction, mol/kg fuel
    
    # Molar flows of air and flue gas
    
    A_C = case.L_bed * case.W_bed - A_G - A_wall
    
    # For generating case without gasification chamber
    if justbed:
        A_C = case.L_bed * case.W_bed
    
    V_fm = case.u_air * A_C 
    
    
    # Molar flow of fluidising media from ideal gas law
    n_fm_tot = V_fm * (p_atm + dP) / (8.3144 * (case.TC+TK))
    
    # molar flow of air needed for this amount of fuel
    n_air_tot = (m_cfuel + m_bonuschar)*yAir # mol/s
    
    #  amount of primary air through empirical relationship
    n_primary = primaryAirFraction(n_air_tot)*n_air_tot
    
    # amount of flue gas recirculated
    n_flue = n_fm_tot - n_primary
    
    # component amounts in fluidising gas
    n_fm = n_primary*x_air + n_flue*x_flue
    
    x_fm = n_fm/n_fm_tot # %-mol composition 
    
    # fraction of total fluidising media that is fluegas
    Xflue = n_flue/n_fm_tot
    
    
    # massflow if fluidising media
    m_fm = n_fm_tot*mixtureMW(comps_flue,x_fm)
    
    # ------------------ Energy demands in combustion chamber ----------------- 
    
        
    # Sensible heat of of char kj /kg daf
    hchar = case.y_char * cpChar(case.TC) * (case.TC - Tref) 
    
    # Sensible heat of volatiles kj/kg daf
    hvol = case.y_vol * meanCpmix(comps_gas,case.x_gas,Tref,case.TC) * (case.TC - Tref)   
    
    # Sensible heat of entering fuel kj/ kgdaf
    hfuel = cpFuel(case.Tfuel) * (case.Tfuel - Tref)
    
    # energy for heating fuel kj / kgdaf
    qfuel_c = hchar + hvol - hfuel
    
    # energy for heating and drying water kj/kgdaf
    qh2o_c = case.xh2o_c*(enthalpy(h2o,case.TC) - enthalpy(h2ol,case.Tfuel))
    
    # Energy for heating fluidising media (kj/mol fluidising media)
    qfm = meanCpmix(comps_flue,x_fm,case.Tair,case.TC,True) * (case.TC - case.Tair)
    
    Q_c = (m_cfuel*(qfuel_c + qh2o_c) + n_fm_tot*qfm)/1000 
    
    # ------------------ Energy balance tally ----------------
    Qtot = Q_c + Q_g
    
    # For generating case without gasification chamber
    if justbed:
        Qtot = Q_c
    
    # heat released from combustion
    Q_comb = n_fm[2]*dHr_c
    
    # mismatch between heat released and heat needed.
    Q_corr = Qtot - Q_comb
    
    # ------------------------------------------------------------------------- 
    # Postprocessing
    # -------------------------------------------------------------------------
    # Formulation of the source terms
    S_c = Qtot
    
    sink_c = m_cfuel*(qfuel_c + qh2o_c) /1e3  # MW, heat needed to dry combusting fuel
    
    S_g = m_gfuel*qchar/1e3            # MW, heat needed for char gasification
    
    sink_g = m_gfuel*(qh2o_g + qfuel_g)/1e3                 # MW, heat needed to dry gasifying fuel
    
    # For generating case without gasification chamber
    if justbed:
        S_g=0
        sink_g=0
        case.L_chamber=0
        case.W_chamber=0
        case.wall_thickness =0
    
    S_tot = S_c - S_g - sink_g - sink_c             # summation of source terms to correspond to a volume integral over domain in star.
    
    # Calculation of k_eff from case data
    k_eff = case.D*case.rho_solid*cp_solid((case.TC + case.TG)/2+TK)*(1-case.porosity)
    
    # Fluid properties at 800 - 850 deg C
    mu_air = 44.0e-6 # Pa s 
    
    tc_air = 68.0e-3 # W/(m K)
    
    cp_air =  meanCpmix(comps_flue,x_fm,case.Tair,case.TC,True) / mixtureMW(comps_flue,x_fm) * 1e3 # J/(kg K)
    
    mu_steam = np.interp(case.TG + TK,steam_props.Temp,steam_props.mu) * 1e-6 # Pa s
    
    tc_steam = np.interp(case.TG + TK,steam_props.Temp,steam_props.tc) * 1e-3 # W/(m K) 
    
    cp_steam = np.interp(case.TG + TK,steam_props.Temp,steam_props.Cp) * 1000 # J/(kg K) 
    
    
    # Formulation of heating source term based on the other source terms.
    
    cp_air = meanCpmix(comps_flue,x_fm,case.Tair,case.TC) * 1e3
    S_c_steam = m_steam * cp_steam * (case.TG - case.Tsteam) /1e6
    
    S_c = m_fm * cp_air * (case.TC - case.Tair) /1e6 + sink_c + S_c_steam + S_g
    
    #############################
    #################################
    
    #--------------------------------------------------------------------------
    # Write sourcefile for STAR CCM+ to read
    
    
    data = [1, case.case_index, S_c, -sink_c, -S_g, const_g, K_g, m_fm, m_steam, k_eff,
            mu_air, mu_steam,tc_air,tc_steam, cp_air, cp_steam, case.rho_solid,
            case.TG+TK,case.TC+TK,case.Tair+TK,case.Tsteam+TK,
            case.porosity, case.L_chamber, case.W_chamber, case.H_gap,case.wall_thickness, case.L_bed,
            case.W_bed, case.H_bed]
    header = 'row,index,S_c,sink_c,S_g,const_g,K_g,m_fm,m_steam,k_eff,mu_air,mu_steam,tc_air,tc_steam,cp_air,cp_steam,rho_solid,TG,TC,Tair,Tsteam,porosity,L_chamber,W_chamber,H_gap,chamber_thickness,L,W,H'
    
    exceldata = [case.case_index, m_cfuel, qh2o_c*m_cfuel/1e3, qfuel_c*m_cfuel/1e3, m_fm * cp_air * (case.TC - case.Tair) /1e6, Q_c, Q_comb, Q_corr, Xflue, m_gfuel, qh2o_g*m_gfuel/1e3, qfuel_g*m_gfuel/1e3, m_steam * cp_steam * (case.TG - case.Tsteam) /1e6, qchar*m_gfuel/1e3, Q_g, m_cfuel, m_cfuel*case.fuel_LHV, SB, SFR, X_ch, X_steam, S_c, sink_c, S_g, sink_g,S_tot ]
    excelheader = 'index,m_fuel_c,qh2o_c,q_fuel_c,q_fm_c,Q_c,Qcomb,Qcorr,Xflue,m_fuel_g,q_h2o_g,q_fuel_g,q_steam_g,q_charconv,Q_g,load_kg,load_MW,SB,SFR,X_ch,X_steam,S_c,sink_c,S_g,sink_g,S_tot'
    #excelunits = 


    return data, header, exceldata, excelheader
    


sourcepath='/scratch/gabgus/geometry/redone/sourcefiles/' 


# Function to write calculated source terms to STAR input file
def writeSource(case,data,header,casename,sourcepath=sourcepath):
    with open(sourcepath + casename, 'w') as f:
        f.write(header)
        f.write('\n')
        for i,val in enumerate(data):    
            f.write(str(val) + ',')
        f.close()


excellist = []
## loop to generate all cases for this batch
for index in range(len(cases)):
    case = Case(cases,index)
        
    casestring = str(index)
    
    data, header, exceldata, excelheader = calcSource(case)
    
    casename = 'source_' + casestring + '.csv'
    writeSource(case,data,header,casename)
    
    exceldata[0] = casestring
    excellist.append(exceldata)





### Test matrix source generation

## Test matrix formulation
#param_order = 'P_gas, D, Xbed, xH2O_G_wet, Tsteam, L_chamber, W_chamber, H_gap'
#paramrange = [ [12.0, 8.0], [0.03,0.01], [1.0,0.4], [0.55,0.2], [450.0,200.0], [3.0,2.2], [1.2,0.8], [0.4,0.3] ]
#n = len(paramrange)    
#N = 2**n
#
#ic=[]
#indc=[]
#for i in range(n):
#    b = [0]*2**i + [1]*2**i
#    paramrange[i].sort()
#    a= np.choose(b,paramrange[i])
#    ic.append(np.resize(a,N))
#    indc.append(np.resize(b,N))
#    
#    
#    
#ic.reverse()
#indc.reverse()
#param_order = param_order.split(', ')
#param_order.reverse()
#
#matrix = np.column_stack(tuple(ic))
#imatrix = np.column_stack(tuple(indc))
#
#
#
#
#
#excellist = []
#
#for i in range(N):
#    case = Case(cases,0)
#    case.P_gas = ic[7][i]
#    case.D = ic[6][i]
#    case.Xbed = ic[5][i]
#    case.xH2O_G_wet = ic[4][i]
#    case.Tsteam = ic[3][i]
#    case.L_chamber=ic[2][i]
#    case.W_chamber = ic[1][i]
#    case.H_gap = ic[0][i]
#    
#    casestring=''
#    for k in range(len(imatrix[i])):
#        casestring = casestring + str(imatrix[i][k])
#    
#    case.case_index=int(casestring)
#    
#    
#    data, header, exceldata, excelheader = calcSource(case)
#    
#    casename = 'source_' + casestring + '.csv'
#    #writeSource(case,data,header,casename)
#    
#    exceldata[0] = casestring
#    excellist.append(exceldata)


exceldf = pd.DataFrame(excellist,columns=excelheader.split(','))

book = load_workbook(fnam)

writer = pd.ExcelWriter(fnam, engine='openpyxl') 
writer.book = book
writer.sheets = dict((ws.title, ws) for ws in book.worksheets)

exceldf.to_excel(writer, "Sheet2", startrow=1)

writer.save()
#exceldf.to_excel(fnam,'Sheet2',startrow=2)


#
#P_gas = np.resize(a[0],N)
#D = np.resize(a[1],N)
#X_vola_mois = np.resize(a[2],N)
#xH2O_G_wet = np.resize(a[3],N)
#Tsteam = np.resize(a[4],N)
#L_chamber = np.resize(a[5],N)
#W_chamber = np.resize([0,1],N)
#H_gap = np.resize([0,1],N)


#Fukt = np.choose(Fukt,[0.08,0.7])
#char = np.choose(char,[0,1])
#D = np.choose(D,[20,100])

#matrix = np.column_stack((Fukt,char,D))




# scrap calculation to get loads for empiric data


CHONS = np.array([0.035, 0.518, 0.06, 0.381, 0.0054, 0.0009])
CHONS_daf = CHONS[1:]/CHONS[1:].sum()

lhv_fuel = 20.456 # MW/kgdaf

MW_CHONS=np.array([0.01201 , 0.001007, 0.015998, 0.0140067, 0.032])


X_CHONS = CHONS_daf/MW_CHONS # mol element/kgdaf 

n_co2 = X_CHONS[0] # mol co2 /kgdaf

n_h2o = X_CHONS[1]/2 # mol h2o/kgdaf

n_so2 = X_CHONS[4]

n_o2 = (X_CHONS[0] + X_CHONS[1]/4 - X_CHONS[2]/2 + X_CHONS[4]) # mol o2 /kgdaf

n_n2 = n_o2*3.76 + X_CHONS[3]/2    # mol n2/kgdaf

lt = (n_o2 + n_n2)* R*(TK)/p_atm # theoretical Nm3 air/kgdaf


#SEEMS WRONG? need o2conc in dry fluegas -> need total amount of flue gas
o2conc=np.array([0.058, 0.048]) # vol % o2 in moist flue gas?


xh2o_wet = 0.55
xh2o_daf = xh2o_wet /(1- xh2o_wet) 

n_moisture = xh2o_daf/ h2o.molecular_weight

gt = (n_co2 + n_n2 + n_so2)* R*(TK)/p_atm # theoretical Nm3 gas / kgdaf

# need to remake with dry o2conc
m = 1.0 + gt/lt * o2conc/(0.21 - o2conc) # air factor 

lv = lt*m  # real Nm3 air/kgdaf 

ltot = np.array([39.2241, 69.4282])/3.6 # Nm3 total air 

mfuel= ltot/lv # kgdaf 

Qfuel = mfuel * lhv_fuel


print(m)
print(Qfuel)










    
    
    