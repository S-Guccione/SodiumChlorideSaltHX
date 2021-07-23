import math
import numpy as np
#Author: "Salvatre Guccione"
#First Version Released: "16 April 2021"

state=np.empty(5)
p_default=1e5

def setState_pTX(p, T):
    state=np.empty(5)
    d = rho_T(T)
    h = h_T(T)
    u = h - p/d
    state[0] = T
    state[1] = p
    state[2] = d
    state[3] = h
    state[4] = u
    #MM = 0.07862523372 #Composition: NaCl–KCl–MgCl2 24.5–20.5–55 (% in weight)
    #R = 8.3144 / MM
    #T = T_h(h)
    #h = state.h
    return state

def setState_phX(p, h):
    state=np.empty(5)
    T=T_h(h)
    d = rho_T(T)
    u = h - p/d
    state[0] = T
    state[1] = p
    state[2] = d
    state[3] = h
    state[4] = u
    #MM = 0.07862523372 #Composition: NaCl–KCl–MgCl2 24.5–20.5–55 (% in weight)
    #R = 8.3144 / MM
    #T = T_h(h)
    #h = state.h
    return state

def h_T(T):
    a = -0.2241891775
    b = 1411.9325606493
    c = -832203.52407882
    #h is obtained by integrating (cp dT). The integration constant was added such that the h value at T = 298.15K (i.e. 25 degC) becomes zero.
    h = a * T ** 2 + b * T + c
    return h

def T_h(h):
    a = 0.0009793714
    b = 651.8423145829
    T = a*h + b
    return T

def rho_T(T):
    a = -0.5786666667
    b = 2124.1516888889
    rho = a*T + b #Constant bulk density
    return rho

def cp_T(T):
    a = -0.448
    b = 1411.6156444445
    cp = a * T + b
    return cp

def eta_T(T):
    a = 3.12354312353038e-14
    b = -1.281501061896e-10
    c = 2.01613494565798e-7
    d = -0.0001466121
    e = 0.044260166
    #From interpolation of NREL data
    eta = a*T**4 + b*T**3 + c*T**2 + d*T + e
    return eta

def k_T(T):
    a = 7.06493506493493e-7
    b = -0.0014404163
    c = 1.1483003444
    k = a*T**2 + b*T + c #Constant thermal conductivity
    return k

def temperature(state):
    T=state[0]
    return T

def pressure(state):
    p=state[1]
    return p

def density(state):
    d=state[2]
    return d

def specificEnthalpy(state):
    h=state[3]
    return h

def specificInternalEnergy(state):
    u=state[4]
    return u

def specificHeatCapacityCp(state):
    T=state[0]
    cp=cp_T(T)
    return cp

def dynamicViscosity(state):
    T=state[0]
    eta=eta_T(T)
    return eta

def thermalConductivity(state):
    T=state[0]
    k=k_T(T)
    return k



