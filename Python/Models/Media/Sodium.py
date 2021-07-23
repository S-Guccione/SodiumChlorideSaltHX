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
    #MM=0.02298977
    #R = 8.3144/MM
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
    #MM=0.02298977
    #R = 8.3144/MM
    return state

def h_T(T):
    #Ref. ANL/RE-95/2, pp. 4
    #371K to 2000K liquid on saturation curve:
    h = 1000 * (-365.77 + 1.6582 * T - 4.2395e-4 * T ** 2 + 1.4847e-7 * T ** 3 + 2992.6 / T - 104.90817873321107)
    return h

def T_h(h):
    p = [-0.009242628114334224, -0.2256606342363793, 0.6441022075255126, 2.609381245065044, -12.17897142894775, -21.99755753201832, 498.300362679127, 1204.309372258783]
    p1 = p[0]
    p2 = p[1]
    p3 = p[2]
    p4 = p[3]
    p5 = p[4]
    p6 = p[5]
    p7 = p[6]
    p8 = p[7]
    h_mean=1173255.929754884
    h_std=638348.3064134347
    #correlation based on Ref. ANL/RE-95/2
    x = (h - h_mean) / h_std #h_norm
    T = p1 * x ** 7 + p2 * x ** 6 + p3 * x ** 5 + p4 * x ** 4 + p5 * x ** 3 + p6 * x ** 2 + p7 * x + p8
    return T

def rho_T(T):
    #Ref. ANL/RE-95/2, pp. 20
    # 371K to 2503.7K liquid on saturation curve:
    rho = 219 + 275.32 * (1 - T / 2503.7) + 511.58 * math.sqrt(1 - T / 2503.7)
    return rho

def cp_T(T):
    #Ref. ANL/RE-95/2, pp. 29
    # 371K to 2000K liquid on saturation curve
    cp = 1000 * (1.6582 - 8.4790e-4 * T + 4.4541e-7 * T ** 2 - 2992.6 * T ** (-2))
    return cp

def eta_T(T):
    #Ref. ANL/RE-95/2, pp. 207
    eta = math.exp(-6.4406 - 0.3958 * math.log(max(1,T)) + 556.835 / T)
    return eta

def k_T(T):
    #Ref. ANL/RE-95/2, pp. 181
    k = 124.67 - 0.11381 * T + 5.5226e-5 * T ** 2 - 1.1842e-8 * T ** 3
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