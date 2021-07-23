import numpy as np
import math as MA
#Author: "Salvatre Guccione"
#First Version Released: "16 April 2021"

def from_degC(T_degC):
    T=T_degC+273.15
    return T

def to_degC(T_K):
    T_degC=T_K-273.15
    return T_degC

def from_deg(deg):
    rad=deg*np.pi/180
    return rad

def to_deg(rad):
    deg=rad*180/np.pi
    return deg

def Not(y):
    if y:
        z=False
    else:
        z=True
    return z

def Or(y1, y2):
    if y1 or y2:
        z=True
    else:
        z=False
    return z

def And(y1, y2):
    if y1 and y2:
        z=True
    else:
        z=False
    return z

def LMTD_calc(T_hot_in, T_hot_out, T_cold_in, T_cold_out):
    DT1 = T_hot_in - T_cold_out
    DT2 = T_hot_out - T_cold_in
    if DT1/DT2 <= 0:
        LMTD = 0
    elif abs(DT1 - DT2) < 1e-3:
        LMTD = DT1
    else:
        LMTD = (DT1 - DT2) / MA.log(DT1 / DT2)
    return LMTD
