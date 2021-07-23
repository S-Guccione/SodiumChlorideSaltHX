import pytest
import Models.Utilities.U_HeatExchanger.Design_HX as DHX
import Models.Utilities.U_HeatExchanger.Optimize_HX as OHX
import Models.HeatExchangers.HX_SodiumSalt as HX
import numpy as np

#Design Parameters
Q_d_des = 543e6 #Design Heat Flow Rate
T_Na1_des = 740 + 273.15 #Desing Sodium Hot Fluid Temperature#
T_Na2_des = 520 + 273.15 #Optimal outlet sodium temperature#
T_MS1_des = 500 + 273.15 #Desing Molten Salt Cold Fluid Temperature#
T_MS2_des = 720 + 273.15 #Desing Molten Salt Hot Fluid Temperature#
p_Na1_des = 101325 #Design Sodium Inlet Pressure#
p_MS1_des = 101325 #Design Molten Salt Inlet Pressure#

#Common Input Parameters
ratio_cond = True #Activate ratio constraint#  #Default value = true
ratio_max = 10 #Maximum L/D_s ratio# #If ratio_cond = true provide a value (default value = 10)
L_max_cond = False #Activate maximum HX length constraint# #Default value = false
L_max_input = 1 #Maximum HX length# #If L_max_cond = true provide a value (default value = 10)
r = 0.05
H_y = 5600
n = 30
EUR = False #if "EUR" is equal to "True" Euro is the currency adopted, otherwise USD is adopted
M_conv=0.84 #EUR/USD
if EUR:
    c_e = 0.073*M_conv #EUR/kWh
    material_sc = 84*M_conv #EUR/kg
else:
    c_e = 0.073 #USD/kWh
    material_sc = 84 #USD/kg


#Results
m_flow_Na_design = 1971.2931863283993
m_flow_MS_design = 2429.3978314677624
F_design = 1.0
UA_design = 27150000.0
A_HX = 9368.123925650818
U_design = 2897.1877257697315
N_t = 23509.0
Dp_tube_design = 61814.855510959045
Dp_shell_design = 194296.76825766417
h_s_design = 4785.5802078296565
h_t_design = 66120.0530886629
D_s = 1.8149358997901759
D_s_out = 1.830735899790176
N_baffles = 3.0
l_b = 3.3083819361332742
v_Na_design = 1.9136448299258384
v_max_MS_design = 1.204473332834658
V_HX = 35.036198243609626
m_HX = 86716.8806281723
m_material_HX = 53643.412529412584
C_BEC_HX = 15060892.931712547
C_pump_design = 242103.98320429446
TAC = 1221836.6820025896
ex_eff_design = 0.9865612745786051
en_eff_design = 0.9166666666666666
L = 13.309927744533097
ratio_HX = 7.270260962304051
penalty = 0.0
d_o = 0.00953
N_p = 1.0
N_sp = 1.0
layout = 2.0
N_t_min = 18451.0
N_t_max = 36901.0

# Reorganize Results
HX_design=[m_flow_Na_design, m_flow_MS_design, F_design, UA_design, A_HX, U_design, N_t, Dp_tube_design, Dp_shell_design, h_s_design, h_t_design, D_s, D_s_out, N_baffles, l_b, v_Na_design, v_max_MS_design, V_HX, m_HX, m_material_HX, C_BEC_HX, C_pump_design, TAC, ex_eff_design, en_eff_design, L, ratio_HX, penalty, N_t_min, N_t_max] 
HX_optimize=[m_flow_Na_design, m_flow_MS_design, F_design, UA_design, A_HX, U_design, N_t, Dp_tube_design, Dp_shell_design, h_s_design, h_t_design, D_s, D_s_out, N_baffles, l_b, v_Na_design, v_max_MS_design, V_HX, m_HX, m_material_HX, C_BEC_HX, C_pump_design, TAC, ex_eff_design, en_eff_design, L, ratio_HX, penalty, d_o, N_p, N_sp, layout, N_t_min, N_t_max] 

#Design HX Extra Inputs
N_t_input_on = True #Activate fixed number of tubes#
N_t_input = N_t #Input Number of tubes#

# Inputs HX
optimize_and_run=False

def test_optimize_HX():
    assert (np.array(OHX.Optimize_HX(Q_d_des, T_Na1_des, T_Na2_des, T_MS1_des, T_MS2_des, p_Na1_des, p_MS1_des, c_e, r, H_y, n, material_sc, ratio_max, ratio_cond, L_max_cond, L_max_input, EUR)) - np.array(HX_optimize)).all() < 1e-6

def test_design_HX():
    assert (np.array(DHX.Design_HX(Q_d_des, T_Na1_des, T_MS1_des, T_MS2_des, d_o, N_p, N_sp, layout, T_Na2_des, p_MS1_des, p_Na1_des, c_e, r, H_y, n, material_sc, ratio_max, ratio_cond, L_max_cond, L_max_input, N_t_input_on, N_t_input, EUR)) - np.array(HX_design)).all() < 1e-6

def test_HX():
    assert (np.array(HX.DesignHeatExchanger(Q_d_des, T_Na1_des, T_Na2_des, T_MS1_des, T_MS2_des, p_Na1_des, p_MS1_des, optimize_and_run, d_o, N_p, layout, ratio_cond, ratio_max, L_max_cond, L_max_input, N_t_input_on, N_t_input, c_e, r, H_y, n, material_sc, EUR)) - np.array(HX_optimize)).all() < 1e-6