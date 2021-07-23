import math as MA
import numpy as np
import Models.Utilities.U_HeatExchanger.Design_HX as DHX
import Models.Utilities.U_HeatExchanger.Optimize_HX as OHX
import Models.Utilities.GeneralUtilities as GU
import Models.Media.Sodium as Medium1
import Models.Media.ChlorideSalt as Medium2
import Models.Utilities.U_HeatExchanger.HeatTransfer as HT
import Models.Utilities.U_HeatExchanger.PressureLosses as PL
#Author: "Salvatre Guccione"
#First Version Released: "16 April 2021"

################################################################################################### Design Mode ###################################################################################################

def DesignHeatExchanger(Q_d_des, T_Na1_des, T_Na2_des, T_MS1_des, T_MS2_des, p_Na1_des, p_MS1_des, optimize_and_run, d_o_input, N_p_input, layout_input, ratio_cond, ratio_max, L_max_cond, L_max_input, N_t_input_on, N_t_input, c_e, r, H_y, n, material_sc, v_max_MS_lim_min, v_max_MS_lim_max, v_Na_lim_min, v_Na_lim_max, M_conv_currency_to_USD):
    if optimize_and_run == True:
        [m_flow_Na_design, m_flow_MS_design, F_design, UA_design, A_HX, U_design, N_t, Dp_tube_design, Dp_shell_design, h_s_design, h_t_design, D_s, D_s_out, N_baffles, l_b, v_Na_design, v_max_MS_design, V_HX, m_HX, m_material_HX, C_BEC_HX, C_pump_design, TAC, ex_eff_design, en_eff_design, L, ratio_HX, penalty, d_o, N_p, N_sp, layout, N_t_min, N_t_max] = OHX.Optimize_HX(Q_d_des, T_Na1_des, T_Na2_des, T_MS1_des, T_MS2_des, p_Na1_des, p_MS1_des, c_e, r, H_y, n, material_sc, ratio_max, ratio_cond, L_max_cond, L_max_input, v_max_MS_lim_min, v_max_MS_lim_max, v_Na_lim_min, v_Na_lim_max, M_conv_currency_to_USD)
    else:
        d_o=d_o_input
        N_p=N_p_input
        N_sp=N_p_input
        layout=layout_input
        [m_flow_Na_design, m_flow_MS_design, F_design, UA_design, A_HX, U_design, N_t, Dp_tube_design, Dp_shell_design, h_s_design, h_t_design, D_s, D_s_out, N_baffles, l_b, v_Na_design, v_max_MS_design, V_HX, m_HX, m_material_HX, C_BEC_HX, C_pump_design, TAC, ex_eff_design, en_eff_design, L, ratio_HX, penalty, N_t_min, N_t_max] = DHX.Design_HX(Q_d_des, T_Na1_des, T_MS1_des, T_MS2_des, d_o, N_p, N_sp, layout, T_Na2_des, p_MS1_des, p_Na1_des, c_e, r, H_y, n, material_sc, ratio_max, ratio_cond, L_max_cond, L_max_input, N_t_input_on, N_t_input, v_max_MS_lim_min, v_max_MS_lim_max, v_Na_lim_min, v_Na_lim_max, M_conv_currency_to_USD)
    HX_design=[m_flow_Na_design, m_flow_MS_design, F_design, UA_design, A_HX, U_design, N_t, Dp_tube_design, Dp_shell_design, h_s_design, h_t_design, D_s, D_s_out, N_baffles, l_b, v_Na_design, v_max_MS_design, V_HX, m_HX, m_material_HX, C_BEC_HX, C_pump_design, TAC, ex_eff_design, en_eff_design, L, ratio_HX, penalty, d_o, N_p, N_sp, layout, N_t_min, N_t_max] 
    return HX_design

################################################################################################### Operating Mode ###################################################################################################

def Operating_HeatExchanger(SF_on, Q_out_rec, m_flow_Na, m_flow_CS, h1_Na, h1_CS, HX_design, k_loss_Na, k_loss_CS):
    d_o=HX_design[28]
    N_p=HX_design[29]
    layout=HX_design[31]
    L=HX_design[25]
    N_t=HX_design[6]
    l_b=HX_design[14]
    N_baffles=HX_design[13]
    A_HX=HX_design[4]
    T_Na2_max = GU.from_degC(700)
    T_Na2_min = GU.from_degC(501)
    state_Na2_max = Medium1.setState_pTX(Medium1.p_default, T_Na2_max)
    state_Na2_min = Medium1.setState_pTX(Medium1.p_default, T_Na2_min)
    h2_Na_max = Medium1.specificEnthalpy(state_Na2_max)
    h2_Na_min = Medium1.specificEnthalpy(state_Na2_min)
    eff_pump=0.75

    if SF_on:
        h2_Na=max(h2_Na_min, min(h2_Na_max, h1_Na-Q_out_rec/max(1e-3,m_flow_Na)))
        h2_CS=h1_CS+Q_out_rec/max(1e-3,m_flow_CS)
        state_input_CS = Medium2.setState_phX(Medium2.p_default, h1_CS)
        T_CS1 = Medium2.temperature(state_input_CS)
        state_input_Na = Medium1.setState_phX(Medium1.p_default, h1_Na)
        T_Na1 = Medium1.temperature(state_input_Na)
        state_output_CS = Medium2.setState_phX(Medium2.p_default, h2_CS)
        T_CS2 = Medium2.temperature(state_output_CS)
        state_output_Na = Medium1.setState_phX(Medium1.p_default, h2_Na)
        T_Na2 = Medium1.temperature(state_output_Na)
        Tm_Na=(T_Na1 + T_Na2) / 2
        state_mean_Na=Medium1.setState_pTX(Medium1.p_default, Tm_Na)
        Tm_CS = (T_CS1 + T_CS2) / 2
        state_mean_CS = Medium2.setState_pTX(Medium2.p_default, Tm_CS)
        state_wall_CS = Medium2.setState_pTX(Medium2.p_default, Tm_Na)
        DT1 = T_Na1 - T_CS2
        DT2 = T_Na2 - T_CS1
        if DT1/DT2 <= 0:
            LMTD = 0
        elif abs(DT1 - DT2) < 1e-3:
            LMTD = DT1
        else:
            LMTD = (DT1 - DT2) / MA.log(DT1 / DT2)
        F=1
        [U, h_s, h_t]=HT.HTCs(d_o, N_p, N_p, layout, N_t, state_mean_Na, state_mean_CS, state_wall_CS, m_flow_Na, m_flow_CS, l_b)
        Q_HX = U * A_HX * F * LMTD
        [Dp_tube, Dp_shell, v_Na, v_max_CS]=PL.Dp_losses(d_o, N_p, N_p, layout, N_t, L, state_mean_Na, state_mean_CS, state_wall_CS, m_flow_Na, m_flow_CS, l_b, N_baffles)
        W_loss_Na=m_flow_Na*k_loss_Na+m_flow_Na*Dp_tube/Medium1.density(state_mean_Na)/eff_pump
        W_loss_CS=m_flow_CS*k_loss_CS+m_flow_CS*Dp_shell/Medium2.density(state_mean_CS)/eff_pump
    else:
        h2_Na=h1_Na
        h2_CS=h1_CS
        U=0
        h_s=0
        h_t=0
        Dp_tube=0
        Dp_shell=0
        v_Na=0
        v_max_CS=0
        Q_HX=0
        W_loss_Na=0
        W_loss_CS=0
    return Q_HX, h2_Na, h2_CS, U, h_s, h_t, Dp_tube, Dp_shell, v_Na, v_max_CS, W_loss_Na, W_loss_CS