import numpy as np
import Models.Utilities.U_HeatExchanger.Design_HX as DHX
#Author: "Salvatre Guccione"
#First Version Released: "16 April 2021"

def Optimize_HX(Q_d_des, T_Na1_des, T_Na2_des, T_MS1_des, T_MS2_des, p_Na1_des, p_MS1_des, c_e, r, H_y, n, material_sc, ratio_max, ratio_cond, L_max_cond, L_max_input, v_max_MS_lim_min, v_max_MS_lim_max, v_Na_lim_min, v_Na_lim_max, M_conv_currency_to_USD):

    ############################################################################################### Sweep Parameters ####################################################################################################################################
    # 1) Outer Tube Diameter
    # d_o=[6.35e-3, 9.53e-3, 12.70e-3, 15.88e-3, 19.05e-3, 22.23e-3, 25.40e-3, 28.58e-3, 31.75e-3, 34.93e-3, 38.10e-3, 41.28e-3, 44.45e-3, 47.63e-3, 50.80e-3, 53.98e-3, 57.15e-3, 60.33e-3, 63.50e-3] #Outer Tube Diameter
    d_o=[6.35e-3, 9.53e-3, 19.05e-3, 25.40e-3, 31.75e-3, 41.28e-3, 53.98e-3, 63.50e-3] #Outer Tube Diameter
    
    # 2) Tube passes number
    # N_p = [1,2] #Tube passes number
    N_p = [1] #Tube passes number

    # 3) Tube layout
    # layout = [1,2] #Tube layout
    layout= [2] #Tube layout

    #Auxiliary parameters
    num_dim = 3
    dim1 = len(d_o)
    dim2 = len(N_p)
    dim3 = len(layout)
    dim_tot = dim1 * dim2 * dim3
    N_t_input_on=False #Activate maximum HX length constraint
    N_t_input=1 #Number of tubes
    vec=np.empty((dim_tot, num_dim))

    #Design Parameters
    m_flow_Na_vec=np.empty(dim_tot)
    m_flow_MS_vec=np.empty(dim_tot)
    F_vec=np.empty(dim_tot)
    UA_vec=np.empty(dim_tot)
    A_vec=np.empty(dim_tot)
    U_vec=np.empty(dim_tot)
    N_t_vec=np.empty(dim_tot)
    Dp_tube_vec=np.empty(dim_tot)
    Dp_shell_vec=np.empty(dim_tot)
    h_s_vec=np.empty(dim_tot)
    h_t_vec=np.empty(dim_tot)
    D_s_vec=np.empty(dim_tot)
    D_s_out_vec=np.empty(dim_tot)
    N_baffles_vec=np.empty(dim_tot)
    l_b_vec=np.empty(dim_tot)
    v_Na_vec=np.empty(dim_tot)
    v_max_MS_vec=np.empty(dim_tot)
    V_HX_vec=np.empty(dim_tot)
    m_HX_vec=np.empty(dim_tot)
    m_material_HX_vec=np.empty(dim_tot)
    C_BEC_vec=np.empty(dim_tot)
    C_pump_vec=np.empty(dim_tot)
    TAC_vec=np.empty(dim_tot)
    ex_eff_vec=np.empty(dim_tot)
    en_eff_vec=np.empty(dim_tot)
    L_vec=np.empty(dim_tot)
    ratio_vec=np.empty(dim_tot)
    penalty_vec=np.empty(dim_tot)
    N_t_min_vec=np.empty(dim_tot)
    N_t_max_vec=np.empty(dim_tot)
    iter = 0
    for ww in range(0,dim3,1):
        for ii in range(0,dim2,1):
            for kk in range(0,dim1,1):
                vec[iter, 0] = d_o[kk]
                vec[iter, 1] = N_p[ii]
                vec[iter, 2] = layout[ww]
                (m_flow_Na_vec[iter], m_flow_MS_vec[iter], F_vec[iter], UA_vec[iter], A_vec[iter], U_vec[iter], N_t_vec[iter], Dp_tube_vec[iter], Dp_shell_vec[iter], h_s_vec[iter], h_t_vec[iter], D_s_vec[iter], D_s_out_vec[iter], N_baffles_vec[iter], l_b_vec[iter], v_Na_vec[iter], v_max_MS_vec[iter], V_HX_vec[iter], m_HX_vec[iter], m_material_HX_vec[iter], C_BEC_vec[iter], C_pump_vec[iter], TAC_vec[iter], ex_eff_vec[iter], en_eff_vec[iter], L_vec[iter], ratio_vec[iter], penalty_vec[iter], N_t_min_vec[iter], N_t_max_vec[iter]) = DHX.Design_HX(Q_d_des, T_Na1_des, T_MS1_des, T_MS2_des, d_o[kk], N_p[ii], N_p[ii], layout[ww], T_Na2_des, p_MS1_des, p_Na1_des, c_e, r, H_y, n, material_sc, ratio_max, ratio_cond, L_max_cond, L_max_input, N_t_input_on, N_t_input, v_max_MS_lim_min, v_max_MS_lim_max, v_Na_lim_min, v_Na_lim_max, M_conv_currency_to_USD)
                iter = iter + 1
    
    TAC=np.amin(TAC_vec)
    result=np.amin(np.where(TAC_vec == TAC))
    A_tot = A_vec[result]
    U = U_vec[result]
    N_t = N_t_vec[result]
    m_flow_Na = m_flow_Na_vec[result]
    m_flow_MS = m_flow_MS_vec[result]
    F = F_vec[result]
    UA = UA_vec[result]
    L=L_vec[result]
    Dp_tube = Dp_tube_vec[result]
    Dp_shell = Dp_shell_vec[result]
    h_s = h_s_vec[result]
    h_t = h_t_vec[result]
    D_s = D_s_vec[result]
    D_s_out = D_s_out_vec[result]
    N_baffles = N_baffles_vec[result]
    l_b=l_b_vec[result]
    v_Na = v_Na_vec[result]
    v_max_MS = v_max_MS_vec[result]
    V_HX = V_HX_vec[result]
    m_HX = m_HX_vec[result]
    m_material_HX = m_material_HX_vec[result]
    C_BEC = C_BEC_vec[result]
    C_pump = C_pump_vec[result]
    ex_eff= ex_eff_vec[result]
    en_eff = en_eff_vec[result]
    ratio= ratio_vec[result]
    penalty=penalty_vec[result]
    d_o_opt = vec[result, 0]
    N_p_opt = (vec[result, 1])
    N_sp_opt = (vec[result, 1])
    layout_opt = (vec[result, 2])
    N_t_min=N_t_min_vec[result]
    N_t_max=N_t_max_vec[result]
    return m_flow_Na, m_flow_MS, F, UA, A_tot, U, N_t, Dp_tube, Dp_shell, h_s, h_t, D_s, D_s_out, N_baffles, l_b, v_Na, v_max_MS, V_HX, m_HX, m_material_HX, C_BEC, C_pump, TAC, ex_eff, en_eff, L, ratio, penalty, d_o_opt, N_p_opt, N_sp_opt, layout_opt, N_t_min, N_t_max

