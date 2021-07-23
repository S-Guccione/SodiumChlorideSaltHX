import numpy as np
import math as MA
import Models.Utilities.U_HeatExchanger.Utilities_HX as UHX
import Models.Media.Sodium as Medium1
import Models.Media.ChlorideSalt as Medium2
import Models.Utilities.U_HeatExchanger.HeatTransfer as HT
import Models.Utilities.U_HeatExchanger.PressureLosses as PL
#Author: "Salvatre Guccione"
#First Version Released: "16 April 2021"

def Design_HX(Q_d, T_Na1, T_MS1, T_MS2, d_o, N_p, N_sp, layout, T_Na2, p_Na1, p_MS1, c_e, r, H_y, n, material_sc, ratio_max, ratio_cond, L_max_cond, L_max_input, N_t_input_on, N_t_input, v_max_MS_lim_min, v_max_MS_lim_max, v_Na_lim_min, v_Na_lim_max, M_conv_currency_to_USD):
    #Parameters
    tol=1e-3 #Heat transfer coefficient tollerance# 
    tol2=1e-2 #Heat transfer coefficient tollerance# 
    iter_max=10 
    iter_max_2=10

    #Velocity limits
    # v_max_MS_lim_min=0.50
    # v_max_MS_lim_max=1.50
    # v_Na_lim_min=4/3.281
    # v_Na_lim_max=8/3.281
    max_length_N_b=10
    l_vec2=5 

    vec_Nb=np.empty(max_length_N_b)
    vec_lb=np.empty(max_length_N_b)
    vec_tb=np.empty(max_length_N_b)
    vec_U=np.empty(max_length_N_b)
    vec_h_s=np.empty(max_length_N_b)
    vec_h_t=np.empty(max_length_N_b)
    vec_condition=np.empty(max_length_N_b)

    #Design Parameters_des_L
    m_flow_Na_vec_des_L=np.empty(l_vec2)
    m_flow_MS_vec_des_L=np.empty(l_vec2)
    F_vec_des_L=np.empty(l_vec2)
    UA_vec_des_L=np.empty(l_vec2)
    A_vec_des_L=np.empty(l_vec2)
    U_vec_des_L=np.empty(l_vec2)
    N_t_vec_des_L=np.empty(l_vec2)
    Dp_tube_vec_des_L=np.empty(l_vec2)
    Dp_shell_vec_des_L=np.empty(l_vec2)
    h_s_vec_des_L=np.empty(l_vec2)
    h_t_vec_des_L=np.empty(l_vec2)
    D_s_vec_des_L=np.empty(l_vec2)
    D_s_out_vec_des_L=np.empty(l_vec2)
    N_baffles_vec_des_L=np.empty(l_vec2)
    l_b_vec_des_L=np.empty(l_vec2)
    v_Na_vec_des_L=np.empty(l_vec2)
    v_max_MS_vec_des_L=np.empty(l_vec2)
    V_HX_vec_des_L=np.empty(l_vec2)
    m_HX_vec_des_L=np.empty(l_vec2)
    m_material_HX_vec_des_L=np.empty(l_vec2)
    C_BEC_vec_des_L=np.empty(l_vec2)
    C_pump_vec_des_L=np.empty(l_vec2)
    TAC_vec_des_L=np.empty(l_vec2)
    ex_eff_vec_des_L=np.empty(l_vec2)
    en_eff_vec_des_L=np.empty(l_vec2)
    L_vec_des_L=np.empty(l_vec2)
    ratio_vec_des_L=np.empty(l_vec2)
    penalty_vec_des_L=np.empty(l_vec2)

    #Shell-side
    B=0.25

    #Cost Function
    mmm=0.37 
    c2=11 
    #M_conv=1 #Conversion factor# 
    eta_pump=0.75 #Pump efficiency# 
    F_ma_min=1.65 #Manufacturing Factor# 

    #Turton Cost Function
    CEPCI_01=397 #CEPCI 2001
    CEPCI_18=603.1 #CEPCI 2018

    if N_t_input_on:
        l_vec1=1
    else:
        l_vec1=2000
    #l_vec=l_vec1*l_vec2
    
    #Design Parameters
    m_flow_Na_vec_des=np.empty(l_vec1)
    m_flow_MS_vec_des=np.empty(l_vec1)
    F_vec_des=np.empty(l_vec1)
    UA_vec_des=np.empty(l_vec1)
    A_vec_des=np.empty(l_vec1)
    U_vec_des=np.empty(l_vec1)
    N_t_vec_des=np.empty(l_vec1)
    Dp_tube_vec_des=np.empty(l_vec1)
    Dp_shell_vec_des=np.empty(l_vec1)
    h_s_vec_des=np.empty(l_vec1)
    h_t_vec_des=np.empty(l_vec1)
    D_s_vec_des=np.empty(l_vec1)
    D_s_out_vec_des=np.empty(l_vec1)
    N_baffles_vec_des=np.empty(l_vec1)
    l_b_vec_des=np.empty(l_vec1)
    v_Na_vec_des=np.empty(l_vec1)
    v_max_MS_vec_des=np.empty(l_vec1)
    V_HX_vec_des=np.empty(l_vec1)
    m_HX_vec_des=np.empty(l_vec1)
    m_material_HX_vec_des=np.empty(l_vec1)
    C_BEC_vec_des=np.empty(l_vec1)
    C_pump_vec_des=np.empty(l_vec1)
    TAC_vec_des=np.empty(l_vec1)
    ex_eff_vec_des=np.empty(l_vec1)
    en_eff_vec_des=np.empty(l_vec1)
    L_vec_des=np.empty(l_vec1)
    ratio_vec_des=np.empty(l_vec1)
    penalty_vec_des=np.empty(l_vec1)

    sc_A=material_sc*8.6*1.1

    #Tube-side
    t_tube=UHX.TubeThickness(d_o)
    d_i=d_o - 2 * t_tube 
    A_cs=np.pi * d_i**2/4
    P_t=1.25 * d_o

    #Temperatures
    Tm_Na=(T_Na1 + T_Na2)/2 
    Tm_MS=(T_MS1 + T_MS2)/2 
    Tm_wall=(Tm_MS + Tm_Na)/2
    F=1
    [_, rho_wall]=UHX.Haynes230_BaseProperties(Tm_wall)

    #Sodium Properties
    state_mean_Na=Medium1.setState_pTX(p_Na1, Tm_Na)
    state_input_Na=Medium1.setState_pTX(p_Na1, T_Na1)
    state_output_Na=Medium1.setState_pTX(p_Na1, T_Na2)
    rho_Na=Medium1.density(state_mean_Na) 
    cp_Na=Medium1.specificHeatCapacityCp(state_mean_Na) 
    #mu_Na=Medium1.dynamicViscosity(state_mean_Na) 
    #mu_Na_wall=mu_Na 
    #k_Na=Medium1.thermalConductivity(state_mean_Na) 
    h_Na1=Medium1.specificEnthalpy(state_input_Na) 
    h_Na2=Medium1.specificEnthalpy(state_output_Na)

    #Chloride Salt properties
    state_mean_MS=Medium2.setState_pTX(Medium2.p_default, Tm_MS) 
    state_wall_MS=Medium2.setState_pTX(Medium2.p_default, Tm_Na) 
    state_input_MS=Medium2.setState_pTX(p_Na1, T_MS1)
    state_output_MS=Medium2.setState_pTX(p_Na1, T_MS2)
    rho_MS=Medium2.density(state_mean_MS) 
    cp_MS=Medium2.specificHeatCapacityCp(state_mean_MS) 
    #mu_MS=Medium2.dynamicViscosity(state_mean_MS) 
    #mu_MS_wall=Medium2.dynamicViscosity(state_wall_MS) 
    #k_MS=Medium2.thermalConductivity(state_mean_MS) 
    h_MS1=Medium2.specificEnthalpy(state_input_MS) 
    h_MS2=Medium2.specificEnthalpy(state_output_MS)

    #Temperature Differences
    DT1=T_Na1 - T_MS2
    DT2=T_Na2 - T_MS1
    if abs(DT1 - DT2)<1e-6:
        LMTD=DT1
    else:
        LMTD=(DT1 - DT2) / MA.log(DT1 / DT2)
    m_flow_Na=Q_d / (cp_Na * (T_Na1 - T_Na2)) 
    m_flow_MS=Q_d / (cp_MS * (T_MS2 - T_MS1)) 
    UA=Q_d / (F * LMTD) 
    ex_eff=m_flow_MS * (h_MS2 - h_MS1 - (25 + 273.15) * cp_MS * MA.log(T_MS2 / T_MS1)) / (m_flow_Na * (h_Na1 - h_Na2 - (25 + 273.15) * cp_Na * MA.log(T_Na1 / T_Na2))) 
    if cp_Na*m_flow_Na > cp_MS * m_flow_MS:
        en_eff=(T_MS2-T_MS1)/(T_Na1-T_MS1)
    else:
        en_eff=(T_Na1 - T_Na2)/(T_Na1 - T_MS1)

    #Minimum Nt
    N_t_min=MA.ceil(m_flow_Na*N_p/(rho_Na*A_cs*v_Na_lim_max))
    N_t_max=MA.floor(m_flow_Na*N_p/(rho_Na*A_cs*v_Na_lim_min))
    N_t_input_corr=min(N_t_max,max(N_t_input,N_t_min))*np.ones(l_vec1)
    if N_t_input_on==True:
        N_t_real=N_t_input_corr
    else:
        N_t_real=(np.linspace(N_t_min,N_t_max,l_vec1))

    for qq in range(0,l_vec1,1):
        N_t=round(N_t_real[qq])
        if N_p>1:
            Tep=MA.floor(N_t/N_p)
            N_t=Tep*N_p 
        else:
            Tep=MA.ceil(N_t/N_p)
        [L_bb, D_b, D_s, D_s_out]=UHX.ShellDiameter(d_o, N_t, layout, N_p) 
        l_b_min=m_flow_MS/(rho_MS*v_max_MS_lim_max)/(L_bb+(D_b/P_t)*(P_t-d_o)) 
        l_b_max=m_flow_MS/(rho_MS*v_max_MS_lim_min)/(L_bb+(D_b/P_t)*(P_t-d_o)) 
        l_b=l_b_max 
        t_baffle_max=UHX.BaffleThickness(D_s, l_b_max) 
        t_baffle_min=UHX.BaffleThickness(D_s, l_b_min) 
        L_max_constrain=D_s_out*ratio_max
        if (ratio_cond and L_max_cond):
            L_max=min(L_max_constrain, L_max_input)
        elif L_max_cond:
            L_max=L_max_input
        elif ratio_cond:
            L_max=L_max_constrain
        else:
            L_max=100
        L_input=np.linspace(2,L_max,l_vec2)

        for iu in range(0,l_vec2,1):
            LL=L_input[iu]
            A_st=LL*np.pi*d_o 
            A_tot=A_st*N_t 
            U_calc=UA/A_tot 
            condition=10
            iter_var=0
            while condition>tol and iter_var<iter_max:
                A_tot=UA/U_calc 
                A_st=A_tot/N_t
                L=A_st/(np.pi*d_o) 
                A_st=np.pi*d_o*L 
                A_tot=A_st*N_t 
                [L_bb, D_b, D_s, D_s_out]=UHX.ShellDiameter(d_o, N_t, layout, N_p)
                N_baffles_min=max(1,np.ceil((L/(l_b_max+t_baffle_max)-1)))
                N_baffles_max=max(1,MA.floor((L/(l_b_min+t_baffle_min)-1)))
                N_baffles_vec=np.around(np.linspace(N_baffles_min,N_baffles_max,max_length_N_b))
                skip=0
                for yy in range(0,len(N_baffles_vec),1):
                    N_baffles=N_baffles_vec[yy]
                    geom_error=np.ones(1)*10
                    iter_2=0 
                    l_b_approx=L/(N_baffles+1)
                    while ((geom_error>tol2 or iter_2>iter_max_2)):
                        if N_baffles<1:
                            t_baffle=0
                        else:
                            t_baffle=UHX.BaffleThickness(D_s,l_b_approx)
                        l_b=L/(N_baffles+1)-t_baffle 
                        geom_error=abs(l_b-l_b_approx)/l_b 
                        l_b_approx=l_b 
                        iter_2=iter_2+1 
                    if geom_error<tol2:
                        vec_Nb[yy-skip]=N_baffles
                        vec_lb[yy-skip]=l_b 
                        vec_tb[yy-skip]=t_baffle 
                    else:
                        skip=skip+1
                length_vec_U=len(vec_tb)-skip
                for uu in range(0,length_vec_U,1):
                    N_baffles=vec_Nb[uu] 
                    l_b=vec_lb[uu] 
                    t_baffle=vec_tb[uu] 
                    [U, h_s, h_t]=HT.HTCs(d_o, N_p, N_p, layout, N_t, state_mean_Na, state_mean_MS, state_wall_MS, m_flow_Na, m_flow_MS, l_b) 
                    condition=abs(U*A_tot-UA)/UA
                    vec_U[uu]=U 
                    vec_h_s[uu]=h_s 
                    vec_h_t[uu]=h_t 
                    vec_condition[uu]=condition
                min_condition=np.amin(vec_condition)
                pos_min_condition=np.amin(np.where(vec_condition == min_condition))
                N_baffles=vec_Nb[pos_min_condition] 
                l_b=vec_lb[pos_min_condition] 
                t_baffle=vec_tb[pos_min_condition] 
                U=vec_U[pos_min_condition] 
                h_s=vec_h_s[pos_min_condition] 
                h_t=vec_h_t[pos_min_condition]
                U_calc=U
                condition = vec_condition[pos_min_condition]
                iter_var=iter_var+1
            [Dp_tube, Dp_shell, v_Na, v_max_MS]=PL.Dp_losses(d_o, N_p, N_sp, layout, N_t, L, state_mean_Na, state_mean_MS, state_wall_MS, m_flow_Na, m_flow_MS, l_b, N_baffles) 
            V_ShellThickness=(D_s_out**2-(D_s**2))*np.pi/4*L 
            V_tubes=np.pi*(d_o**2-d_i**2)/4*L*N_t 
            V_baffles=(np.pi*D_s**2)/4*(1-B)*N_baffles*t_baffle+t_baffle*D_s*L*(N_sp-1) 
            V_material=V_ShellThickness+V_tubes+V_baffles 
            V_Na=np.pi/4*(d_i**2)*L*N_t 
            V_MS=(D_s**2-(d_o**2)*N_t)*np.pi/4*L-V_baffles 
            V_HX=V_material+V_MS+V_Na 
            m_Na=V_Na*rho_Na 
            m_MS=V_MS*rho_MS 
            m_material_HX=V_material*rho_wall 
            m_HX=m_material_HX+m_MS+m_Na

            #Turton Cost function
            P_shell=p_MS1 
            P_tubes=p_Na1 
            P_tube_cost=(P_tubes/10**5)-1 
            P_shell_cost=(P_shell/10**5)-1
            if ((P_tube_cost>5 and P_shell_cost>5)or(P_tube_cost<5 and P_shell_cost>5)):
                both=True 
                P_cost=max(P_tube_cost,P_shell_cost) 
            else:
                both=False 
                P_cost=P_tube_cost 
            k1=4.3247 
            k2=-0.3030 
            k3=0.1634 
            if both:
                C1=0.03881 
                C2=-0.11272 
                C3=0.08183 
            else:
                if P_cost<5:
                    C1=0 
                    C2=0 
                    C3=0 
                else:
                    C1=-0.00164 
                    C2=-0.00627 
                    C3=0.0123
            Fp=10**(C1+C2*MA.log10(P_cost)+C3*(MA.log10(P_cost))**2) 
            Fm=3.7 
            B1=1.63 
            B2=1.66 
            if A_tot>1000:
                A_cost=1000 
            elif A_tot<10:
                A_cost=10
            else:
                A_cost=A_tot
            C_p0=10**(k1+k2*MA.log10(A_cost)+k3*(MA.log10(A_cost))**2) 
            C_BM=C_p0*(CEPCI_18/CEPCI_01)*(B1+B2*Fm*Fp) 
            C_BEC_Turton = C_BM*(A_tot/A_cost)**0.7 * M_conv_currency_to_USD

            #Cost Fucntion
            F_ma=F_ma_min+c2*A_tot**(-mmm)
            C_BEC=max(sc_A*A_tot*F_ma,C_BEC_Turton)
            C_pump=c_e*H_y/eta_pump*(m_flow_MS*Dp_shell/rho_MS+m_flow_Na*Dp_tube/rho_Na)/(1000)
            f=(r*(1+r)**n)/((1+r)**n-1) 
            ratio=L/D_s_out
            if L>L_max:
                l_constrain=1
            else:
                l_constrain=0
            if v_max_MS<v_max_MS_lim_min:
                v_constrain1=1
            else:
                v_constrain1=0
            if v_max_MS>v_max_MS_lim_max:
                v_constrain2=1
            else:
                v_constrain2=0
            if v_Na<v_Na_lim_min:
                v_constrain3=1
            else:
                v_constrain3=0
            if v_Na>v_Na_lim_max:
                v_constrain4=1
            else:
                v_constrain4=0
            v_constrain=v_constrain1+v_constrain2+v_constrain3+v_constrain4
            if condition>0.05:
                UA_constrain1=1
            else:
                UA_constrain1=0
            if geom_error>tol2:
                UA_constrain2=1
            else:
                UA_constrain2=0
            UA_constrain=UA_constrain1+UA_constrain2
            if (v_constrain+UA_constrain+l_constrain)>0:
                TAC=10e10 
                C_BEC=10e10 
                penalty=10e10
            elif condition<tol:
                if (C_BEC>0 and C_pump>0):
                    C_BEC=max(sc_A*A_tot*F_ma,C_BEC_Turton)
                    TAC=f*C_BEC+C_pump 
                    penalty=0
                else:
                    TAC=10e10 
                    C_BEC=10e10 
                    penalty=10e10
            else:
                if (C_BEC>0 and C_pump>0):
                    TAC=(f*C_BEC+C_pump) 
                    penalty=(condition*50)*TAC 
                    C_BEC=max(sc_A*A_tot*F_ma,C_BEC_Turton) 
                    TAC=(f*C_BEC+C_pump)+penalty 
                else:
                    TAC=10e10 
                    C_BEC=10e10 
                    penalty=10e10

            m_flow_Na_vec_des_L[iu]=m_flow_Na 
            m_flow_MS_vec_des_L[iu]=m_flow_MS 
            F_vec_des_L[iu]=F 
            UA_vec_des_L[iu]=UA 
            A_vec_des_L[iu]=A_tot 
            U_vec_des_L[iu]=U 
            N_t_vec_des_L[iu]=N_t
            Dp_tube_vec_des_L[iu]=Dp_tube 
            Dp_shell_vec_des_L[iu]=Dp_shell 
            h_s_vec_des_L[iu]=h_s 
            h_t_vec_des_L[iu]=h_t 
            D_s_vec_des_L[iu]=D_s 
            D_s_out_vec_des_L[iu]=D_s_out 
            N_baffles_vec_des_L[iu]=N_baffles 
            l_b_vec_des_L[iu]=l_b 
            v_Na_vec_des_L[iu]=v_Na 
            v_max_MS_vec_des_L[iu]=v_max_MS 
            V_HX_vec_des_L[iu]=V_HX 
            m_HX_vec_des_L[iu]=m_HX 
            m_material_HX_vec_des_L[iu]=m_material_HX 
            C_BEC_vec_des_L[iu]=C_BEC 
            C_pump_vec_des_L[iu]=C_pump 
            TAC_vec_des_L[iu]=TAC 
            ex_eff_vec_des_L[iu]=ex_eff 
            en_eff_vec_des_L[iu]=en_eff 
            L_vec_des_L[iu]=L 
            ratio_vec_des_L[iu]=ratio 
            penalty_vec_des_L[iu]=penalty

        TAC_des_L=np.amin(TAC_vec_des_L)
        result_des_L=np.amin(np.where(TAC_vec_des_L == TAC_des_L))
        A_tot_des_L=A_vec_des_L[result_des_L] 
        U_des_L=U_vec_des_L[result_des_L] 
        N_t_des_L=N_t_vec_des_L[result_des_L] 
        m_flow_Na_des_L=m_flow_Na_vec_des_L[result_des_L] 
        m_flow_MS_des_L=m_flow_MS_vec_des_L[result_des_L] 
        F_des_L=F_vec_des_L[result_des_L] 
        UA_des_L=UA_vec_des_L[result_des_L] 
        L_des_L=L_vec_des_L[result_des_L] 
        Dp_tube_des_L=Dp_tube_vec_des_L[result_des_L] 
        Dp_shell_des_L=Dp_shell_vec_des_L[result_des_L] 
        h_s_des_L=h_s_vec_des_L[result_des_L] 
        h_t_des_L=h_t_vec_des_L[result_des_L] 
        D_s_des_L=D_s_vec_des_L[result_des_L] 
        D_s_out_des_L=D_s_out_vec_des_L[result_des_L] 
        N_baffles_des_L=N_baffles_vec_des_L[result_des_L] 
        l_b_des_L=l_b_vec_des_L[result_des_L] 
        v_Na_des_L=v_Na_vec_des_L[result_des_L] 
        v_max_MS_des_L=v_max_MS_vec_des_L[result_des_L] 
        V_HX_des_L=V_HX_vec_des_L[result_des_L] 
        m_HX_des_L=m_HX_vec_des_L[result_des_L] 
        m_material_HX_des_L=m_material_HX_vec_des_L[result_des_L] 
        C_BEC_des_L=C_BEC_vec_des_L[result_des_L] 
        C_pump_des_L=C_pump_vec_des_L[result_des_L] 
        ex_eff_des_L=ex_eff_vec_des_L[result_des_L] 
        en_eff_des_L=en_eff_vec_des_L[result_des_L] 
        ratio_des_L=ratio_vec_des_L[result_des_L] 
        penalty_des_L=penalty_vec_des_L[result_des_L]

        m_flow_Na_vec_des[qq]=m_flow_Na_des_L
        m_flow_MS_vec_des[qq]=m_flow_MS_des_L
        F_vec_des[qq]=F_des_L 
        UA_vec_des[qq]=UA_des_L 
        A_vec_des[qq]=A_tot_des_L 
        U_vec_des[qq]=U_des_L 
        N_t_vec_des[qq]=N_t_des_L 
        Dp_tube_vec_des[qq]=Dp_tube_des_L 
        Dp_shell_vec_des[qq]=Dp_shell_des_L 
        h_s_vec_des[qq]=h_s_des_L 
        h_t_vec_des[qq]=h_t_des_L 
        D_s_vec_des[qq]=D_s_des_L 
        D_s_out_vec_des[qq]=D_s_out_des_L 
        N_baffles_vec_des[qq]=N_baffles_des_L 
        l_b_vec_des[qq]=l_b_des_L 
        v_Na_vec_des[qq]=v_Na_des_L 
        v_max_MS_vec_des[qq]=v_max_MS_des_L 
        V_HX_vec_des[qq]=V_HX_des_L 
        m_HX_vec_des[qq]=m_HX_des_L 
        m_material_HX_vec_des[qq]=m_material_HX_des_L 
        C_BEC_vec_des[qq]=C_BEC_des_L 
        C_pump_vec_des[qq]=C_pump_des_L 
        TAC_vec_des[qq]=TAC_des_L 
        ex_eff_vec_des[qq]=ex_eff_des_L 
        en_eff_vec_des[qq]=en_eff_des_L 
        L_vec_des[qq]=L_des_L 
        ratio_vec_des[qq]=ratio_des_L 
        penalty_vec_des[qq]=penalty_des_L
    
    TAC=np.amin(TAC_vec_des)
    result=np.amin(np.where(TAC_vec_des == TAC))
    A_tot=A_vec_des[result] 
    U=U_vec_des[result] 
    N_t=N_t_vec_des[result] 
    m_flow_Na=m_flow_Na_vec_des[result] 
    m_flow_MS=m_flow_MS_vec_des[result] 
    F=F_vec_des[result] 
    UA=UA_vec_des[result] 
    L=L_vec_des[result] 
    Dp_tube=Dp_tube_vec_des[result] 
    Dp_shell=Dp_shell_vec_des[result] 
    h_s=h_s_vec_des[result] 
    h_t=h_t_vec_des[result] 
    D_s=D_s_vec_des[result] 
    D_s_out=D_s_out_vec_des[result] 
    N_baffles=N_baffles_vec_des[result] 
    l_b=l_b_vec_des[result] 
    v_Na=v_Na_vec_des[result] 
    v_max_MS=v_max_MS_vec_des[result] 
    V_HX=V_HX_vec_des[result] 
    m_HX=m_HX_vec_des[result] 
    m_material_HX=m_material_HX_vec_des[result] 
    C_BEC=C_BEC_vec_des[result] 
    C_pump=C_pump_vec_des[result] 
    ex_eff=ex_eff_vec_des[result] 
    en_eff=en_eff_vec_des[result] 
    ratio=ratio_vec_des[result] 
    penalty=penalty_vec_des[result]

    return m_flow_Na, m_flow_MS, F, UA, A_tot, U, N_t, Dp_tube, Dp_shell, h_s, h_t, D_s, D_s_out, N_baffles, l_b, v_Na, v_max_MS, V_HX, m_HX, m_material_HX, C_BEC, C_pump, TAC, ex_eff, en_eff, L, ratio, penalty, N_t_min, N_t_max

    # return N_t_real, N_t_min, N_t_max, en_eff, UA, m_flow_MS, m_flow_Na, LMTD, h_MS1, h_MS2, h_Na1, h_Na2, state_input_MS, state_output_MS, L_input, t_baffle_max, t_baffle_min, Tep, U_calc, A_tot, LL, N_baffles_min, N_baffles_max, N_baffles_vec, vec_Nb, vec_lb, vec_tb, length_vec_U, L, vec_U, vec_h_s, vec_h_t, vec_condition, min_condition, pos_min_condition, N_baffles, U_calc, condition, Dp_tube, Dp_shell, v_Na, v_max_MS
    # return m_flow_Na_des_L, m_flow_MS_des_L, F_des_L, UA_des_L, A_tot_des_L, U_des_L, N_t_des_L, Dp_tube_des_L, Dp_shell_des_L, h_s_des_L, h_t_des_L, D_s_des_L, D_s_out_des_L, N_baffles_des_L, l_b_des_L, v_Na_des_L, v_max_MS_des_L, V_HX_des_L, m_HX_des_L, m_material_HX_des_L, C_BEC_des_L, C_pump_des_L, TAC_des_L, ex_eff_des_L, en_eff_des_L, L_des_L, ratio_des_L, penalty_des_L, N_t_min, N_t_max
