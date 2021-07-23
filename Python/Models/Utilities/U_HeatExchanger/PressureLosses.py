import numpy as np
import math as MA
import Models.Utilities.U_HeatExchanger.Utilities_HX as UHX
import Models.Media.Sodium as Medium1
import Models.Media.ChlorideSalt as Medium2
#Author: "Salvatre Guccione"
#First Version Released: "16 April 2021"

def Dp_losses(d_o, N_p, N_sp, layout, N_t, L, state_mean_Na, state_mean_MS, state_wall_MS, m_flow_Na, m_flow_MS, l_b, N):
    t_tube=UHX.TubeThickness(d_o=d_o)
    d_i = d_o - 2 * t_tube
    P_t=1.25*d_o
    B=0.25
    N_ss=0.2
    L_tb=0.0008
    
    #Tm_MS = Medium2.temperature(state_mean_MS)
    #Tm_Na = Medium1.temperature(state_mean_Na)
    #Tm_wall=(Tm_MS+Tm_Na)/2
    #(k_wall, rho_wall) = UHX.Haynes230_BaseProperties(Tm_wall)

    #Sodium properties
    rho_Na = Medium1.density(state_mean_Na)
    #cp_Na = Medium1.specificHeatCapacityCp(state_mean_Na)
    mu_Na = Medium1.dynamicViscosity(state_mean_Na)
    mu_Na_wall = mu_Na
    #k_Na = Medium1.thermalConductivity(state_mean_Na)

    #Chloride Salt properties
    rho_MS = Medium2.density(state_mean_MS)
    #cp_MS = Medium2.specificHeatCapacityCp(state_mean_MS)
    mu_MS = Medium2.dynamicViscosity(state_mean_MS)
    #mu_MS_wall = Medium2.dynamicViscosity(state_wall_MS)
    #k_MS = Medium2.thermalConductivity(state_mean_MS)

    #Tube-side pressure drop:
    Tep=MA.ceil(N_t/N_p)
    A_cs=MA.pi/4*d_i**2
    A_cs_tot=Tep*A_cs
    M_Na=m_flow_Na/A_cs_tot
    Re_Na=M_Na*d_i/mu_Na
    v_Na=M_Na/rho_Na
    if (Re_Na>0):
        if (Re_Na<=855): 
            j_f=8.1274*Re_Na**(-1.011)
        else:
            j_f=0.046*Re_Na**(-0.244)
        
        if (Re_Na<=2100): 
            m=0.25
        else:
            m=0.14
        Dp_tube=(N_p*(2.5+8*j_f*(L/d_i)*(mu_Na/mu_Na_wall)**(-m)))*0.5*rho_Na*v_Na**2
    else:
        Dp_tube=0

    #Shell-side heat transfer coefficient:
    [L_bb, D_b, D_s, _] = UHX.ShellDiameter(d_o=d_o, N_t=N_t, layout=layout, N_p=N_p)  
    L_c=B*D_s
    S_m=(l_b/N_sp)*(L_bb+(D_b/P_t)*(P_t-d_o))
    v_max_MS=m_flow_MS/rho_MS/S_m
    Re_MS=rho_MS*d_o*v_max_MS/mu_MS
    if (Re_MS>0):
        if layout==1:
            N_c=MA.ceil(D_s*(1-2*L_c/D_s)/P_t)
        else:
            N_c=MA.ceil(D_s*(1-2*L_c/D_s)/P_t/0.866)
        if layout==1:
            N_cw=MA.ceil(0.8/P_t*(L_c-(D_s-D_b)/2))
        else:
            N_cw=MA.ceil(0.8/(0.866*P_t)*(L_c-(D_s-D_b)/2))
        if layout==1 :
            if (Re_MS<2300):  
                K_f=0.272+(0.207e3/Re_MS)+(0.102e3/Re_MS**2)-(0.286e3/Re_MS**3)
            else:
                K_f=0.267+(0.249e4/Re_MS)-(0.927e7/Re_MS**2)+(10**10/Re_MS**3)
        else:
            if (Re_MS>4000):
                K_f=0.245+(0.339e4/Re_MS)-(0.984e7/Re_MS**2)+(0.133e11/Re_MS**3)-(0.599e13/Re_MS**4)
            else:
                K_f=11.474*Re_MS**(-0.34417)
        Dp_bi=N_c*K_f*0.5*rho_MS*v_max_MS**2
        S_b=L_bb*l_b
        L_sb=(3.1+0.004*D_s)/1000
        theta_ds=2*MA.acos(1-2*B)
        S_sb=(D_s/N_sp)*(MA.pi/2)*L_sb*((2*MA.pi-theta_ds)/(2*MA.pi))
        theta_ctl=2*MA.acos((D_s-2*L_c)/D_b)
        F_w=theta_ctl/(2*MA.pi)-MA.sin(theta_ctl)/(2*MA.pi)
        #F_c=1-2*F_w
        S_tb=(1/N_sp)*(MA.pi/4)*((d_o+L_tb)**2-d_o**2)*N_t*(1-F_w)
        r_s=S_sb/(S_sb+S_tb)
        r_lm=(S_sb+S_tb)/S_m
        xx=-0.15*(1+r_s)+0.8
        R_L=MA.exp(-1.23*(1+r_s))*r_lm**xx
        F_bp=S_b/S_m
        r_ss=N_ss/N_c
        R_B=MA.exp(-3.7*F_bp*(1-r_ss**(1/3)))
        #Dp_c=Dp_bi*(N-1)*R_B*R_L
        S_wg=(MA.pi/4)*(D_s**2/N_sp)*(theta_ds/(2*MA.pi)-MA.sin(theta_ds)/(2*MA.pi))
        S_wt=(N_t/N_sp)*F_w*(MA.pi/4)*d_o**2
        S_w=S_wg-S_wt
        Dp_w=(0.2+0.6*N_cw)/(2*S_m*S_w*rho_MS)*m_flow_MS**2
        Dp_e=2*Dp_bi*R_B*(1+N_cw/N_c)
        Dp_shell=N_sp*(((N-1)*Dp_bi*R_B+N*Dp_w)*R_L+Dp_e)
    else:
        Dp_shell=0

    return Dp_tube, Dp_shell, v_Na, v_max_MS