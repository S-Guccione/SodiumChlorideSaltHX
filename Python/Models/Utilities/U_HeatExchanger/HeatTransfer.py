import numpy as np
import math as MA
import Models.Utilities.U_HeatExchanger.Utilities_HX as UHX
import Models.Media.Sodium as Medium1
import Models.Media.ChlorideSalt as Medium2
#Author: "Salvatre Guccione"
#First Version Released: "16 April 2021"

def HTCs(d_o, N_p, N_sp, layout, N_t, state_mean_Na, state_mean_MS, state_wall_MS, m_flow_Na, m_flow_MS, l_b):
    R_ss = 8.808e-5
    t_tube=UHX.TubeThickness(d_o=d_o)
    d_i = d_o - 2 * t_tube
    P_t=1.25*d_o
    B=0.25
    #N_ss=0.2

    Tm_MS = Medium2.temperature(state_mean_MS)
    Tm_Na = Medium1.temperature(state_mean_Na)
    Tm_wall = (Tm_MS+Tm_Na)/2
    [k_wall, _] = UHX.Haynes230_BaseProperties(Tm_wall)
    
    #Sodium properties
    #rho_Na = Medium1.density(state_mean_Na)
    cp_Na = Medium1.specificHeatCapacityCp(state_mean_Na)
    mu_Na = Medium1.dynamicViscosity(state_mean_Na)
    #mu_Na_wall = mu_Na
    k_Na = Medium1.thermalConductivity(state_mean_Na)

    #Chloride Salt properties
    rho_MS = Medium2.density(state_mean_MS)
    cp_MS = Medium2.specificHeatCapacityCp(state_mean_MS)
    mu_MS = Medium2.dynamicViscosity(state_mean_MS)
    mu_MS_wall = Medium2.dynamicViscosity(state_wall_MS)
    k_MS = Medium2.thermalConductivity(state_mean_MS)

    #Tube Side Heat Transfer Coefficient
    Tep=MA.ceil(N_t/N_p)
    A_cs=np.pi/4*d_i**2
    A_cs_tot=Tep*A_cs
    M_Na=m_flow_Na/A_cs_tot
    Re_Na=M_Na*d_i/mu_Na
    Pr_Na=mu_Na*cp_Na/k_Na
    Pe_Na=Re_Na*Pr_Na
    if (Pe_Na>0):
        if Pe_Na<=1000:
            A=4.5
        elif Pe_Na>=2000:
            A=3.6
        else:
            A=5.4-9e-4*Pe_Na
        Nu_Na=A+0.018*Pe_Na
        h_t=Nu_Na*k_Na/d_i
    else:
        h_t=0
    
    #Shell-side heat transfer coefficient:
    [L_bb, D_b, D_s, _] = UHX.ShellDiameter(d_o=d_o, N_t=N_t, layout=layout, N_p=N_p)  
    S_m=(l_b/N_sp)*(L_bb+(D_b/P_t)*(P_t-d_o))
    v_max_MS=m_flow_MS/rho_MS/S_m
    Re_MS=rho_MS*d_o*v_max_MS/mu_MS
    Pr_MS=mu_MS*cp_MS/k_MS
    if (Re_MS>0) and (Pr_MS>0):
        if layout==1:
            if (Re_MS<=300): 
                aa=0.742
                mm=0.431
            elif (Re_MS>300) and (Re_MS<2e5):
                aa=0.211
                mm=0.651
            elif (Re_MS>2e5) and (Re_MS<2e6) :
                aa=0.116
                mm=0.7
            #new
            else:
                aa=0.116
                mm=0.7
        else:
            if (Re_MS<=300) :
                aa=1.309
                mm=0.36
            elif (Re_MS>300) and (Re_MS<2e5) : 
                aa=0.273
                mm=0.635
            elif (Re_MS>2e5) and (Re_MS<2e6) :
                aa=0.124
                mm=0.7
            #new
            else:
                aa=0.124
                mm=0.7
        Nu_MS=aa*(Re_MS**mm)*(Pr_MS**0.34)*((mu_MS/mu_MS_wall)**0.26)
        h_s_id=Nu_MS*k_MS/d_o
        L_c=B*D_s
        theta_ctl=2*MA.acos((D_s-2*L_c)/D_b)
        F_w=theta_ctl/(2*np.pi)-MA.sin(theta_ctl)/(2*np.pi)
        F_c=1-2*F_w
        J_C=0.55+0.72*F_c
        L_sb=(3.1+0.004*D_s)/1000
        theta_ds=2*MA.acos(1-2*B)
        S_sb=(D_s/N_sp)*(np.pi/2)*L_sb*((2*np.pi-theta_ds)/(2*np.pi))
        L_tb=0.0008
        S_tb=(1/N_sp)*(np.pi/4)*((d_o+L_tb)**2-d_o**2)*N_t*(1-F_w)
        r_lm=(S_sb+S_tb)/S_m
        r_s=S_sb/(S_sb+S_tb)
        #xx=-0.15*(1+r_s)+0.8
        J_L=0.44*(1-r_s)+(1-0.44*(1-r_s))*MA.exp(-2.2*r_lm)
        S_b=L_bb*l_b
        F_bp=S_b/S_m
        # if layout==1:
        #     N_c=MA.ceil(D_s*(1-2*L_c/D_s)/P_t)
        # else:
        #     N_c=MA.ceil(D_s*(1-2*L_c/D_s)/P_t/0.866)
        #r_ss=N_ss/N_c
        J_B=MA.exp(-1.35*F_bp*(1-(2*r_s)**(1/3)))
        h_s=h_s_id*J_C*J_L*J_B  
    else:
        h_s=0

    #Global heat transfer coefficient:
    if (h_s==0) or (h_t==0):
        U=0
    else:
        U=(1/h_s+R_ss+1/h_t*d_o/d_i+d_o*0.5/k_wall*MA.log(d_o/d_i))**(-1)

    return U, h_s, h_t