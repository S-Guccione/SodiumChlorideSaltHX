function [Dp_tube, Dp_shell, v_Na, v_max_MS]=Dp_losses(d_o, N_p, N_sp, layout, N_t, L, m_flow_Na, m_flow_MS, Tm_MS, Tm_Na, l_b, N);

%% Parameters
t_tube=TubeThickness(d_o);
d_i=d_o-2*t_tube;
P_t=1.25*d_o;
B=0.25;
L_tb=0.0008;
N_ss=0.2;

%% Tube-side pressure drop:
Tm_wall=(Tm_MS+Tm_Na)/2;
rho_Na = 219 + 275.32 * (1 - Tm_Na / 2503.7) + 511.58 * sqrt(1 - Tm_Na / 2503.7);
cp_Na = 1000 * (1.6582 - 8.4790e-4 * Tm_Na + 4.4541e-7 * Tm_Na^2 - 2992.6*Tm_Na^(-2));
mu_Na = exp(-6.4406 - 0.3958 * log(Tm_Na) + 556.835 / Tm_Na);
k_Na = 124.67 - 0.11381 *Tm_Na + 5.5226e-5 * Tm_Na^2 - 1.1842e-8*Tm_Na^3;
mu_Na_wall=mu_Na;

Tep=ceil(N_t/N_p);
A_cs=pi/4*d_i^2;
A_cs_tot=Tep*A_cs;
M_Na=m_flow_Na/A_cs_tot;
v_Na=M_Na/rho_Na;
Re_Na=M_Na*d_i/mu_Na;
  if (Re_Na>0)  
    if (Re_Na<=855)  
      j_f=8.1274*Re_Na^(-1.011);
    else
      j_f=0.046*Re_Na^(-0.244);
    end  
    if (Re_Na<=2100)  
      m=0.25;
    else
      m=0.14;
    end  
    Dp_tube=(N_p*(2.5+8*j_f*(L/d_i)*(mu_Na/mu_Na_wall)^(-m)))*0.5*rho_Na*v_Na^2;
  else
    Dp_tube=0;
  end  
  
  %% Shell-side heat transfer coefficient:
cp_MS = -0.448*Tm_MS+1411.6156444445;
rho_MS = 2124.1516888889-0.5786666667*Tm_MS;
k_MS = 7.06493506493493e-7*Tm_MS^2 -0.0014404163*Tm_MS + 1.1483003444;
mu_a = 3.12354312353038e-14;
mu_b = -1.281501061896e-10;
mu_c = 2.01613494565798e-7;
mu_d = -0.0001466121;
mu_e = 0.044260166;
mu_MS= mu_a*Tm_MS^4 + mu_b*Tm_MS^3 + mu_c*Tm_MS^2 + mu_d*Tm_MS + mu_e;
mu_MS_wall= mu_a*Tm_Na^4 + mu_b*Tm_Na^3 + mu_c*Tm_Na^2 + mu_d*Tm_Na + mu_e;

[L_bb, D_b, D_s] = ShellDiameter(d_o, N_t, layout, N_p);  
L_c=B*D_s;
S_m=(l_b/N_sp)*(L_bb+(D_b/P_t)*(P_t-d_o));
v_max_MS=m_flow_MS/rho_MS/S_m;
Re_MS=rho_MS*d_o*v_max_MS/mu_MS;
if (Re_MS>0)  
    if layout==1  
        N_c=ceil(D_s*(1-2*L_c/D_s)/P_t);
    else
        N_c=ceil(D_s*(1-2*L_c/D_s)/P_t/0.866);
    end  
if layout==1  
    N_cw=ceil(0.8/P_t*(L_c-(D_s-D_b)/2));
else
    N_cw=ceil(0.8/(0.866*P_t)*(L_c-(D_s-D_b)/2));
end  
if layout==1  
  if (Re_MS<2300)  
    K_f=0.272+(0.207e3/Re_MS)+(0.102e3/Re_MS^2)-(0.286e3/Re_MS^3);
  else
    K_f=0.267+(0.249e4/Re_MS)-(0.927e7/Re_MS^2)+(10^10/Re_MS^3);
  end  
else
  if (Re_MS>4000)  
    K_f=0.245+(0.339e4/Re_MS)-(0.984e7/Re_MS^2)+(0.133e11/Re_MS^3)-(0.599e13/Re_MS^4);
  else
    K_f=11.474*Re_MS^(-0.34417);
  end   
end

    Dp_bi=N_c*K_f*0.5*rho_MS*v_max_MS^2;
    S_b=L_bb*l_b;
    L_sb=(3.1+0.004*D_s)/1000;
    theta_ds=2*acos(1-2*B);
    S_sb=(D_s/N_sp)*(pi/2)*L_sb*((2*pi-theta_ds)/(2*pi));
    theta_ctl=2*acos((D_s-2*L_c)/D_b);
    F_w=theta_ctl/(2*pi)-sin(theta_ctl)/(2*pi);
    F_c=1-2*F_w;
    S_tb=(1/N_sp)*(pi/4)*((d_o+L_tb)^2-d_o^2)*N_t*(1-F_w);
    r_s=S_sb/(S_sb+S_tb);
    r_lm=(S_sb+S_tb)/S_m;
    xx=-0.15*(1+r_s)+0.8;
    R_L=exp(-1.23*(1+r_s))*r_lm^xx;
    F_bp=S_b/S_m;
    r_ss=N_ss/N_c;
    R_B=exp(-3.7*F_bp*(1-r_ss^(1/3)));
    Dp_c=Dp_bi*(N-1)*R_B*R_L;
    S_wg=(pi/4)*(D_s^2/N_sp)*(theta_ds/(2*pi)-sin(theta_ds)/(2*pi));
    S_wt=(N_t/N_sp)*F_w*(pi/4)*d_o^2;
    S_w=S_wg-S_wt;
    Dp_w=(0.2+0.6*N_cw)/(2*S_m*S_w*rho_MS)*m_flow_MS^2;
    Dp_e=2*Dp_bi*R_B*(1+N_cw/N_c);
    Dp_shell=N_sp*(((N-1)*Dp_bi*R_B+N*Dp_w)*R_L+Dp_e);
  else
    Dp_shell=0;
end


  end
  