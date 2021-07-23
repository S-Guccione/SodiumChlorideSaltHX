function [U, h_s, h_t] = HTCs(d_o, N_p, N_sp, layout, N_t, m_flow_Na, m_flow_MS, Tm_MS, Tm_Na, l_b);

%% Parameters
R_ss=8.808e-5;
t_tube=TubeThickness(d_o);
d_i=d_o-2*t_tube;
P_t=1.25*d_o;
B=0.25;
N_ss=0.2;
Tm_wall=(Tm_MS+Tm_Na)/2;

%% Tube Side Heat Transfer Coefficient
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
Pr_Na=mu_Na*cp_Na/k_Na;
Pe_Na=Re_Na*Pr_Na;

if (Re_Na>0)  
    if Pe_Na<=1000  
      A=4.5;
    elseif Pe_Na>=2000  
      A=3.6;
    else
      A=5.4-9e-4*Pe_Na;
    end
    Nu_Na=A+0.018*Pe_Na;
    h_t=Nu_Na*k_Na/d_i;
  else
    h_t=0;
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
mu_MS = mu_a*Tm_MS^4 + mu_b*Tm_MS^3 + mu_c*Tm_MS^2 + mu_d*Tm_MS + mu_e;
mu_MS_wall= mu_a*Tm_Na^4 + mu_b*Tm_Na^3 + mu_c*Tm_Na^2 + mu_d*Tm_Na + mu_e;

[L_bb, D_b, D_s] = ShellDiameter(d_o, N_t, layout, N_p);  
S_m=(l_b/N_sp)*(L_bb+(D_b/P_t)*(P_t-d_o));
v_max_MS=m_flow_MS/rho_MS/S_m;
Re_MS=rho_MS*d_o*v_max_MS/mu_MS;
Pr_MS=mu_MS*cp_MS/k_MS;
if (Re_MS>0)  
    if layout==1  
        if (Re_MS<=300)  
            aa=0.742;
            mm=0.431;
        elseif (Re_MS>300) && (Re_MS<2e5)   
            aa=0.211;
            mm=0.651;
        elseif (Re_MS>2e5) && (Re_MS<2e6)  
            aa=0.116;
            mm=0.7;
        end  
    else
        if (Re_MS<=300)  
            aa=1.309;
            mm=0.36;
        elseif (Re_MS>300) && (Re_MS<2e5)   
            aa=0.273;
            mm=0.635;
        elseif (Re_MS>2e5) && (Re_MS<2e6)  
            aa=0.124;
            mm=0.7;
        end  
    end
    Nu_MS=aa*(Re_MS^mm)*(Pr_MS^0.34)*((mu_MS/mu_MS_wall)^0.26);
    h_s_id=Nu_MS*k_MS/d_o;
    L_c=B*D_s;
    theta_ctl=2*acos((D_s-2*L_c)/D_b);
    F_w=theta_ctl/(2*pi)-sin(theta_ctl)/(2*pi);
    F_c=1-2*F_w;
    J_C=0.55+0.72*F_c;
    L_sb=(3.1+0.004*D_s)/1000;
    theta_ds=2*acos(1-2*B);
    S_sb=(D_s/N_sp)*(pi/2)*L_sb*((2*pi-theta_ds)/(2*pi));
    L_tb=0.0008;
    S_tb=(1/N_sp)*(pi/4)*((d_o+L_tb)^2-d_o^2)*N_t*(1-F_w);
    r_lm=(S_sb+S_tb)/S_m;
    r_s=S_sb/(S_sb+S_tb);
    xx=-0.15*(1+r_s)+0.8;
    J_L=0.44*(1-r_s)+(1-0.44*(1-r_s))*exp(-2.2*r_lm);
    S_b=L_bb*l_b;
    F_bp=S_b/S_m;
    if layout==1  
        N_c=ceil(D_s*(1-2*L_c/D_s)/P_t);
    else
        N_c=ceil(D_s*(1-2*L_c/D_s)/P_t/0.866);
    end  
    r_ss=N_ss/N_c;
    J_B=exp(-1.35*F_bp*(1-(2*r_s)^(1/3)));
    h_s=h_s_id*J_C*J_L*J_B;  
else
    h_s=0;
end  
  
  %% Global heat transfer coefficient:
  rho_wall = 8970;
  k_wall = 0.01996*Tm_wall + 2.981;
  if (Re_Na==0) || (Re_MS==0)  
    U=0;
  else
    U=(1/h_s+R_ss+1/h_t*d_o/d_i+d_o*0.5/k_wall*log(d_o/d_i))^(-1);
  end

end