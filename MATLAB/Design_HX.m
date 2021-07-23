function [res, l_vec2]=Design_HX(Q_d, T_Na1, T_MS1, T_MS2, d_o, L, N_p, N_sp, layout, T_Na2, p_MS1, p_Na1, c_e, r, H_y, n, ratio_max, ratio_cond, L_max_cond, L_max_input)

%% Problem
tol=1e-3;
tol2=1e-2;
iter_max=100;
iter_max_2=10;
t_tube=TubeThickness(d_o);
d_i=d_o-2*t_tube;
B=0.25; 
A_cs=pi*(d_i^2)/4;
P_t=1.25*d_o;
l_vec=2000;

%% Cost Function
material_sc=1*84;
moA_HX=8.6*1.1;
mmm=0.37;
c2=11;
eta_pump=0.75;
F_ma_min=1.65;

%% Velocity limits
v_max_MS_lim_min=0.50;
v_max_MS_lim_max=1.5;
v_Na_lim_min=4/3.281;
v_Na_lim_max=8/3.281;

%% Sodium Properties
Tm_Na=(T_Na1+T_Na2)/2;
rho_Na = 219 + 275.32 * (1 - Tm_Na / 2503.7) + 511.58 * sqrt(1 - Tm_Na / 2503.7);
cp_Na = 1000 * (1.6582 - 8.4790e-4 * Tm_Na + 4.4541e-7 * Tm_Na^2 - 2992.6*Tm_Na^(-2));
mu_Na = exp(-6.4406 - 0.3958 * log(Tm_Na) + 556.835 / Tm_Na);
k_Na = 124.67 - 0.11381 *Tm_Na + 5.5226e-5 * Tm_Na^2 - 1.1842e-8*Tm_Na^3;

%% Chloride Salt properties
Tm_MS=(T_MS1+T_MS2)/2;
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

%% Wall Properties
Tm_wall=(Tm_MS+Tm_Na)/2;
rho_wall =8970;
k_wall = 0.01996*Tm_wall + 2.981;
mu_Na_wall=mu_Na;

%% Temperature Differences
DT1=T_Na1-T_MS2;
DT2=T_Na2-T_MS1;
if abs(DT1-DT2)<1e-6
    LMTD=DT1;
else
    LMTD=(DT1-DT2)/log(DT1 / DT2);
end
m_flow_Na=Q_d/(cp_Na*(T_Na1-T_Na2));
m_flow_MS=Q_d/(cp_MS*(T_MS2 - T_MS1));
C_min=min(m_flow_Na*cp_Na,m_flow_MS*cp_MS);
F=1;
UA=Q_d/(F*LMTD);
NTU=UA/C_min;
if (cp_Na*m_flow_Na)>(cp_MS*m_flow_MS)
    en_eff=(T_MS2-T_MS1)./(T_Na1-T_MS1);
else
    en_eff=(T_Na1-T_Na2)./(T_Na1-T_MS1);
end
R=(T_Na1-T_Na2)/(T_MS2-T_MS1);


%% Minimum Nt
N_t_min=ceil(m_flow_Na*N_p/(rho_Na*A_cs*v_Na_lim_max));
N_t_max=floor(m_flow_Na*N_p/(rho_Na*A_cs*v_Na_lim_min));
N_t_vec=round(linspace(N_t_min,N_t_max,l_vec))';
%N_t_vec=(N_t_min:1:N_t_max)';
l_vec=length(N_t_vec);
indice_vecchio=1;
for qq=1:l_vec
    N_t=N_t_vec(qq);
    if N_p>1
        Tep=floor(N_t/N_p);
        N_t=Tep*N_p;
    else
        Tep=ceil(N_t/N_p);
    end
    [L_bb, D_b, D_s] = ShellDiameter(d_o, N_t, layout, N_p);
    L_max_ratio=D_s*ratio_max;
    if ratio_cond && L_max_cond
        L_max=min(L_max_ratio,L_max_input);
    else
        if L_max_cond
            L_max=L_max_input;
        else
            if ratio_cond
                L_max=L_max_ratio; 
            else
                L_max=100;
            end
        end
    end
    l_b_min=m_flow_MS/(rho_MS*v_max_MS_lim_max)/(L_bb+(D_b/P_t)*(P_t-d_o));
    l_b_max=m_flow_MS/(rho_MS*v_max_MS_lim_min)/(L_bb+(D_b/P_t)*(P_t-d_o));
    l_b=l_b_max;
    t_baffle_max= BaffleThickness(D_s,l_b_max);
    t_baffle_min= BaffleThickness(D_s,l_b_min);
    L_input=(linspace(2,L_max,5))';
    clear res_l
    for iu=1:1:length(L_input)
        LL=L_input(iu);
        A_st=LL*(pi*d_o);
        A_tot=A_st*N_t;
        U_calc=UA/A_tot;
        condition=10;
        iter=0;
        while condition>tol && iter<iter_max
            clear config
            A_tot=UA/U_calc;
            A_st=A_tot/N_t;
            L=A_st/(pi*d_o);
            A_st=pi*d_o*L;
            A_tot=A_st*N_t;
            [L_bb, D_b, D_s] = ShellDiameter(d_o, N_t, layout, N_p);
            N_baffles_min=max(1,ceil((L/(l_b_max+t_baffle_max)-1)));
            N_baffles_max=max(1,floor((L/(l_b_min+t_baffle_min)-1)));
            N_baffles_vec=(N_baffles_min:1:N_baffles_max)';
            skip=0;
            check=0;
            for yy=1:length(N_baffles_vec)
                N_baffles=N_baffles_vec(yy);
                geom_error=10;
                iter_2=0;
                l_b_approx=L/(N_baffles+1);
                while geom_error>tol2 || iter_2>iter_max_2
                    if N_baffles<1
                        t_baffle=0;
                    else
                        t_baffle=BaffleThickness(D_s,l_b_approx);
                    end
                    l_b=L/(N_baffles+1)-t_baffle;
                    geom_error=abs(l_b-l_b_approx)/l_b;
                    l_b_approx=l_b;
                    iter_2=iter_2+1;
                end
                if geom_error<tol2
                    config(yy-skip,:)=[N_baffles, l_b, t_baffle];
                    check=check+1;
                else
                    skip=skip+1;
                end
            end
            clear res_for
            if check<1
                U=0;
                h_s=0;
                h_t=0;
                condition=1;
            else
                for uu=1:length(config(:,1))
                    N_baffles=config(uu,1);
                    l_b=config(uu,2);
                    t_baffle=config(uu,3);
                    [U, h_s, h_t]=HTCs(d_o, N_p, N_sp, layout, N_t, m_flow_Na, m_flow_MS, Tm_MS, Tm_Na, l_b);
                    condition=abs(U*A_tot-UA)/UA;
                    res_for(uu,:)=[N_baffles, l_b, t_baffle, U, h_s, h_t, condition];
                end
                [min_condition pos_min_condition]=min(res_for(:,end));
                N_baffles=res_for(pos_min_condition,1);
                l_b=res_for(pos_min_condition,2);
                t_baffle=res_for(pos_min_condition,3);
        %                     [max_baffles pos_max_baffles]=max(config(:,1));
        %                     N_baffles=config(pos_max_baffles,1);
        %                     l_b=config(pos_max_baffles,2);
        %                     t_baffle=config(pos_max_baffles,3);
        %                     [U, h_s, h_t]=HTCs(d_o, N_p, N_sp, layout, N_t, m_flow_Na, m_flow_MS, Tm_MS, Tm_Na, l_b);
        %                     condition=abs(U*A_tot-UA)/UA;
                U=res_for(pos_min_condition,4);
                h_s=res_for(pos_min_condition,5);
                h_t=res_for(pos_min_condition,6);
                condition=res_for(pos_min_condition,7);
            end
            U_calc=U;
            iter=iter+1;
        end
        [Dp_tube, Dp_shell, v_Na, v_max_MS]=Dp_losses(d_o, N_p, N_sp, layout, N_t, L, m_flow_Na, m_flow_MS, Tm_MS, Tm_Na, l_b, N_baffles);
        t_shell=ShellThickness(D_s);
        D_s_out=D_s+2*t_shell;
        V_ShellThickness=(D_s_out^2-(D_s^2))*pi/4*L;
        V_tubes=pi*(d_o^2-d_i^2)/4*L*N_t;
        V_baffles=(pi*D_s^2)/4*(1-B)*N_baffles*t_baffle+t_baffle*D_s*L*(N_sp-1);
        V_material=V_ShellThickness+V_tubes+V_baffles;
        V_Na=pi/4*(d_i^2)*L*N_t;
        V_MS=(D_s^2-(d_o^2)*N_t)*pi/4*L-V_baffles;
        V_HX=V_material+V_MS+V_Na;
        m_Na=V_Na*rho_Na;
        m_MS=V_MS*rho_MS;
        m_material_HX=V_material*rho_wall;
        m_HX=m_material_HX+m_MS+m_Na;

        %% Cost Fucntion
        F_ma=F_ma_min+c2.*A_tot.^(-mmm);
        CEPCI_01=397;
        CEPCI_18=603.1;
        k1=4.3247;
        k2=-0.3030;
        k3=0.1634;
        Fp=1;
        Fm=3.7;
        % Fm=10;
        B1=1.63;
        B2=1.66;
        if A_tot>1000
            A_cost=1000;    
        elseif A_tot<10
            A_cost=10; 
        else
            A_cost=A_tot;    
        end
        C_p0=10.^(k1+k2.*log10(A_cost)+k3.*(log10(A_cost)).^2);
        C_BM=C_p0.*(CEPCI_18/CEPCI_01)*(B1+B2.*Fm.*Fp);
        C_BEC_HX_turton=C_BM;
        C_BEC=max(material_sc*moA_HX*A_tot*F_ma,C_BEC_HX_turton);
        C_pump=c_e*H_y/eta_pump*(m_flow_MS*Dp_shell/rho_MS+m_flow_Na*Dp_tube/rho_Na)/(1000);
        f=(r*(1+r)^n)/((1+r)^n-1);
        ratio=L/D_s_out;
        if L>L_max
            l_constrain=true;
        else
            l_constrain=false;
        end
        v_constrain=(v_max_MS<v_max_MS_lim_min || v_max_MS>v_max_MS_lim_max || v_Na<v_Na_lim_min || v_Na>v_Na_lim_max);
        if (condition>0.01) || (v_constrain) || (geom_error>tol2) || (l_constrain)
            TAC=10e10;
            C_BEC=10e10;
            penalty=10e10;
        elseif condition<tol 
            if C_BEC>0 && C_pump>0
              C_BEC=max(material_sc*moA_HX*A_tot*F_ma,C_BEC_HX_turton);
              TAC=f*C_BEC+C_pump;
              penalty=0;
            else
              TAC=10e10;
              C_BEC=10e10;
              penalty=10e10;
            end 
        else
            if C_BEC>0 && C_pump>0 
              TAC=(f*C_BEC+C_pump);
              penalty=(condition*50)*TAC;
              C_BEC=max(material_sc*moA_HX*A_tot*F_ma,C_BEC_HX_turton);
              TAC=(f*C_BEC+C_pump)+penalty;
            else
              TAC=10e10;
              C_BEC=10e10;
              penalty=10e10;
            end 
        end
        res_l(iu,:)=[m_flow_Na, m_flow_MS, R, C_min, F, UA, en_eff, N_t, N_t_min, N_t_max, Tep, A_tot, A_st, D_s_out, D_s, L, N_baffles, N_baffles_min, N_baffles_max, l_b, t_baffle, t_tube, t_shell, condition, U_calc, h_s, h_t, Dp_tube, Dp_shell, v_Na, v_max_MS, V_HX, m_HX, m_material_HX, F_ma, C_BEC, C_pump, penalty, rho_MS, rho_Na, ratio, TAC];
        ind_res_new=length(res_l(:,1));
    end
    indice_nuovo=(indice_vecchio:1:indice_vecchio-1+ind_res_new)';
    res(indice_nuovo,:)=res_l;
    indice_vecchio=length(res(:,1))+1;
    l_vec2=length(res(:,1));
%     res(qq,:)=[m_flow_Na, m_flow_MS, R, C_min, F, UA, en_eff, N_t, N_t_min, N_t_max, Tep, A_tot, A_st, D_s_out, D_s, L, N_baffles, N_baffles_min, N_baffles_max, l_b, t_baffle, t_tube, t_shell, condition, U_calc, h_s, h_t, Dp_tube, Dp_shell, v_Na, v_max_MS, V_HX, m_HX, m_material_HX, F_ma, C_BEC, C_pump, penalty, rho_MS, rho_Na, ratio, TAC];
end
