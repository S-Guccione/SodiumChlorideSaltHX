close all
clear
clc

%% Input Parameters
Q_d_vec = 543203279.460279;
T_Na1 = 740+273.15;
T_MS1 = 500+273.15;
T_MS2 = 720+273.15;
p_MS1 = 101325;
p_Na1 = 101325;
c_e = 0.073;
% 0.13/0.9175;
r = 0.05;
H_y = 5600;
% 4500;
n = 30;
ratio_max=10;
ratio_cond=true;
L_max_cond=false;
L_max_input=1;

for qq=1:length(Q_d_vec)
% %% Sweep Parameters
%d_o = [6.35e-3, 9.53e-3, 12.70e-3, 15.88e-3, 19.05e-3, 22.23e-3, 25.40e-3, 28.58e-3, 31.75e-3, 34.93e-3, 38.10e-3, 41.28e-3, 44.45e-3, 47.63e-3, 50.80e-3, 53.98e-3, 57.15e-3, 60.33e-3, 63.50e-3]';
d_o = 9.53e-3;
N_p = [1]';
layout = [2]';
T_Na2 = [520 + 273.15]';
num_dim = 5;
dim1 = size(d_o, 1);
dim2 = size(N_p, 1);
dim3 = size(layout, 1);
dim4 = size(T_Na2, 1);
iter=1;
N_sp=N_p;
results=[];
for ll=1:dim4
    for ww=1:dim3
        for ii=1:dim2
            for kk=1:dim1
                clear vec res l_vec
                d_o_input=d_o(kk);
                N_p_input=N_p(ii);
                N_sp_input=N_sp(ii);
                layout_input=layout(ww);
                T_Na2_input=T_Na2(ll);
                [res,l_vec]=Design_HX(Q_d_vec(qq), T_Na1, T_MS1, T_MS2, d_o_input, 1, N_p_input, N_sp_input, layout_input, T_Na2_input, p_MS1, p_Na1, c_e, r, H_y, n, ratio_max, ratio_cond, L_max_cond, L_max_input);
                vec =[d_o(kk)*ones(l_vec,1), N_p(ii)*ones(l_vec,1), layout(ww)*ones(l_vec,1), T_Na2(ll)*ones(l_vec,1)];
                results=[results; vec,res];
                iter = iter + 1;
            end
        end
    end
end

C_BEC=results(:,end-6);
ind=find(C_BEC<10e10);
results_inter=results(ind,:);

name_results={'d_o'; 'N_p'; 'layout'; 'T_Na2'; 'm_flow_Na'; 'm_flow_MS'; 'R';'C_min'; 'F'; 'UA'; 'en_eff'; 'N_t';'N_t_min';'N_t_max'; 'Tep'; 'A_tot'; 'A_st'; 'D_s_out'; 'D_s'; 'L'; 'N_baffles'; 'N_baffles_min'; 'N_baffles_max'; 'l_b'; 't_baffle'; 't_tube'; 't_shell'; 'condition'; 'U_calc'; 'h_s'; 'h_t'; 'Dp_tube'; 'Dp_shell'; 'v_Na'; 'v_max_MS'; 'V_HX'; 'm_HX'; 'm_material_HX'; 'F_ma'; 'C_BEC'; 'C_pump'; 'penalty'; 'rho_MS'; 'rho_Na'; 'ratio'; 'TAC'};
    for jj=1:length(N_p)
%     for jj=length(N_p)
        ind_Np=find(results_inter(:,2)==N_p(jj));
        mat_Np=results_inter(ind_Np,:);
        TAC=mat_Np(:,end);
        C_BEC=mat_Np(:,end-6);
        [TAC_min pos_min]=min(TAC);
        [C_BEC_min pos_C_BEC_min]=min(C_BEC);
        all_opt(jj,:)=mat_Np(pos_min,:)';
    end
end


%% Extra 1
for ii=1:length(all_opt(1,:))
    Summary(:,ii)=table(all_opt(:,ii));
end
Summary.Properties.VariableNames ={'d_o', 'N_p', 'layout', 'T_Na2', 'm_flow_Na', 'm_flow_MS', 'R', 'C_min', 'F', 'UA', 'en_eff', 'N_t', 'N_t_min', 'N_t_max', 'Tep', 'A_tot', 'A_st', 'D_s_out', 'D_s', 'L', 'N_baffles', 'N_baffles_min', 'N_baffles_max', 'l_b', 't_baffle', 't_tube', 't_shell', 'condition', 'U_calc', 'h_s', 'h_t', 'Dp_tube', 'Dp_shell', 'v_Na', 'v_max_MS', 'V_HX', 'm_HX', 'm_material_HX', 'F_ma', 'C_BEC', 'C_pump', 'penalty', 'rho_MS', 'rho_Na', 'ratio', 'TAC'};
writetable(Summary);

for ii=1:length(results_inter(1,:))
    Summary2(:,ii)=table(results_inter(:,ii));
end
Summary2.Properties.VariableNames ={'d_o', 'N_p', 'layout', 'T_Na2', 'm_flow_Na', 'm_flow_MS', 'R', 'C_min', 'F', 'UA', 'en_eff', 'N_t', 'N_t_min', 'N_t_max', 'Tep', 'A_tot', 'A_st', 'D_s_out', 'D_s', 'L', 'N_baffles', 'N_baffles_min', 'N_baffles_max', 'l_b', 't_baffle', 't_tube', 't_shell', 'condition', 'U_calc', 'h_s', 'h_t', 'Dp_tube', 'Dp_shell', 'v_Na', 'v_max_MS', 'V_HX', 'm_HX', 'm_material_HX', 'F_ma', 'C_BEC', 'C_pump', 'penalty', 'rho_MS', 'rho_Na', 'ratio', 'TAC'};

%% Extra 2
hh=1;
for jj=1:length(N_p)
    ind_N_p=find(results_inter(:,2)==N_p(jj));
    mat_N_p=results_inter(ind_N_p,:);
    for ii=1:length(d_o)
        ind=find(mat_N_p(:,1)==d_o(ii));
        L_origin=mat_N_p(ind,20);
        mat=mat_N_p(ind,:);
        [L_sort_origin, poss]=sort(L_origin); %Ordina
        mat_sort=mat(poss,:); %Portati i valori giusti
        L_sort_approx=round(L_sort_origin,0); %Approssima L solo per scegliere
        iter=1;
        while iter<length(L_sort_approx)
            ind_L=find(L_sort_approx(:)==L_sort_approx(iter)); %Trova valori uguali di L
            mat_sort_partial=mat_sort(ind_L,:);
            [TAC_min_partial pos_TAC_min_partial]=min(mat_sort_partial(:,end));
            res_inter_new(hh,:)=[mat_sort_partial(pos_TAC_min_partial,1:20),L_sort_approx(iter),mat_sort_partial(pos_TAC_min_partial,21:end)];
            hh=hh+1;
            iter=ind_L(end)+1;
        end
    end
end
for ii=1:length(res_inter_new(1,:))
    Summary3(:,ii)=table(res_inter_new(:,ii));
end
Summary3.Properties.VariableNames ={'d_o', 'N_p', 'layout', 'T_Na2', 'm_flow_Na', 'm_flow_MS', 'R', 'C_min', 'F', 'UA', 'en_eff', 'N_t', 'N_t_min', 'N_t_max', 'Tep', 'A_tot', 'A_st', 'D_s_out', 'D_s', 'L', 'L_round', 'N_baffles', 'N_baffles_min', 'N_baffles_max', 'l_b', 't_baffle', 't_tube', 't_shell', 'condition', 'U_calc', 'h_s', 'h_t', 'Dp_tube', 'Dp_shell', 'v_Na', 'v_max_MS', 'V_HX', 'm_HX', 'm_material_HX', 'F_ma', 'C_BEC', 'C_pump', 'penalty', 'rho_MS', 'rho_Na', 'ratio', 'TAC'};

[TAC_min pos]=min(Summary3.TAC);

ind_stampa=[1,2,3,12,21,46,47,41,16,42,30,38,37,35,36,22];
writetable(Summary3(:,ind_stampa),'OptimumConfigurations.txt');

