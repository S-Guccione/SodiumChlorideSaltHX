function [L_bb, D_b, D_s] = ShellDiameter(d_o, N_t, layout, N_p)

%Shell Diameter
  if layout==1 
    if N_p==1 
      KK1=0.215;
      nn1=2.207;
    elseif N_p==2 
      KK1=0.156;
      nn1=2.291;
    elseif N_p==4 
      KK1=0.158;
      nn1=2.263;
    elseif N_p==6 
      KK1=0.0402;
      nn1=2.617;
    elseif N_p==8 
      KK1=0.0331;
      nn1=2.643;
    end
  elseif layout==2
    if N_p==1 
      KK1=0.319;
      nn1=2.142;
    elseif N_p==2 
      KK1=0.249;
      nn1=2.207;
    elseif N_p==4 
      KK1=0.175;
      nn1=2.285;
    elseif N_p==6 
      KK1=0.0743;
      nn1=2.499;
    elseif N_p==8 
      KK1=0.0365;
      nn1=2.675;
    end
  end
  D_b=(N_t/KK1)^(1/nn1)*d_o;
  L_bb=(12+5*(D_b+d_o))/995;
  D_s=L_bb+D_b+d_o;

end