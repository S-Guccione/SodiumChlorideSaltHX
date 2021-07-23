function t_shell = ShellThickness(D_s)

% TEMA Standards - Minimum Shell Thickness
  if D_s<=0.15   
    t_shell=3.2e-3;
  elseif D_s>0.15 && D_s<=0.3  
    t_shell=3.2e-3;
  elseif D_s>0.3 && D_s<=0.58  
    t_shell=3.2e-3;
  elseif D_s>0.58 && D_s<=0.74  
    t_shell=4.8e-3;
  elseif D_s>0.74 && D_s<=0.99  
    t_shell=6.4e-3;
  elseif D_s>0.99 && D_s<=1.52  
    t_shell=6.4e-3;
  elseif D_s>1.52 && D_s<=2.03  
    t_shell=7.9e-3;
  elseif D_s>2.03  
    t_shell=9.5e-3;
  end  

end