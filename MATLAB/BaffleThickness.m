function t_baffle = BaffleThickness(D_s,l_b)

% TEMA Standards - Minimum Baffle Thickness
if l_b<=0.61  
  if D_s<=0.356   
    t_baffle=3.2e-3;
  elseif D_s>0.356 && D_s<=0.711  
    t_baffle=4.8e-3;
  elseif D_s>0.711 && D_s<=0.965  
    t_baffle=6.4e-3;
  elseif D_s>0.965 && D_s<=1.524  
    t_baffle=6.4e-3;
  elseif D_s>1.524  
    t_baffle=9.5e-3;
  end  
elseif l_b>0.61 && l_b<=0.914  
  if D_s<=0.356   
    t_baffle=4.8e-3;
  elseif D_s>0.356 && D_s<=0.711  
    t_baffle=6.4e-3;
  elseif D_s>0.711 && D_s<=0.965  
    t_baffle=7.5e-3;
  elseif D_s>0.965 && D_s<=1.524  
    t_baffle=9.5e-3;
  elseif D_s>1.524  
    t_baffle=12.7e-3;
  end  
elseif l_b>0.914 && l_b<=1.219  
  if D_s<=0.356   
    t_baffle=6.4e-3;
  elseif D_s>0.356 && D_s<=0.711  
    t_baffle=9.5e-3;
  elseif D_s>0.711 && D_s<=0.965  
    t_baffle=9.5e-3;
  elseif D_s>0.965 && D_s<=1.524  
    t_baffle=12.7e-3;
  elseif D_s>1.524  
    t_baffle=15.9e-3;
  end  
elseif l_b>1.219 && l_b<=1.524  
  if D_s<=0.356   
    t_baffle=9.5e-3;
  elseif D_s>0.356 && D_s<=0.711  
    t_baffle=9.5e-3;
  elseif D_s>0.711 && D_s<=0.965  
    t_baffle=12.7e-3;
  elseif D_s>0.965 && D_s<=1.524  
    t_baffle=15.9e-3;
  elseif D_s>1.524  
    t_baffle=19.1e-3;
  end  
elseif l_b>1.524  
  if D_s<=0.356   
    t_baffle=9.5e-3;
  elseif D_s>0.356 && D_s<=0.711  
    t_baffle=12.7e-3;
  elseif D_s>0.711 && D_s<=0.965  
    t_baffle=15.9e-3;
  elseif D_s>0.965 && D_s<=1.524  
    t_baffle=15.9e-3;
  elseif D_s>1.524  
    t_baffle=19.1e-3;
  end  
end  

end