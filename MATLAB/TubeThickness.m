function t_tube = TubeThickness(d_o)
  
%TEMA Standard minimum thickness + NREL info maximum thickness 1.2mm
  if d_o<=7.9e-3   
    t_tube=0.5e-3;
  elseif d_o>7.9e-3 && d_o<=11.1e-3  
    t_tube=0.6e-3;
  elseif d_o>11.1e-3 && d_o<=14.3e-3  
    t_tube=0.7e-3;
  elseif d_o>14.3e-3 && d_o<=34.9e-3  
    t_tube=0.9e-3;
  elseif d_o>34.9e-3  
    t_tube=1.2e-3;
  end  

end