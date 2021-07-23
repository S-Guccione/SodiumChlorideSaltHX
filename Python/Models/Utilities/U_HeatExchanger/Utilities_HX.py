#Author: "Salvatre Guccione"
#First Version Released: "16 April 2021"

def TubeThickness(d_o):
    #TEMA Standards average thickness
    #t_tube=43.91*d_o**3 - 5.18*d_o**2 + 0.2235*d_o - 0.0007558
    #t_tube approximated with function
    #t_tube= 5.1873151374*d_o**3 - 0.93792823782*d_o**2 + 0.055292752767*d_o + 0.00014137538821
    #t_tube fixed
    #t_tube=1.2e-3
    #TEMA Standard minimum thickness + NREL info maximum thickness 1.2mm
    if d_o<=7.9e-3:
        t_tube=0.5e-3
    elif (d_o>7.9e-3 and d_o<=11.1e-3):
        t_tube=0.6e-3
    elif (d_o>11.1e-3 and d_o<=14.3e-3):
        t_tube=0.7e-3
    elif (d_o>14.3e-3 and d_o<=34.9e-3):
        t_tube=0.9e-3
    elif d_o>34.9e-3:
        t_tube=1.2e-3
    return t_tube

def Inconel800H_BaseProperties(T):
    k_wall=-0.000025*T**2 + 0.05816*T - 8.52 #Interpolated from datasheet table in the temperature range 600-800 °C
    rho_wall=8030
    #YY=-0.06952*T+219.7
    #YY(unit="GPa") - Young's Modulus
    return k_wall, rho_wall


def Haynes230_BaseProperties(T):
    k_wall=0.01996*T + 2.981 #Interpolated from datasheet table in the temperature range 20-1000 °C
    rho_wall=8970
    # YY=0.0000000000005663*T**5 - 0.000000002271*T**4 + 0.000003446*T**3 - 0.002469*T**2 + 0.768*T + 124.8 #Interpolated from datasheet table in the temperature range 20-1000 °C
    # cp_H230=0.00000000000001913*T**6 - 0.0000000000899*T**5 + 0.000000167*T**4 - 0.0001557*T**3 + 0.07654*T**2 - 18.57*T + 2137 #Interpolated from datasheet table in the temperature range 20-1000 °C
    #YY(unit="GPa") - Young's Modulus
    #cp_H230 - Specific Heat
    return k_wall, rho_wall

def ShellThickness(D_s):
    #TEMA Standards - Minimum Shell Thickness - Approximated Function
    #  if D_s<=0.15: 
    #    t_shell=3.2e-3
    #  elif D_s>0.15 and D_s<=2.03:
    #    t_shell=0.006560749*D_s^0.3821246
    #  else
    #    t_shell=9.5e-3
    #  end if
    #TEMA Standards - Minimum Shell Thickness
    if (D_s<=0.15):
        t_shell=3.2e-3
    elif (D_s>0.15 and D_s<=0.3):
        t_shell=3.2e-3
    elif (D_s>0.3 and D_s<=0.58):
        t_shell=3.2e-3
    elif (D_s>0.58 and D_s<=0.74):
        t_shell=4.8e-3
    elif (D_s>0.74 and D_s<=0.99):
        t_shell=6.4e-3
    elif (D_s>0.99 and D_s<=1.52):
        t_shell=6.4e-3
    elif (D_s>1.52 and D_s<=2.03):
        t_shell=7.9e-3
    elif (D_s>2.03):
        t_shell=9.5e-3
    return t_shell

def ShellDiameter(d_o, N_t, layout, N_p):
    #Shell Diameter
    if layout==1:
        if N_p==1:
            KK1=0.215
            nn1=2.207
        elif N_p==2:
            KK1=0.156
            nn1=2.291
        elif N_p==4:
            KK1=0.158
            nn1=2.263
        elif N_p==6:
            KK1=0.0402
            nn1=2.617
        elif N_p==8:
            KK1=0.0331
            nn1=2.643
    else:
        if N_p==1:
            KK1=0.319
            nn1=2.142
        elif N_p==2:
            KK1=0.249
            nn1=2.207
        elif N_p==4:
            KK1=0.175
            nn1=2.285
        elif N_p==6:
            KK1=0.0743
            nn1=2.499
        elif N_p==8:
            KK1=0.0365
            nn1=2.675
    D_b=(N_t/KK1)**(1/nn1)*d_o
    L_bb=(12+5*(D_b+d_o))/995
    D_s=L_bb+D_b+d_o
    t_shell=ShellThickness(D_s)
    D_s_out=D_s+2*t_shell
    return L_bb, D_b, D_s, D_s_out

def BaffleThickness(D_s, l_b):
    #TEMA Standards - Minimum Baffle Thickness
    if (l_b<=0.61):
        if (D_s<=0.356): 
            t_baffle=3.2e-3
        elif (D_s>0.356 and D_s<=0.711):
            t_baffle=4.8e-3
        elif (D_s>0.711 and D_s<=0.965):
            t_baffle=6.4e-3
        elif (D_s>0.965 and D_s<=1.524):
            t_baffle=6.4e-3
        elif (D_s>1.524):
            t_baffle=9.5e-3
    elif (l_b>0.61 and l_b<=0.914):
        if (D_s<=0.356): 
            t_baffle=4.8e-3
        elif (D_s>0.356 and D_s<=0.711):
            t_baffle=6.4e-3
        elif (D_s>0.711 and D_s<=0.965):
            t_baffle=7.5e-3
        elif (D_s>0.965 and D_s<=1.524):
            t_baffle=9.5e-3
        elif (D_s>1.524):
            t_baffle=12.7e-3
    elif (l_b>0.914 and l_b<=1.219):
        if (D_s<=0.356): 
            t_baffle=6.4e-3
        elif (D_s>0.356 and D_s<=0.711):
            t_baffle=9.5e-3
        elif (D_s>0.711 and D_s<=0.965):
            t_baffle=9.5e-3
        elif (D_s>0.965 and D_s<=1.524):
            t_baffle=12.7e-3
        elif (D_s>1.524):
            t_baffle=15.9e-3
    elif (l_b>1.219 and l_b<=1.524):
        if (D_s<=0.356): 
            t_baffle=9.5e-3
        elif (D_s>0.356 and D_s<=0.711):
            t_baffle=9.5e-3
        elif (D_s>0.711 and D_s<=0.965):
            t_baffle=12.7e-3
        elif (D_s>0.965 and D_s<=1.524):
            t_baffle=15.9e-3
        elif (D_s>1.524):
            t_baffle=19.1e-3
    elif (l_b>1.524):
        if (D_s<=0.356): 
            t_baffle=9.5e-3
        elif (D_s>0.356 and D_s<=0.711):
            t_baffle=12.7e-3
        elif (D_s>0.711 and D_s<=0.965):
            t_baffle=15.9e-3
        elif (D_s>0.965 and D_s<=1.524):
            t_baffle=15.9e-3
        elif (D_s>1.524):
            t_baffle=19.1e-3
    
    #Fixed Baffle Thickness
    #t_baffle=19.1e-3

    #TEMA Standards - Minimum Baffle Thickness - Approximated Function
    #if l_b<=0.61:
    #  if D_s<=0.356: 
    #    t_baffle=3.2e-3
    #  elif D_s>0.356 and D_s<=1.254:
    #    t_baffle=0.004332389*log(D_s) + 0.007674598
    #  else
    #    t_baffle=9.5e-3
    #  end if
    #elif l_b>0.61 and l_b<=0.914:
    #  if D_s<=0.356: 
    #    t_baffle=0.0048
    #  elif D_s>0.356 and D_s<=1.254:
    #    t_baffle=0.005432678*log(D_s) + 0.010411
    #  else
    #    t_baffle=0.0127
    #  end if
    #elif l_b>0.914 and l_b<=1.219:
    #  if D_s<=0.356: 
    #    t_baffle=0.0066
    #  elif D_s>0.356 and D_s<=1.254:
    #    t_baffle=0.0064457*log(D_s) + 0.01326375
    #  else
    #    t_baffle=0.0160
    #  end if
    #elif l_b>1.219 and l_b<=1.524:
    #  if D_s<=0.356: 
    #    t_baffle=0.0095
    #  elif D_s>0.356 and D_s<=1.254:
    #    t_baffle=0.006601736*log(D_s) + 0.01631843
    #  else
    #    t_baffle=0.0191
    #  end if
    #elif l_b>1.524:
    #  if D_s<=0.356: 
    #    t_baffle=0.0098
    #  elif D_s>0.356 and D_s<=1.254:
    #    t_baffle=0.006591127*log(D_s) + 0.01663412
    #  else
    #    t_baffle=0.0194
    #  end if
    #end if
    return t_baffle