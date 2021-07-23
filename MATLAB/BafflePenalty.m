function F=BafflePenalty(N_b,NTU)

    a= 2.474007E-01*N_b^(-1.701419E+00);
    c=-0.00007108634*N_b^2 + 0.003993761*N_b + 0.2251199;
    if N_b>3
        b=1.928;
    elseif N_b==1
        b=1.875;
    elseif N_b==2
        b=1.924;
    else
        b=1.921;
    end
    F=1/(1+a*NTU^b)^c;

end