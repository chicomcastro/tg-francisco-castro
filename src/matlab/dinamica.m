function dY = dinamica(t,Y, parameters)
    
    w = parameters.w;
    G = parameters.G;
    m_1 = parameters.m_1;
    m_2 = parameters.m_2;
    P_1 = parameters.P_1;
    P_2 = parameters.P_2;
    P_3 = Y(1:3);     % P_3 = [x_3, y_3, z_3]'
    Pp_3 = Y(4:6);
    P_31 = P_3 - P_1;
    P_32 = P_3 - P_2;
    Ypp = w^2*eye(3)*P_3 + ...
        2*[ 0 w 0;
           -w 0 0;
            0 0 0 ]*Pp_3 - ...
        G*(m_1/norm(P_31)^3*P_31 +...
        m_2/norm(P_32)^3*P_32);
    
    dY = [
        Pp_3;
        Ypp
    ];