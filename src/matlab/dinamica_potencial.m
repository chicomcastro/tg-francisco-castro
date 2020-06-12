function dY = dinamica_potencial(t,Y, parameters)
    
    mu = parameters.m_2;
    
    x = Y(1);
    y = Y(2);
    z = Y(3);
    xp = Y(4);
    yp = Y(5);
    zp = Y(6);
    
    r1 = sqrt((x+mu)^2+y^2+z^2);
    r2 = sqrt((x-(1-mu))^2+y^2+z^2);
    
    dUdx = x - (1-mu)*power(r1,-3)*(x+mu) - mu*power(r2,-3)*(x-(1-mu));
    dUdy = y - (1-mu)*power(r1,-3)*y - mu*power(r2,-3)*y;
    dUdz = (1-mu)*power(r1,-3)*z - mu*power(r2,-3)*z;
    
    Ypp = [
        2*yp + dUdx;
        -2*xp + dUdy;
        dUdz
    ];
    
    dY = [
        Y(4:6);
        Ypp
    ];