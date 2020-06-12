function [C, erro] = constante_jacobi(Y, parameters)

    mu = parameters.m_2;
    
    C0 = calculate_C(Y(1,:), mu);
    
    for i = 1:size(Y,1)
        C(i) = calculate_C(Y(i,:), mu);
        erro(i) = C0 - C(i);
    end
    
    C = C';
    erro = erro';
end

function C_i = calculate_C(Y, mu)
    x = Y(1);
    y = Y(2);
    z = Y(3);
    xp = Y(4);
    yp = Y(5);
    zp = Y(6);

    r1 = sqrt((x+mu)^2+y^2+z^2);
    r2 = sqrt((x-(1-mu))^2+y^2+z^2);
    
    C_i = x^2 + y^2 + 2*(1-mu)/r1 + 2*mu/r2 - (xp^2 + yp^2 + zp^2);    
end
    