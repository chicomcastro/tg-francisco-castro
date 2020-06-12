function dz = dinamica2(t,z, parameters)
% function to be integrated
% The constants required in the calculation
Me= parameters.Me;               % Mass of earth [kg]
Ms= parameters.Ms;               % Mass of sun [kg]
G= parameters.G;                 % Gravitational Constant
w= parameters.w;                 % Angular velocity
de= parameters.de;               % Distance between earth and centre of mass
ds= parameters.ds;               % Distance between sun and centre of mass

dz = zeros(4,1);
% Equation for coordinates of sun and earth
xs= ds*cos(w*t);
ys= ds*sin(w*t);
xe=-de*cos(w*t);
ye=-de*sin(w*t);
% Differential equation for the coordinates of moon
a= sqrt(((xs - z(1))^2) + ((ys - z(3))^2));
b= sqrt(((xe - z(1))^2) + ((ye - z(3))^2));
Gs= G*Ms/power(a,3);
Ge= G*Me/power(b,3);
p = Gs*(xs-z(1)) + Ge*(xe-z(1));
q = Gs*(ys-z(3)) + Ge*(ye-z(3));
% % Writing Equation for gravitational acceleration on moon
% 
% a= sqrt(((ds*cos(w*t) - z(1))^2) + ((ds*sin(w*t) - z(3))^2));
% b= sqrt(((de*cos(w*t) + z(1))^2) + ((de*sin(w*t) + z(3))^2));
% p = (G*Ms*(ds*cos(w*t)-z(1)))/power(a,3) - (G*Me*(de*cos(w*t)+z(1)))/power(b,3);
% q = (G*Ms*(ds*sin(w*t)-z(3)))/power(a,3) - (G*Me*(de*sin(w*t)+z(3)))/power(b,3);
% Giving the values to function dz
dz(1) = z(2);
dz(2) = p;
dz(3) = z(4);
dz(4) = q;
