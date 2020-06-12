clear all
close all

%% The constants required in the calculation
m_terra= 5.972E24;                % Mass of earth [kg]
m_sol= 1.988E30;                % Mass of sun [kg]
m_lua= 7.3459E22;               % Mass of Moon [kg]
d_terra_sol = 1.495E11;                 % Orbital radius of earth [m]
d_terra_lua = 3.85E8;                   % Orbital radius of moon
G = 6.67E-11;                 % Gravitational Constant
T_terra = 365.2563*24*3600;        % orbital time of earth [s]

%% Normalization
[m_unit, d_unit, t_unit, mu] = setup(m_sol, m_terra, d_terra_sol);
%[m_unit, d_unit, t_unit, mu] = setup(m_terra, m_lua, d_terra_lua);

%% Defining tspan, initial condition.
tspan = linspace(0, 6*2*pi, 1E3);
parameters.w = 1;
parameters.G = 1;
parameters.m_1 = 1 - mu;
parameters.m_2 = mu;
parameters.P_1 = [-mu; 0; 0];
parameters.P_2 = [1 - mu; 0; 0];

P_3_0 = [0.994; 0; 0];
Pp_3_0 = [0; -0.5; 0];
Y0 = [P_3_0; Pp_3_0];

ic = 1;
opts = odeset('RelTol',1e-14,'AbsTol',1e-16);
[t,Y] = ode45(@(t,Y)dinamica_potencial(t,Y,parameters), tspan, Y0, opts);

%% Plotting
fig = figure;
ax = axes;
plot(parameters.P_1(1),parameters.P_1(2),'yo');
hold on;
plot(parameters.P_2(1),parameters.P_2(2),'bo');
grid on;
plot(Y(:,1),Y(:,2),'b-');
legend('Sol', 'Terra', 'Lua');
axis equal;
hold off

%% Analise de erro
[C, erro] = constante_jacobi(Y, parameters);

figure;
title('Constante de Jacobi')
plot(tspan, C);

figure;
title('Erro absoluto')
plot(tspan, erro);


%%


















