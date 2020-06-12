function [m_unit, d_unit, t_unit, mu] = setup(m_1, m_2, d_12)

    % A unidade de comprimento será a distância entre os corpos massivos
    d_unit = d_12;                                            % [m]

    % A unidade de massa será a soma das massas dos corpos massivos
    m_unit = m_1 + m_2;                                       % [kg]

    % A unidade de tempo será tal que a velocidade angular de rotação dos massivos
    % em torno de seu centro de massa seja igual a 1
    % Lembrando que $/omega = /frac{2/pi}{T} = /frac{G(m_1+m_2)}{d^3}$
    G = 6.67408e-11;                                          % [m3 kg-1 s-2]
    t_unit = 1/power((G * m_unit / power(d_unit, 3)), 0.5);   % [s]

    % E isso nos dá
    G_unit = 1;

    % Temos, além disso
    mu = m_2 / (m_1 + m_2);

    % Que nos dá, no sistema normalizado
    m_1_norm = mu;
    m_2_norm = 1 - mu;