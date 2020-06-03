# Initial tryings for solving a R3BP

# -------------
# Normalizações
# -------------
# A unidade de comprimento será a distância entre os corpos massivos
d_unit = 3.844e8                    # [m]

# A unidade de massa será a soma das massas dos corpos massivos
m_1 = 5.972e24                      # Massa Terra [kg]
m_2 = 7.35e22                       # Massa Lua [kg]
m_unit = m_1 + m_2                  # [kg]

# A unidade de tempo será tal que a velocidade angular de rotação dos massivos
# em torno de seu centro de massa seja igual a 1
# Lembrando que $/omega = /frac{2/pi}{T} = /frac{G(m_1+m_2)}{d^3}$
G = 6.67408e-11                     # [m3 kg-1 s-2]
t_unit = G * m_unit / d_unit        # [s]

# E isso nos dá
G_unit = 1