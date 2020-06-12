# -*- coding: utf-8 -*-
# Initial tryings for solving a R3BP

# -------------
# Normalizações
# -------------
import numpy as np


def setup(m_1, m_2, d_12):

    # A unidade de comprimento será a distância entre os corpos massivos
    d_unit = d_12                                       # [m]

    # A unidade de massa será a soma das massas dos corpos massivos
    m_unit = m_1 + m_2                                  # [kg]

    # A unidade de tempo será tal que a velocidade angular de rotação dos massivos
    # em torno de seu centro de massa seja igual a 1
    # Lembrando que $/omega = /frac{2/pi}{T} = /frac{G(m_1+m_2)}{d^3}$
    G = 6.67408e-11                                     # [m3 kg-1 s-2]
    t_unit = 1/pow((G * m_unit / pow(d_unit, 3)), 0.5)   # [s]

    # E isso nos dá
    G_unit = 1

    # Temos, além disso
    mu = m_2 / (m_1 + m_2)

    # Que nos dá, no sistema normalizado
    m_1_norm = mu
    m_2_norm = 1 - mu

    return m_unit, t_unit, mu


# --------------------
# Equação de movimento
# --------------------

# y = (x, y)^T
# mu = m1 / (m1 + m2)

# Escalar
# def r1(Y, mu):
#     x = Y[0]
#     y = Y[1]
#     return pow(
#         pow(x + mu, 2) +
#         pow(y, 2),
#         0.5
#     )

# def r2(Y, mu):
#     x = Y[0]
#     y = Y[1]
#     return pow(
#         pow(x - (1 - mu), 2) +
#         pow(y, 2),
#         0.5
#     )

# def U(Y, mu):
#     x = Y[0]
#     y = Y[1]
#     return (
#         1 / 2 * (pow(x, 2) + pow(y, 2)) +
#         (1 - mu) / r1(Y, mu) +
#         mu / r2(Y, mu)
#     )

def dUdx(Y, mu):
    x = Y[0]
    y = Y[1]
    return (
        x -
        (mu*(mu + x - 1)) / pow(pow(mu + x - 1, 2) + pow(y, 2), 1.5) +
        ((mu + x)*(mu - 1)) / pow(pow(mu + x, 2) + pow(y, 2), 1.5)
    )

def dUdy(Y, mu):
    x = Y[0]
    y = Y[1]
    return (
        y -
        (mu*y) / pow(pow(mu + x - 1, 2) + pow(y, 2), 1.5) +
        (y*(mu - 1)) / pow(pow(mu + x, 2) + pow(y, 2), 1.5)
    )

# Vetor
def K(Y, mu):
    return np.array([
        dUdx(Y, mu),
        dUdy(Y, mu)
    ])

# vetor de estado
# Y = [x, y, xp, yp]
# f(t, Y) = Yp = [xp, yp, xpp, ypp]
def f(t, Y, mu):
    x = Y[0]
    y = Y[1]
    xp = Y[2]
    yp = Y[3]
    Y1 = np.array([x, y])
    Y2 = np.array([xp, yp])
    return np.concatenate((
        Y2,
        -2 * np.array([[0, -1], [1, 0]]).dot(Y2) + K(Y1, mu)
    ))
