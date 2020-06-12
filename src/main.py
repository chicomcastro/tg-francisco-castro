import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import p3cr
import numpy as np

# Normalizações
d_terra_lua = 3.844e8                    # [m]
m_terra = 5.972e24                       # Massa Terra [kg]
m_lua = 7.35e22                          # Massa Lua [kg]

m_unit, t_unit, mu = p3cr.setup(m_terra, m_lua, d_terra_lua)
print(m_unit, t_unit, mu)
# Integração
t0, t1 = 0, 100                          # start and end
# Y = [x, y, xp, yp]
y0 = [0.994, 0, 0, -2.0317]                       # initial value

# Default uses RK4
sol = solve_ivp(p3cr.f, [t0, t1], y0, dense_output=True, args=([mu]))
t = np.linspace(t0, t1, 100)            # the points of evaluation of solution
z = sol.sol(t)

# plt.plot(t, z.T)
# plt.xlabel('t')
# plt.legend(['x', 'y', 'xp', 'yp'], shadow=True)
# plt.show()

x = z.T[:, 0]
y = z.T[:, 1]
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')
plt.legend(['Trajetoria'], shadow=True)
plt.show()
