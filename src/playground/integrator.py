# First
print('First')
# From https://stackoverflow.com/questions/48428140/imitate-ode45-function-from-matlab-in-python
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
import time

show_plots = False

def vdp1(t, y):
    return np.array([y[1], (1 - y[0]**2)*y[1] - y[0]])


t0, t1 = 0, 20                # start and end
t = np.linspace(t0, t1, 100)  # the points of evaluation of solution
y0 = [2, 0]                   # initial value
y = np.zeros((len(t), len(y0)))   # array for solution
y[0, :] = y0
tic = time.time()
r = integrate.ode(vdp1).set_integrator("dopri5")  # choice of method
r.set_initial_value(y0, t0)   # initial values

for i in range(1, t.size):
   y[i, :] = r.integrate(t[i]) # get one more value, add it to the array
   if not r.successful():
       raise RuntimeError("Could not integrate")

toc = time.time()
print(toc-tic)
y1 = y
plt.plot(t, y)
if show_plots:
    plt.show()

# ------
# Second
print('Second')
tic = time.time()
results = integrate.solve_ivp(vdp1, [t0, t1], y0, dense_output=True)
toc = time.time()
print(toc-tic)
Y = results
z = Y.sol(t)
y2 = z.T[:, 0]
plt.plot(t, z.T)
if show_plots:
    plt.show()

# Third
print('Third')
# From https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html
# scipy.integrate.solve_ivp(fun, t_span, y0, method='RK45', t_eval=None, dense_output=False, events=None, vectorized=False, args=None, **options)

from scipy.integrate import solve_ivp

# Default uses RK45
tic = time.time()
sol = solve_ivp(vdp1, [t0, t1], y0, dense_output=True)
toc = time.time()
print(toc-tic)
z = sol.sol(t)
y3 = z.T[:, 0]
plt.plot(t, z.T)
plt.xlabel('t')
plt.legend(['x', 'y'], shadow=True)
if show_plots:
    plt.show()


# Fourth
print('Fourth')
tic = time.time()
sol = solve_ivp(vdp1, [t0, t1], y0, method='DOP853', dense_output=True)
toc = time.time()
print(toc-tic)
z = sol.sol(t)
y4 = z.T[:, 0]
plt.plot(t, z.T)
plt.xlabel('t')
plt.legend(['x', 'y'], shadow=True)
if show_plots:
    plt.show()


# Fifth
print('Fifth')
from scipy.integrate import odeint
tic = time.time()
sol = odeint(vdp1, y0, t, tfirst=True)
toc = time.time()
print(toc-tic)
y5 = sol[:, 0]
plt.plot(t, sol[:, 0], 'b', label='x(t)')
plt.plot(t, sol[:, 1], 'g', label='y(t)')
plt.legend(loc='best')
plt.xlabel('t')
plt.grid()
if show_plots:
    plt.show()


def max_diff(a,b):
    a2 = np.array(a)
    b2 = np.array(b)
    c = a2 - b2
    d = [np.abs(elem) for elem in c]
    return np.max(d)

print('Diffs')
print(max_diff(y1,y5))
print(max_diff(y2,y5))
print(max_diff(y3,y5))
print(max_diff(y4,y5))
print(max_diff(y5,y5))