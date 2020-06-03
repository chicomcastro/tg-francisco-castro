# From https://ode-solver.readthedocs.io/en/master/double-pendulum-example.html
import numpy as np
import ode
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
# from IPython.display import HTML

def doublependulum(t,X):
    th1, th2, om1, om2 = X
    g = 9.81
    m1 = 1
    m2 = 1
    l1 = 1
    l2 = 1
    k1 = -g * ((2 * m1) + m2) * np.sin(th1)
    k2 = m2 * g * np.sin(th1 - (2 * th2))
    k3 = 2 * np.sin(th1 - th2) * m2
    k4 = ((om2**2) * l2) + ((om1**2) * l1 * np.cos(th1 - th2))
    k5 = m2 * np.cos((2 * th1) - (2 * th2))
    k6 = 2 * np.sin(th1 - th2)
    k7 = ((om1**2) * l1 * (m1 + m2))
    k8 = g * (m1 + m2) * np.cos(th1)
    k9 = (om2**2) * l2 * m2 * np.cos(th1 - th2)
    dX = np.array([
        om1,
        om2,
        (k1 - k2 - (k3 * k4)) / (l1 * ((2 * m1) + m2 - k5)),
        (k6 * (k7 + k8 + k9)) / (l2 * ((2 * m1) + m2 - k5))
    ])
    return dX

def angles2xy(th1, th2, om1, om2):
    l1 = 1
    l2 = 1
    x1 = l1 * np.sin(th1)
    y1 = -l1 * np.cos(th1)
    x2 = x1 + (l2 * np.sin(th2))
    y2 = y1 - (l2 * np.cos(th2))
    return x1, y1, x2, y2


et, ex = ode.euler(doublependulum, [np.pi/2,np.pi,0,0], [0,5], .025)
th1, th2, om1, om2 = ex
ex1, ey1, ex2, ey2 = zip(*[angles2xy(*xi) for xi in zip(th1, th2, om1, om2)])

bet, bex = ode.backwardeuler(doublependulum, [np.pi/2,np.pi,0,0], [0,5], .025)
bex1, bey1, bex2, bey2 = zip(*[angles2xy(*xi) for xi in zip(*bex)])

fig, ax = plt.subplots()
ax.plot(ex2, ey2, 'r', bex2, bey2, 'g')
ax.set_aspect('equal', 'datalim')
plt.show()
