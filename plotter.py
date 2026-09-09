# Essentially all of this was taken from scipy documentation.

import numpy as np

def SIR_saturated(y,t,b,a,k):
    S,I,R = y
    S_dot = -b*S*I/(k+I)
    I_dot  = b*S*I/(k+I) - a*I
    R_dot  = a*I
    return [S_dot,I_dot,R_dot]

b = 2
a = 0.2
k = 1.5

S_0  = 0.7
I_0  = 0.3
R_0  = 0

y0 = [S_0,I_0,R_0]

t = np.linspace(0, 20, 101)

from scipy.integrate import odeint

sol = odeint(SIR_saturated, y0, t, args=(b, a, k))

import matplotlib.pyplot as plt
plt.plot(t, sol[:, 0], label='S(t)')
plt.plot(t, sol[:, 1], label='I(t)')
plt.plot(t, sol[:, 2], label='R(t)')
plt.legend(loc='best')
plt.xlabel('t')
plt.grid()
plt.show()
