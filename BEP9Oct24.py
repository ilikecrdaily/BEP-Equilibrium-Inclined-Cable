# -*- coding: utf-8 -*-
"""
Created on Sat Sep  7 12:03:30 2024

@author: ilike
"""

from math import sin, cos, sqrt, log, pi
import numpy as np
import matplotlib.pyplot as plt

"""Physical constants"""
g = 9.81

"""Physical parameters"""
M = 486.2*1
A = 1.11e-3
E = 21e10
rho = 48.62/A
X_1 = 4
L = 200
phi = pi/6
T_0 = 4.64e7

omega_0 = T_0/(A*E)
zeta = M*g/(A*E)
epsilon = rho * g * L / E
c = zeta/epsilon
X_1_star = X_1/L

"""Other params"""
no_of_points = 10001

"""Our model (Cable-TMD model) used to obtain the equilibrium"""
def get_v(x):
    X = 200 * (1-x)
    if X < X_1:
        return -1/L*(- 3*A*cos(phi)*X*g*rho*(A*((-1/36 + (X_1_star**2 - 4/3*X_1_star + 1/3)*c**2 +
                                            (1/2*X_1_star**2 - 2/3*X_1_star + 1/6)*c)*L**2 + 
                                           ((-1/4 + c*(X_1_star - 1))*X*L)/3 + X**2/9)*g*rho*sin(phi) 
                                        - (((-1/2 + c*(X_1_star - 1))*L + X/2)*(T_0/(A*E) + 1)*T_0)/3))/T_0**2
    else:
        return -1/L*(3*A*cos(phi)*rho*g*(L - X)*(A*(((X_1_star - 1/3)*c + X_1_star/2 - 1/3)*c*X_1_star*L**2 + 
                                               ((1/12 + (X_1_star + 1/2)*c)*X*L)/3 + X**2/9)*g*rho*sin(phi) - 
                                            ((T_0/(A*E) + 1)*(L*X_1_star*c + X/2)*T_0)/3))/T_0**2
    
    
"""The parabola used in the main paper to estimate the equilibrium"""
def parab(x):
    d_c = rho * A * g * L * cos(phi) / (8*T_0)
    return 4 * d_c*x*(1-x)


"""The approximation for the cable without mass derived in caswita"""
def caswform(x):
    v = rho * g * A * L * cos(phi)/(2*T_0)*x*(1-x)* \
        (1+omega_0-rho*g*A*L*sin(phi)/(6*T_0)*(5-4*x))
    return v
    

x = np.linspace(0, 1, no_of_points)
V = np.array([get_v(i) for i in x])

plt.figure(dpi=666)
plt.plot(x, np.zeros(len(x)), color='black')
plt.plot(x, parab(x), '-', linewidth=1, color='green', label='Su et al.')
casw = caswform(x)
plt.plot(x, casw, '-', linewidth=3.0, color='red', label='Caswita')
plt.plot(x, V, '-', color='blue', linewidth=1.0, alpha=0.75, label='Cable-TMD model')
plt.legend(loc='lower left', bbox_to_anchor=(0,0))
plt.xlabel('$x$')
plt.ylabel('$y_c(x)$')
plt.show()

plt.figure(dpi=666)
plt.plot(x, np.zeros(len(x)), color='black')
plt.plot(x[1:], np.diff(parab(x)*no_of_points), '-', linewidth=1, color='green', label='Su et al.')
casw = caswform(x)
plt.plot(x[1:], np.diff(casw)*no_of_points, '-', linewidth=3.0, color='red', label='Caswita')
plt.plot(x[1:], np.diff(V)*no_of_points, '-', color='blue', alpha=0.75, label='Cable-TMD model')
plt.legend(loc='lower left', bbox_to_anchor=(1,0))
plt.xlabel('$x$')
plt.ylabel('$y\'_c(x)$')