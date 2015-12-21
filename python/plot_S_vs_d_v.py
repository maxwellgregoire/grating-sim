#!/usr/bin/python

import numpy as np
import cmath
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import calc_sensitivity as cs
import nominal_param_values as npv

# plot S vs d and v

v = np.linspace(20.0, 200.0, 10)
d = np.linspace(50.0e-9, 200.0e-9, 10)
V, D = np.meshgrid(v,d*1.0e9)
S = np.zeros([v.size, d.size], dtype = np.float)
f1 = np.zeros([v.size, d.size], dtype = np.float)
f2 = np.zeros([v.size, d.size], dtype = np.float)

for i in range(v.size):
    for j in range(d.size):
        print i,j

        result = cs.calc_sensitivity_vdw(npv.I_inc_nom, npv.l_nom, npv.L_nom, v[i], npv.C3_nom, fixed_d = d[j])
        S[i][j] = result[0]
        f1[i][j] = result[2]
        f2[i][j] = result[3]

fig = plt.figure(figsize=plt.figaspect(0.3))

ax = fig.add_subplot(1, 3, 1, projection='3d')
surf1 = ax.plot_surface(V, D, S, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=True)
ax.set_xlabel("v (m/s)")
ax.set_ylabel("grating period (um)")
ax.set_zlabel("sensitivity (rad/s / sqrt(Hz))")

print "fig 2"
ax = fig.add_subplot(1, 3, 2, projection='3d')
surf2 = ax.plot_surface(V, D, f1, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=True)
ax.set_xlabel("v (m/s)")
ax.set_ylabel("grating period (um)")
ax.set_zlabel("optimum g1 open fraction")

print "fig 3"
ax = fig.add_subplot(1, 3, 3, projection='3d')
surf3 = ax.plot_surface(V, D, f2, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=True)
ax.set_xlabel("v (m/s)")
ax.set_ylabel("grating period (um)")
ax.set_zlabel("optimum g2 open fraction")

print "show"
plt.show()
