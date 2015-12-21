#!/usr/bin/python

import numpy as np
import cmath
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import calc_sensitivity as cs
import nominal_param_values as npv

# plot the optimal g1 and g2 open fractions as a function of velocity and grating period

v = np.linspace(100.0, 1200.0, 5)
d = np.linspace(50.0e-9, 200.0e-9, 5)
V, D = np.meshgrid(v,d*1.0e9)
f1 = np.zeros([v.size, d.size], dtype = np.float)
f2 = np.zeros([v.size, d.size], dtype = np.float)

for i in range(v.size):
    for j in range(d.size):

        result = cs.calc_signal_vdw(npv.I_inc_nom, d[j], npv.l_nom, v[i], npv.C3_nom)
        f1[i][j] = result[1]
        f2[i][j] = result[2]

fig = plt.figure(figsize=plt.figaspect(0.5))

ax = fig.add_subplot(1, 2, 1, projection='3d')
surf1 = ax.plot_surface(V, D, f1, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=True)
ax.set_xlabel("v (m/s)")
ax.set_ylabel("grating period (um)")
ax.set_zlabel("optimum g1 open fraction")
ax.set_zlim3d(-0.01, 1.01)

ax = fig.add_subplot(1, 2, 2, projection='3d')
surf2 = ax.plot_surface(V, D, f2, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=True)
ax.set_xlabel("v (m/s)")
ax.set_ylabel("grating period (um)")
ax.set_zlabel("optimum g2 open fraction")
ax.set_zlim3d(-0.01, 1.01)

plt.show()

