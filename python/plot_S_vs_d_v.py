#!/usr/bin/python

import numpy as np
import cmath
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import calc_sensitivity as cs
import nominal_param_values as npv

# plot S vs d and v

v = np.linspace(20.0, 100.0, 5)
d = np.linspace(80.0e-9, 400.0e-9, 15)
V, D = np.meshgrid(v,d*1.0e9)
S = np.zeros([d.size, v.size], dtype = np.float)
S_kd = np.zeros([d.size, v.size], dtype = np.float)
f1 = np.zeros([d.size, v.size], dtype = np.float)
f2 = np.zeros([d.size, v.size], dtype = np.float)

for i in range(v.size):
    for j in range(d.size):
        print i,j

        result = cs.calc_sensitivity_vdw(npv.I_inc_nom, npv.l_nom, npv.L_nom, v[i], npv.C3_nom, fixed_d = d[j])
        S[j][i] = result[0]
        f1[j][i] = result[2]
        f2[j][i] = result[3]
        S_kd[j][i] = cs.calc_sensitivity_kd(npv.I_inc_nom, npv.L_nom, v[i], d[j])

fig = plt.figure(figsize=plt.figaspect(0.8))

ax = fig.add_subplot(2, 2, 1, projection='3d')
surf1 = ax.plot_surface(V, D, S, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=True)
ax.set_xlabel("v (m/s)")
ax.set_ylabel("grating period (um)")
ax.set_zlabel("S (rad/s / sqrt(Hz))")
ax.set_title("a. material gratings")

ax = fig.add_subplot(2, 2, 2, projection='3d')
surf1 = ax.plot_surface(V, D, S_kd, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=True)
ax.set_xlabel("v (m/s)")
ax.set_ylabel("grating period (um)")
ax.set_zlabel("S (rad/s / sqrt(Hz))")
ax.set_title("b. Kapitza-Dirac gratings")

ax = fig.add_subplot(2, 2, 3, projection='3d')
surf2 = ax.plot_surface(V, D, f1, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=True)
ax.set_xlabel("v (m/s)")
ax.set_ylabel("grating period (um)")
ax.set_title("c. best g1 open fraction for material gratings")

ax = fig.add_subplot(2, 2, 4, projection='3d')
surf3 = ax.plot_surface(V, D, f2, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=True)
ax.set_xlabel("v (m/s)")
ax.set_ylabel("grating period (um)")
ax.set_title("d. best g2 open fraction for material gratings")

plt.show()
