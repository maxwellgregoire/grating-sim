#!/usr/bin/python

import numpy as np
import cmath
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import calc_sensitivity as cs
import nominal_param_values as npv

# plot S vs C3

f = np.linspace(0.1, 0.99, 20)
S = np.zeros(f.size)

for i in range(f.size):

    S[i] = cs.calc_sensitivity_vdw(npv.I_inc_nom, npv.l_nom, npv.L_nom, npv.v_nom, npv.C3_nom, npv.d_nom, f[i])[0]

plt.semilogy(f, S, 'bo')
plt.xlabel("g1 and g2 open fraction")
plt.ylabel("sensitivity (rad/s / sqrt(Hz))")
plt.show()

