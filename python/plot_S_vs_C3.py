#!/usr/bin/python

import numpy as np
import cmath
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import calc_sensitivity as cs
import nominal_param_values as npv

hartree = 4.36e-18 # J
bohr_radius = 5.29e-11 # m

# plot S vs C3

C3 = np.linspace(0.1, 1.5, 20)
S_max = np.zeros(C3.size)
S_kd = np.zeros(C3.size)

s_result_kd = cs.calc_sensitivity_kd(npv.I_inc_nom, npv.L_nom, npv.v_nom, npv.d_nom)

for i in range(C3.size):

    result = cs.calc_sensitivity_vdw(npv.I_inc_nom, npv.l_nom, npv.L_nom, npv.v_nom, C3[i]*hartree*bohr_radius**3, npv.d_nom)
    S_max[i] = result[0]
    S_kd[i] = s_result_kd
    print result


l1, l2 = plt.plot(C3, S_max, 'b-', C3, S_kd, 'r-')
plt.xlabel("C3 (atomic units)")
plt.ylabel("sensitivity (rad/s / sqrt(Hz))")
l1_str = "material gratings, optimal open fraction = " + str(cs.f_max(npv.d_nom))
plt.legend((l1, l2), ("material gratings, optimal open fraction", "Kapitza-Dirac gratings (open fraction = 0.5)"), loc='upper right')
plt.show()

