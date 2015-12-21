#!/usr/bin/python

import numpy as np
import cmath
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import calc_sensitivity as cs
import nominal_param_values as npv

# plot S vs C3

C3 = np.linspace(0.0, 10.0e-49, 30)
S_max = np.zeros(C3.size)
S_5 = np.zeros(C3.size)
S_kd = np.zeros(C3.size)

s_result_kd = cs.calc_sensitivity_kd(npv.I_inc_nom, npv.L_nom, npv.v_nom, npv.d_nom)

for i in range(C3.size):

    S_max[i] = cs.calc_sensitivity_vdw(npv.I_inc_nom, npv.l_nom, npv.L_nom, npv.v_nom, C3[i], npv.d_nom, cs.f_max)[0]
    S_5[i] = cs.calc_sensitivity_vdw(npv.I_inc_nom, npv.l_nom, npv.L_nom, npv.v_nom, C3[i], npv.d_nom, 0.5)[0]
    S_kd[i] = s_result_kd


C3_plot = C3*10.0e49

l1, l2, l3 = plt.plot(C3_plot, S_max, 'bo', C3_plot, S_5, 'ro', C3_plot, S_kd, 'k-')
plt.xlabel("C3 (10^-49 J*m^3)")
plt.ylabel("sensitivity (rad/s / sqrt(Hz))")
l1_str = "f = " + str(cs.f_max)
plt.legend((l1, l2, l3), (l1_str, "f = 0.5", "K-D gratings"), loc='upper right')
plt.show()

