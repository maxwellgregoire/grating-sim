#!/usr/bin/python

import time
import numpy as np
import matplotlib.pyplot as plt
import calc_sensitivity as cs
import nominal_param_values as npv


v = 50.0 # m/s

d = np.linspace(50.0e-9,250.0e-9,8)
d_plot = d*1.0e9
S_vdw = np.zeros(d.size, dtype=np.float)
for i in range(d.size):
    print d[i]

    S_vdw[i] = cs.calc_sensitivity_vdw(npv.I_inc_nom, npv.l_nom, npv.L_nom, v, npv.C3_nom, d[i])[0]

# test it
S_vdw_fit = cs.calc_sensitivity_vdw(npv.I_inc_nom, npv.l_nom, npv.L_nom, v, npv.C3_nom)
print S_vdw_fit

leg = plt.plot(d_plot, S_vdw, 'bo')
plt.xlabel('grating period (nm)')
plt.ylabel('sensitivity (rad/s / sqrt(Hz))')
title_str = "v = " + str(v)
plt.title(title_str)
plt.show()




