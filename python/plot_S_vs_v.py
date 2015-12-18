#!/usr/bin/python

import time
import numpy as np
import matplotlib.pyplot as plt
import calc_sensitivity as cs
import nominal_param_values as npv

# NOTE: currently assuming optimal open fractions are 0.5

v = np.arange(25,200,25)
S = np.zeros(v.size, dtype=np.float)
d = np.zeros(v.size, dtype=np.float)
for i in range(v.size):
    print v[i]
    S_result = cs.calc_sensitivity_vdw_optimized(npv.I_inc_nom, npv.l_nom, npv.L_nom, v[i], npv.C3_nom, fixed_f = 0.5)
    S[i] = S_result.fun
    d[i] = S_result.x[0]

d_plot = d*1.0e9

plt.figure(num=1, figsize=(8,12))
plt.subplot(211)
plt.plot(v, S, 'bs')
plt.ylabel("sensitivity (rad/s / sqrt(Hz))")
plt.xlabel("velocity (m/s)")
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.subplot(212)
plt.plot(v, d_plot, 'rs')
plt.ylabel("optimal grating period (nm)")
plt.xlabel("velocity (m/s)")
plt.show()


