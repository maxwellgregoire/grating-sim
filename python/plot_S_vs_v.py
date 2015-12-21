#!/usr/bin/python

import time
import numpy as np
import matplotlib.pyplot as plt
import calc_sensitivity as cs
import nominal_param_values as npv

# NOTE: unreliable

v = np.linspace(20.0,200.0,20)
S = np.zeros(v.size, dtype=np.float)
d = np.zeros(v.size, dtype=np.float)
f1 = np.zeros(v.size, dtype=np.float)
f2 = np.zeros(v.size, dtype=np.float)
for i in range(v.size):
    print v[i]
    S_result = cs.calc_sensitivity_vdw(npv.I_inc_nom, npv.l_nom, npv.L_nom, v[i], npv.C3_nom)
    S[i] = S_result[0]
    d[i] = S_result[1]
    f1[i] = S_result[2]
    f2[i] = S_result[3]
    print S_result

d_plot = d*1.0e9

plt.figure(figsize=(800,1200))
plt.subplot(311)
plt.plot(v, S, 'bo')
plt.ylabel("sensitivity (rad/s / sqrt(Hz))")
plt.xlabel("velocity (m/s)")
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.subplot(312)
plt.plot(v, d_plot, 'bo')
plt.ylabel("optimal grating period (nm)")
plt.xlabel("velocity (m/s)")
plt.subplot(313)
l1, l2 = plt.plot(v, f1, 'rs', v, f2, 'go')
plt.ylabel("optimal grating open fraction")
plt.xlabel("velocity (m/s)")
plt.legend((l1, l2), ('grating 1', 'grating 2'), loc='upper right')
plt.show()


