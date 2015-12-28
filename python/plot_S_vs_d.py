#!/usr/bin/python

import time
import numpy as np
import matplotlib.pyplot as plt
import calc_sensitivity as cs
import nominal_param_values as npv


v = 50.0 # m/s

d = np.linspace(80.0e-9,400.0e-9,30)
d_plot = d*1.0e9
S_vdw = np.zeros(d.size, dtype=np.float)
S_kd = np.zeros(d.size, dtype=np.float)
f1 = np.zeros(d.size, dtype=np.float)
f2 = np.zeros(d.size, dtype=np.float)
f_max = cs.f_max(d)
for i in range(d.size):
    print d[i]

    result = cs.calc_sensitivity_vdw(npv.I_inc_nom, npv.l_nom, npv.L_nom, npv.v_nom, npv.C3_nom, d[i])
    S_vdw[i] = result[0]
    f1[i] = result[2]
    f2[i] = result[3]
    S_kd[i] = cs.calc_sensitivity_kd(npv.I_inc_nom, npv.L_nom, npv.v_nom, d[i])

plt.figure(figsize=plt.figaspect(1.4))

plt.subplot(211)
l1, l2 = plt.plot(d_plot, S_vdw, 'b-', d_plot, S_kd, 'r-')
plt.xlabel('grating period (nm)')
plt.ylabel('sensitivity (rad/s / sqrt(Hz))')
plt.legend((l1, l2), ('material gratings', 'Kapitza-Dirac gratings'), loc='upper right')

plt.subplot(212)
l3, l4, l5 = plt.plot(d_plot, f_max, 'k-', d_plot, f1, 'g-', d_plot, f2, 'm-')
plt.xlabel('grating period (nm)')
plt.ylabel('open fraction')
plt.legend((l3, l4, l5), ('max possible open fraction', 'best g1 open fraction', 'best g2 open fraction'), loc='lower right')

plt.show()







