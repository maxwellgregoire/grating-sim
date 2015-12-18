#!/usr/bin/python

import time
import numpy as np
import matplotlib.pyplot as plt
import calc_sensitivity as cs
import nominal_param_values as npv


# Replicates Fig 5 of the Cronin 2005 paper, "Limitations on nanotechnology..."

d = np.arange(10.0e-9,150.0e-9,5.0e-9)
d_plot = d*1.0e9
S_vdw = np.zeros(d.size, dtype=np.float)
S_kd = np.zeros(d.size, dtype=np.float)
for i in range(d.size):
    S_vdw[i] = cs.calc_sensitivity_vdw_optimized(npv.I_inc_nom, 150.0e-9, 1.0, 1000.0, 4.8e-49, d[i], 0.5)
    S_kd[i] = cs.calc_sensitivity_kd(npv.I_inc_nom, 1.0, 1000.0, d[i]) 

S_vdw_inv_sq = 1/S_vdw**2
S_kd_inv_sq = 1/S_kd**2

S_vdw_inv_sq_norm = S_vdw_inv_sq / np.linalg.norm(S_vdw_inv_sq)
S_kd_inv_sq_norm = S_kd_inv_sq / np.linalg.norm(S_kd_inv_sq)

plt.plot(d_plot, S_vdw_inv_sq_norm, 'bo', d_plot, S_kd_inv_sq_norm, 'ro')
plt.xlabel('grating period (nm)')
plt.ylabel('normalized inverse-squared sensitivity')
plt.show()




