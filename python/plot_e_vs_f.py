#!/usr/bin/python

import numpy as np
import cmath
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import calc_sensitivity as cs
import nominal_param_values as npv

# plot diffraction efficiencies vs open fraction for vdw and non-vdw models

f = np.linspace(0.01, 0.99, 20)
e0 = np.zeros(f.size)
e1 = np.zeros(f.size)
e0_vdw = np.zeros(f.size)
e1_vdw = np.zeros(f.size)

for i in range(f.size):

    e0[i] = f[i]
    e1[i] = f[i] * np.sin(np.pi * f[i]) / (np.pi * f[i])
    e0_vdw[i] = cmath.polar(cs.calc_diffraction_eff_vdw(0, f[i], npv.d_nom, npv.l_nom, npv.v_nom, npv.C3_nom)[0])[0]
    e1_vdw[i] = cmath.polar(cs.calc_diffraction_eff_vdw(1, f[i], npv.d_nom, npv.l_nom, npv.v_nom, npv.C3_nom)[0])[0]

l1, l2, l3, l4 = plt.plot(f,e0,'bo',    f,e1,'ro',    f,e0_vdw,'bs',     f,e1_vdw,'rs')
plt.xlabel("grating open fraction")
plt.ylabel("|diffraction efficiency|")
plt.legend((l1, l2, l3, l4), ("no vdw, n=0", "no vdw, n=1", "vdw, n=0", "vdw, n=1"), loc='upper left')
plt.show()

