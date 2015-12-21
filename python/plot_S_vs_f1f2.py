#!/usr/bin/python

import numpy as np
import cmath
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import calc_sensitivity as cs
import nominal_param_values as npv

# plot sensitivity vs g1 and g2 open fractions for a given velocity and grating period

f1 = np.linspace(0.3,0.95,11)
f2 = np.linspace(0.3,0.95,11)
F1, F2 = np.meshgrid(f1,f2)
S = np.zeros([f1.size,f2.size], dtype = np.float)

for i in range(f1.size):
    for j in range(f2.size):

        # calculate diffraction efficiencies
        e0_g1 = cmath.polar(cs.calc_diffraction_eff_vdw(0, f1[i], npv.d_nom, npv.l_nom, npv.v_nom, npv.C3_nom)[0])[0]
        e1_g1 = cmath.polar(cs.calc_diffraction_eff_vdw(1, f1[i], npv.d_nom, npv.l_nom, npv.v_nom, npv.C3_nom)[0])[0]
        e1_g2 = cmath.polar(cs.calc_diffraction_eff_vdw(1, f2[j], npv.d_nom, npv.l_nom, npv.v_nom, npv.C3_nom)[0])[0]

        # calculate signal
        signal = 4.0*npv.I_inc_nom*0.23*(e0_g1*e1_g1*e1_g2)**2 / (e0_g1**2 + e1_g1**2)

        # calculate sensitivity
        S[i][j] = npv.v_nom*npv.d_nom / (4.0*np.pi*npv.L_nom**2*np.sqrt(signal))


fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(F1, F2, S, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=True)
plt.ticklabel_format(style='sci', axis='z', scilimits=(0,0))
ax.set_xlabel("g1 open fraction")
ax.set_ylabel("g2 open fraction")
ax.set_zlabel("Sensitivity (rad/s / sqrt(Hz))")
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()
