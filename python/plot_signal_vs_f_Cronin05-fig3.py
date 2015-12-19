#!/usr/bin/python

import numpy as np
import cmath
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import calc_sensitivity as cs
import nominal_param_values as npv

# plot the bracketed terms in eqn 14 of Cronin2005 paper vs open fraction,
# i.e. replicate Fig 3 of the Cronin2005 paper

f = np.linspace(0.05,0.95,20)
sig1 = np.zeros(f.size)
sig2 = np.zeros(f.size)
sig1_novdw = np.zeros(f.size)
sig2_novdw = np.zeros(f.size)

for i in range(f.size):

    # WITH vdW 

    # calculate modulus-squared diffraction efficiencies
    e0_g1 = cmath.polar(cs.calc_diffraction_eff_vdw(0, f[i], 100.0e-9, 150.0e-9, 1000.0, 4.8e-49)[0])[0]
    e1_g1 = cmath.polar(cs.calc_diffraction_eff_vdw(1, f[i], 100.0e-9, 150.0e-9, 1000.0, 4.8e-49)[0])[0]
    e1_g2 = e1_g1

    # calculate the bracketted terms
    sig1[i] = (e0_g1*e1_g1)**2 / (e0_g1**2 + e1_g1**2)
    sig2[i] = e1_g2**2

    # WITHOUT vdW

    # calculate modulus-squared diffraction efficiencies
    e0_g1 = cmath.polar(cs.calc_diffraction_eff_vdw(0, f[i], 100.0e-9, 150.0e-9, 1000.0, 0.0)[0])[0]
    e1_g1 = cmath.polar(cs.calc_diffraction_eff_vdw(1, f[i], 100.0e-9, 150.0e-9, 1000.0, 0.0)[0])[0]
    e1_g2 = e1_g1

    # calculate the bracketted terms
    sig1_novdw[i] = (e0_g1*e1_g1)**2 / (e0_g1**2 + e1_g1**2)
    sig2_novdw[i] = e1_g2**2

plt.figure()
plt.subplot(211)
plt.plot(f, sig1, 'bo', f, sig1_novdw, 'ro')
plt.ylabel("bracketed term 1")
plt.xlabel("g1 open fraction")
plt.subplot(212)
plt.plot(f, sig2, 'bo', f, sig2_novdw, 'ro')
plt.ylabel("bracketed term 2")
plt.xlabel("g2 open fraction")
plt.show()




