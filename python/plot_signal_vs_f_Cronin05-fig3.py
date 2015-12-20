#!/usr/bin/python

import numpy as np
import cmath
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy.optimize import curve_fit
import calc_sensitivity as cs
import nominal_param_values as npv

# plot the bracketed terms in eqn 14 of Cronin2005 paper vs open fraction,
# i.e. replicate Fig 3 of the Cronin2005 paper

f = np.linspace(0.1,0.9,10)
sig1 = np.zeros(f.size)
sig2 = np.zeros(f.size)
sig1_novdw = np.zeros(f.size)
sig2_novdw = np.zeros(f.size)

d = 100.0e-9
v = 1000.0

print "making graphs"

# calculate signal vs f
for i in range(f.size):

    # WITH vdW 

    # calculate diffraction efficiencies
    e0_g1 = cmath.polar(cs.calc_diffraction_eff_vdw(0, f[i], d, 150.0e-9, v, 4.8e-49)[0])[0]
    e1_g1 = cmath.polar(cs.calc_diffraction_eff_vdw(1, f[i], d, 150.0e-9, v, 4.8e-49)[0])[0]
    e1_g2 = e1_g1

    # calculate the bracketted terms
    sig1[i] = (e0_g1*e1_g1)**2 / (e0_g1**2 + e1_g1**2)
    sig2[i] = e1_g2**2

    # WITHOUT vdW

    # calculate diffraction efficiencies
    e0_g1 = cmath.polar(cs.calc_diffraction_eff_vdw(0, f[i], d, 150.0e-9, v, 0.0)[0])[0]
    e1_g1 = cmath.polar(cs.calc_diffraction_eff_vdw(1, f[i], d, 150.0e-9, v, 0.0)[0])[0]
    e1_g2 = e1_g1

    # calculate the bracketted terms
    sig1_novdw[i] = (e0_g1*e1_g1)**2 / (e0_g1**2 + e1_g1**2)
    sig2_novdw[i] = e1_g2**2

# gaussian fit the result
def gauss(x, A, x0, sig):
    return A * np.exp(-(x-x0)**2/(2*sig**2))

print "fitting curves"
popt1, pcov1 = curve_fit(gauss, f, sig1)
print popt1
popt2, pcov2 = curve_fit(gauss, f, sig2)
print popt2

# make versions of the fit result to plot
sig1_fit = gauss(f, popt1[0], popt1[1], popt1[2])
sig2_fit = gauss(f, popt2[0], popt2[1], popt2[2])


sig_res = cs.calc_signal_vdw(npv.I_inc_nom, d, 150.0e-9, v, 4.8e-49)
print "optimimal open fractions:"
print "f1: ", sig_res[1]
print "f2: ", sig_res[2]
print "gives signal: ", sig_res[0]


plt.figure()
plt.subplot(211)
plt.plot(f, sig1, 'bo', f, sig1_novdw, 'ro', f, sig1_fit, 'b-')
plt.ylabel("bracketed term 1")
plt.xlabel("g1 open fraction")
plt.subplot(212)
plt.plot(f, sig2, 'bo', f, sig2_novdw, 'ro', f, sig2_fit, 'b-')
plt.ylabel("bracketed term 2")
plt.xlabel("g2 open fraction")
plt.show()




