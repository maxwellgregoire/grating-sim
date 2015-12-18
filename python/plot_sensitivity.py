#!/usr/bin/python

import time
import numpy as np
import matplotlib.pyplot as plt
from calc_sensitivity import calc_sensitivity_vdw_optimized
from calc_sensitivity import calc_sensitivity_kd

# nominal values for plotting
I_inc_nom = 1.0e7 # atoms/s i.e. Hz
l_nom = 150.0e-9 # m
L_nom = 0.05 # m
v_nom = 1000.0 # m/s
C3_nom = 8.85e-49 # J*m^3 (SI units)
f_nom = 0.5 

# Assume optimal open fractions are 50% (So far, this has shown to be exactly true for the first two gratings)
def plot_sensitivity_vs_v():
    v = np.arange(25,200,25)
    S = np.zeros(v.size, dtype=np.float)
    d = np.zeros(v.size, dtype=np.float)
    for i in range(v.size):
        print v[i]
        S_result = calc_sensitivity_vdw_optimized(I_inc_nom, l_nom, L_nom, v[i], C3_nom, fixed_f = 0.5)
        S[i] = S_result.fun
        d[i] = S_result.x[0]

    d_plot = d*1.0e9

    plt.figure(num=1, figsize=(8,12))
    plt.subplot(211)
    plt.plot(v, S, 'bs')
    plt.ylabel("sensitivity (rad/s / sqrt(Hz))")
    plt.xlabel("velocity (m/s)")
    plt.subplot(212)
    plt.plot(v, d_plot, 'rs')
    plt.ylabel("optimal grating period (nm)")
    plt.xlabel("velocity (m/s)")
    plt.show()



def plot_sensitivity_vs_d_fixed_f():
    d = np.arange(10.0e-9,150.0e-9,10.0e-9)
    d_plot = d*1.0e9
    S_fixed_f = np.zeros(d.size, dtype=np.float)
    for i in range(d.size):
        S_fixed_f[i] = calc_sensitivity_vdw_optimized(I_inc_nom, l_nom, L_nom, v_nom, C3_nom, d[i], f_nom)

    plt.plot(d_plot, S_fixed_f, 'bo')
    plt.xlabel('grating period (nm)')
    plt.ylabel('Sensitivity (rad/s) / sqrt(Hz)')
    plt.show()




def plot_sensitivity_vs_d():
    d = np.arange(10.0e-9,150.0e-9,20.0e-9)
    d_plot = d*1.0e9
    S_fixed_f = np.zeros(d.size, dtype=np.float)
    S = np.zeros(d.size, dtype=np.float)
    for i in range(d.size):
        print d[i]
        S_fixed_f[i] = calc_sensitivity_vdw_optimized(I_inc_nom, l_nom, L_nom, v_nom, C3_nom, d[i], f_nom)
        S_result = calc_sensitivity_vdw_optimized(I_inc_nom, l_nom, L_nom, v_nom, C3_nom, d[i])
        S[i] = S_result.fun
        print "f_g1: ", S_result.x[0]
        print "f_g2: ", S_result.x[1]
        
    plt.plot(d_plot, S_fixed_f, 'b-', d_plot, S, 'ro')
    plt.xlabel('grating period (nm)')
    plt.ylabel('Sensitivity (rad/s / sqrt(Hz))')
    plt.show()



plot_sensitivity_vs_v()

#S_min_result = calc_sensitivity_vdw_optimized(I_inc_nom, l_nom, L_nom, v_nom, C3_nom, fixed_f = 0.5)
#print "Min S: ", S_min_result.fun
#print "d: ", S_min_result.x[0]
#
#S_min_result = calc_sensitivity_vdw_optimized(I_inc_nom, l_nom, L_nom, v_nom, C3_nom)
#print "Min S: ", S_min_result.fun
#print "d: ", S_min_result.x

#plot_sensitivity_vs_d_fixed_f()

#plot_sensitivity_vs_d()


