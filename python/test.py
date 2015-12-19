#!/usr/bin/python

import numpy as np
import cmath
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import calc_sensitivity as cs
import nominal_param_values as npv

res = cs.calc_sensitivity_vdw_optimized(npv.I_inc_nom, npv.l_nom, npv.L_nom, npv.v_nom, npv.C3_nom, fixed_d = 200.0e-9)
print "result:"
print res.fun
print res.x

