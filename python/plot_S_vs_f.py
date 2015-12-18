#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import calc_sensitivity as cs
import nominal_param_values as npv

# plot sensitivity vs g1 and g2 open fractions for a given velocity and grating period

d1 = np.linspace(0.3,0.7,9)
d2 = np.linspace(0.3,0.7,9)
D1, D2 = np.meshgrid(d1,d2)
S = np.zeros([d2.size,d2.size], dtype = np.float)

for i in d1:
    for j in d2:
        pass
        #S[i][j] = 




