#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from calc_sensitivity import calc_diffraction_eff_vdw

#calc_diffraction_eff_vdw(1, 0.5, 100.0e-9, 125.0e-9, 2000, 8.01e-49)

f = np.arange(0.4, 0.9, 0.04)
eff_re = np.zeros(f.size, dtype=np.float)
eff_im = np.zeros(f.size, dtype=np.float)
eff_err_re = np.zeros(f.size, dtype=np.float)
eff_err_im = np.zeros(f.size, dtype=np.float)
for i in range(f.size):
    result_re = calc_diffraction_eff_vdw(1, f[i], 100.0e-9, 125.0e-9, 2000, 8.01e-49)
    eff_re[i] = result_re[0].real
    eff_err_re[i] = result_re[1].real
    result_im = calc_diffraction_eff_vdw(1, f[i], 100.0e-9, 125.0e-9, 2000, 8.01e-49)
    eff_im[i] = result_im[0].imag
    eff_err_im[i] = result_im[1].imag

plt.plot(f, eff_re, 'bs', f, eff_im, 'rs')
plt.show()

#plt.errorbar(f, eff_re, eff_err_re)
#plt.show()

