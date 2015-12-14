import numpy as np

# Here we refer to the quantity C^2<I> as "signal"

# Calculate signal given
#   beam intensity incident on the first grating
#   diffraction efficiency, grating 1, order 0
#   diffraction efficiency, grating 1, order 1
#   diffraction efficiency, grating 2, order 1
#   open fraction, grating 3
def calc_signal(I_inc, e0_g1, e1_g1, e1_g2, of_g3)
    return 4*I_inc * (e0_g1*e1_g1)**2/(e0_g1**2+e1_g1**2) * e1_g2 * of_g3*(np.sinc(np.pi*of_g3))**2 

# Calculate a particular diffraction efficiency of a grating assuming vdW interactions given
#   diffraction order
#   open fraction
#   period
#   thickness
#   vdW C3 coefficient
#   atom velocity
