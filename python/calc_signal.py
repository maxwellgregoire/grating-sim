import numpy as np
from scipy import integrate

# In this code, the quantity |C|^2 * <I> is referred to as "signal"

# Constants:
hbar = 1.054e-34

# Calculate sensitivity assuming vdW interactions given:
#   beam intensity incident on the first grating (I_inc)
#   grating thickness (l)
#   longitudinal spacing between gratings (L)
#   atom velocity (v)
#   vdW C_3 coefficient between the atoms and the silicon-nitride gratings (C3)
# This function numerically optimizes:
#   grating open fractions
#   grating period
# TODO: we may also want to write a version of this function that doesn't optimize period and open fractions
def calc_sensitivity_vdw_optimized(I_inc, l, L, v, C3):

    # Optimize the first two bracketted terms in Cronin2005 paper eqn 14
    #   (The third bracketted term depends only on 3rd grating open fraction, and is optimized for f_g3 = 0.37) 


    pass



# Calculate a particular diffraction efficiency of a grating assuming vdW interactions given:
#   diffraction order (n)
#   open fraction (f)
#   period (d)
#   thickness (l)
#   atom velocity (v)
#   vdW C3 coefficient (C3)
# Returns value [0] and error [1]
def calc_diffraction_eff_vdw(n, f, d, l, v, C3):

    p1 = 2.0*np.pi*n
    p2 = f/2.0
    p3 = C3*l/(hbar*v*d**3)

    # integrate real and imaginary parts separately
    result_re = integrate.quad(lambda y: np.cos(p1*y + p3/(p2-y)**3 + p3/(p2+y)**3), -p2, p2, 
            limit=500, full_output=1)
    result_im = integrate.quad(lambda y: np.sin(p1*y + p3/(p2-y)**3 + p3/(p2+y)**3), -p2, p2, 
            limit=500, full_output=1)

    return np.array([complex(result_re[0], result_im[0]), complex(result_re[1], result_im[1])])



# Calculate sensitivity assuming no vdW interactions (i.e. assuming Kapitza-Dirac gratings) given:
#   beam intensity incident on the first grating (I_inc)
#   longitudinal spacing between gratings (L)
#   grating period (d)
#   atom velocity (v)
# If last argument is true, assume all grating open fractions are 50% (as they would be with 
#   monochromatic standing waves of light) rather than using optimal open fractions
# NOTE: For kapitza-dirac gratings, sensitivity depends linearly on grating period, 
#   so we cannot optimize w/r/t grating period as we would assuming vdW interactions
def calc_sensitivity_kd(I_inc, L, d, v, assume_f_half = True):

    # calc signal
    if assume_f_half:
        # use w/d = 0.5 in Cronin2005 paper eqn 16
        signal = 0.0059*I_inc
    else:
        # use optimal open fractions, calculated in Cronin2005 paper, that lead to optimal signal
        signal = 0.0070*I_inc
        
    return v*d / (4.0*np.pi*L**2*np.sqrt(signal))
