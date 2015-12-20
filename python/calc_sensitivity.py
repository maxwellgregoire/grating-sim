import cmath
import numpy as np
from scipy import integrate
from scipy.optimize import minimize
from scipy.optimize import curve_fit


# In this code, the quantity |C|^2 * <I> is referred to as "signal"

# Constants
hbar = 1.054e-34

# Calculate sensitivity assuming vdW interactions given:
#   beam intensity incident on the first grating (I_inc)
#   grating thickness (l)
#   longitudinal spacing between gratings (L)
#   atom velocity (v)
#   vdW C_3 coefficient (C3)
# This function (optionally) numerically optimizes:
#   grating open fractions (unless fixed_f is specified, in which case all gratings have that open fraction)
#   grating period (unless fixed_d is specified)
def calc_sensitivity_vdw_optimized(I_inc, l, L, v, C3, fixed_d = None, fixed_f = None):
    """ Calculate sensitivity assuming vdW interactions """

    # Optimize the sensitivity, given by eqn 8 in then notes. 
    #   The |C|^2<I> term in eqn 8 is given by eqn 14 in the Cronin2005 paper.
    #   The third bracketted term in the Cronin2005 paper eqn 14 depends only on 3rd grating open fraction, 
    #       which itself is not present anywhere else in notes eqn 8,
    #       and is optimized when the 3rd grating open fraction = 0.371.
    #       Therefore, the optimized value of the third bracketted term in the Cronin2005 paper is 0.230.

    if fixed_d is None and fixed_f is None:

        # Arguments (i.e. elements of the input array x) are:
        #   [0]: grating period
        #   [1]: grating 1 open fraction
        #   [2]: grating 2 open fraction 
        def sensitivity_to_optimize(x):

            # calculate diffraction efficiencies
            e0_g1 = cmath.polar(calc_diffraction_eff_vdw(0, x[1], x[0], l, v, C3)[0])[0]
            e1_g1 = cmath.polar(calc_diffraction_eff_vdw(1, x[1], x[0], l, v, C3)[0])[0]
            e1_g2 = cmath.polar(calc_diffraction_eff_vdw(1, x[2], x[0], l, v, C3)[0])[0]

            # calculate signal
            signal = 4.0*I_inc*0.23*(e0_g1*e1_g1*e1_g2)**2 / (e0_g1**2 + e1_g1**2)

            return v*x[0] / (4.0*np.pi*L**2*np.sqrt(signal))

        guesses = np.array([100.0e-9, 0.5, 0.5])
        return minimize(sensitivity_to_optimize, guesses, method='L-BFGS-B', bounds=((10.0e-9, None), (0.01, 1.0), (0.01, 1.0)))

    elif fixed_f is None: # and fixed_d is specified

        # Arguments (i.e. elements of the input array x) are:
        #   [0]: grating 1 open fraction
        #   [1]: grating 2 open fraction 
        def sensitivity_to_optimize(x):

            print x

            # calculate diffraction efficiencies
            e0_g1 = cmath.polar(calc_diffraction_eff_vdw(0, x[0], fixed_d, l, v, C3)[0])[0]
            e1_g1 = cmath.polar(calc_diffraction_eff_vdw(1, x[0], fixed_d, l, v, C3)[0])[0]
            e1_g2 = cmath.polar(calc_diffraction_eff_vdw(1, x[1], fixed_d, l, v, C3)[0])[0]

            # calculate signal
            signal = 4.0*I_inc*0.23*(e0_g1*e1_g1*e1_g2)**2 / (e0_g1**2 + e1_g1**2)

            return v*fixed_d / (4.0*np.pi*L**2*np.sqrt(signal))

        guesses = np.array([0.5, 0.5])
        return minimize(sensitivity_to_optimize, guesses, method='L-BFGS-B', bounds=((0.01, 1.0), (0.01, 1.0)))
    
    elif fixed_d is None: # and fixed_f is specified

        # Arguments (i.e. elements of the input array x) are:
        #   [0]: grating period
        def sensitivity_to_optimize(x):

            # calculate diffraction efficiencies
            e0_g1 = cmath.polar(calc_diffraction_eff_vdw(0, fixed_f, x[0], l, v, C3)[0])[0]
            e1_g1 = cmath.polar(calc_diffraction_eff_vdw(1, fixed_f, x[0], l, v, C3)[0])[0]
            e1_g2 = e1_g1

            # calculate signal
            signal = 4.0*I_inc*0.23*(e0_g1*e1_g1*e1_g2)**2 / (e0_g1**2 + e1_g1**2)

            return v*x[0] / (4.0*np.pi*L**2*np.sqrt(signal))

        guesses = np.array([100.0e-9])
        return minimize(sensitivity_to_optimize, guesses, method='L-BFGS-B', bounds=np.array([(10.0e-9, None)]))

    else: # both fixed_f and fixed_d are specified

        # calculate diffraction efficiencies
        e0_g1 = cmath.polar(calc_diffraction_eff_vdw(0, fixed_f, fixed_d, l, v, C3)[0])[0]
        e1_g1 = cmath.polar(calc_diffraction_eff_vdw(1, fixed_f, fixed_d, l, v, C3)[0])[0]
        e1_g2 = e1_g1

        # calculate signal
        signal = 4.0*I_inc*0.23*(e0_g1*e1_g1*e1_g2)**2 / (e0_g1**2 + e1_g1**2)

        return v*fixed_d / (4.0*np.pi*L**2*np.sqrt(signal))





# Calculate the signal given
#   beam intensity incident on the first grating (I_inc)
#   grating period (d)
#   grating thickness (l)
#   atom velocity (v)
#   vdW C_3 coefficient (C3)
# This function numerically optimizes open fractions
# By default, this function optimizes open fractions by calculating terms in eqn 14 in the Cronin2005 paper
#   and fitting a gaussian to it.
#   It sounds rudimentary, but it's about 10x faster than using scipy.optimize.minimize, 
#   and it gets the answer to within +/- 0.01, which is 
#       a) good enough
#       b) not much better than scipy.optimize.minimize, judging by the different answers you get using 
#           different fit algorithms
#   By setting use_gauss_fit to False, you can force it to use scipy.optimize.minimize
def calc_signal_vdw(I_inc, d, l, v, C3, use_gauss_fit = True):
    """ Calculate signal assuming vdW interactions """

    # Optimize the signal, given by eqn 14 in the Cronin2005 paper
    # The first two bracketted terms in that equation can be maximized independently, 
    #   as the first depends on the grating 1 open fraction and the second depends on the grating 2 open fraction
    # The third bracketted term in the Cronin2005 paper eqn 14 depends only on 3rd grating open fraction, 
    #   which itself is not present anywhere else in notes eqn 8,
    #   and is optimized when the 3rd grating open fraction = 0.371.
    #   Therefore, the optimized value of the third bracketted term in the Cronin2005 paper is 0.230.

    def b1_to_optimize(f1):

        # calculate diffraction efficiencies
        e0_g1 = cmath.polar(calc_diffraction_eff_vdw(0, f1, d, l, v, C3)[0])[0]
        e1_g1 = cmath.polar(calc_diffraction_eff_vdw(1, f1, d, l, v, C3)[0])[0]

        # return the negative so that minimize works
        return -(e0_g1*e1_g1)**2 / (e0_g1**2 + e1_g1**2)

    def b2_to_optimize(f2):

        # calculate diffraction efficiency
        e1_g2 = cmath.polar(calc_diffraction_eff_vdw(1, f2, d, l, v, C3)[0])[0]

        # return the negative so that minimize works
        return -e1_g2**2

    if use_gauss_fit:

        # calculate the bracketted terms vs f
        f = np.linspace(0.1,0.9,10)
        sig1 = np.zeros(f.size)
        sig2 = np.zeros(f.size)

        for i in range(f.size):
            sig1[i] = -b1_to_optimize(f[i])
            sig2[i] = -b2_to_optimize(f[i])

        # define a guassian with no y offset 
        def gauss(x, A, x0, sig):
            return A * np.exp(-(x-x0)**2/(2*sig**2))

        # popt is the result, pcov is the covariance matrix
        popt1, pcov1 = curve_fit(gauss, f, sig1)
        popt2, pcov2 = curve_fit(gauss, f, sig2)

        return np.array([4.0*I_inc*0.23 * popt1[0] * popt2[0], popt1[1], popt2[1]])

    else:

        res_b1 = minimize(b1_to_optimize, 0.5, method='SLSQP', bounds=np.array([(0.1,0.9)]), options={'eps':0.01})
        res_b2 = minimize(b2_to_optimize, 0.5, method='SLSQP', bounds=np.array([(0.1,0.9)]), options={'eps':0.01})
    
        # returns
        #   [0]: signal
        #   [1]: optimal grating 1 open fraction
        #   [1]: optimal grating 2 open fraction
        return np.array([4.0*I_inc*0.23 * res_b1.fun * res_b2.fun, res_b1.x[0], res_b2.x[0]])



# Calculate a particular diffraction efficiency of a grating assuming vdW interactions given:
#   diffraction order (n)
#   open fraction (f)
#   period (d)
#   thickness (l)
#   atom velocity (v)
#   vdW C3 coefficient (C3)
# Returns complex value [0] and error [1]
def calc_diffraction_eff_vdw(n, f, d, l, v, C3):
    """ Calculate a particular diffraction efficiency of a grating assuming vdW interactions """

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
def calc_sensitivity_kd(I_inc, L, v, d, assume_f_half = True):
    """ Calculate sensitivity assuming no vdW interaction """

    # calc signal
    if assume_f_half:
        # use w/d = 0.5 in Cronin2005 paper eqn 16
        signal = 0.0059*I_inc
    else:
        # use optimal open fractions, calculated in Cronin2005 paper, that lead to optimal signal
        signal = 0.0070*I_inc
        
    return v*d / (4.0*np.pi*L**2*np.sqrt(signal))
