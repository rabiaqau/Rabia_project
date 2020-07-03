import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
import math
from sympy import exp

def derivative(formula,parameter):


    derivative_formula = formula.diff(parameter)


    return derivative_formula


# variables define                                                                                                     
x,amplitude1,sigma1,mean1,amplitude2,sigma2,mean2 = sym.symbols(' x amplitude1 sigma1 mean1 amplitude2 sigma2 mean2',real=True)



#cap sigma formula


# no bug in formula "checked"


width_of_beam = 0.4 * (amplitude1 + amplitude2) / ( amplitude1 * (1.0 / ( math.sqrt(2 * math.pi) * sigma1 ) ) * exp( - (1.0 / 2.0) * ( (x - mean1)**2 / ( sigma1 **2) ) ) + amplitude2 * (1.0 / (  math.sqrt(2 * math.pi) * sigma2 ) ) * exp( - (1.0 / 2.0) * ( (x - mean2)**2 / ( sigma2 **2 ) ) ) ) 



#print width_of_beam,"width_of_beam"


derivative_width_beam_amplitude1 = derivative(width_of_beam,amplitude1)
derivative_width_beam_amplitude2= derivative(width_of_beam,amplitude2)
derivative_width_beam_sigma1 =derivative(width_of_beam,sigma1)
derivative_width_beam_sigma2 =derivative(width_of_beam,sigma2)
derivative_width_beam_mean1 =derivative(width_of_beam,mean1)
derivative_width_beam_mean2 =derivative(width_of_beam,mean2)




def set_Capsigma_derivative_value(amplitude1_fit, amplitude2_fit, sigma1_fit, sigma2_fit, mean1_fit, mean2_fit, position):


    #amplitude1(checked)
    amplitude1_value = derivative_width_beam_amplitude1.evalf(subs={amplitude1:amplitude1_fit, sigma1:sigma1_fit, mean1:mean1_fit, amplitude2:amplitude2_fit, sigma2:sigma2_fit, mean2:mean2_fit, x:position})


    # amplitude2(checked)

    amplitude2_value = derivative_width_beam_amplitude2.evalf(subs={amplitude1:amplitude1_fit, sigma1:sigma1_fit, mean1:mean1_fit ,amplitude2:amplitude2_fit, sigma2:sigma2_fit ,mean2:mean2_fit, x:position})


    # sigma1(checked)

    sigma1_value = derivative_width_beam_sigma1.evalf(subs={amplitude1:amplitude1_fit, sigma1:sigma1_fit ,mean1:mean1_fit ,amplitude2:amplitude2_fit ,sigma2:sigma2_fit ,mean2:mean2_fit, x:position})


    # sigma2
    sigma2_value = derivative_width_beam_sigma2.evalf(subs={amplitude1:amplitude1_fit, sigma1:sigma1_fit ,mean1:mean1_fit ,amplitude2:amplitude2_fit ,sigma2:sigma2_fit ,mean2:mean2_fit, x:position})


    # mean1
    mean1_value = derivative_width_beam_mean1.evalf(subs={amplitude1:amplitude1_fit, sigma1:sigma1_fit ,mean1:mean1_fit ,amplitude2:amplitude2_fit, sigma2:sigma2_fit ,mean2:mean2_fit, x:position})


    # mean2
    mean2_value = derivative_width_beam_mean2.evalf(subs={amplitude1:amplitude1_fit,sigma1:sigma1_fit ,mean1:mean1_fit ,amplitude2:amplitude2_fit, sigma2:sigma2_fit ,mean2:mean2_fit, x:position})


    return amplitude1_value, amplitude2_value, sigma1_value, sigma2_value, mean1_value, mean2_value












