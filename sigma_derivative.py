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

width_of_beam = 0.399 * (amplitude1 + amplitude2) / ( amplitude1 * (1.0 / ( math.sqrt(2 * math.pi) * sigma1 ) ) * exp( - (1.0 / 2.0) * ( (x - mean1)**2 / ( sigma1 **2) ) ) + amplitude2 * (1.0 / (  math.sqrt(2 * math.pi) * sigma2 ) ) * exp( - (1.0 / 2.0) * ( (x - mean2)**2 / ( sigma2 **2 ) ) ) ) 



derivative_width_beam_amplitude1 = derivative(width_of_beam,amplitude1)
derivative_width_beam_amplitude2= derivative(width_of_beam,amplitude2)
derivative_width_beam_sigma1 =derivative(width_of_beam,sigma1)
derivative_width_beam_sigma2 =derivative(width_of_beam,sigma2)


derivative_width_beam_mean1 =derivative(width_of_beam,mean1)
derivative_width_beam_mean2 =derivative(width_of_beam,mean2)



def derivative_value(derivative,amp1,sig1,mea1,amp2,sig2,mea2,position):


    required_value = derivative.evalf(subs={amplitude1:amp1,sigma1:sig1,mean1:mea1,amplitude2:amp2,sigma2:sig2,mean2:mea2, x:position})


    return required_value









