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

width_of_beam = 0.4 * (amplitude1 + amplitude2) / ( amplitude1 * (1.0 / ( math.sqrt(2 * math.pi) * sigma1 ) ) * exp( - (1.0 / 2.0) * ( (x - mean1)**2 / ( sigma1 **2) ) ) + amplitude2 * (1.0 / (  math.sqrt(2 * math.pi) * sigma2 ) ) * exp( - (1.0 / 2.0) * ( (x - mean2)**2 / ( sigma2 **2 ) ) ) ) 



derivative_width_beam_amplitude1 = derivative(width_of_beam,amplitude1)
derivative_width_beam_amplitude2= derivative(width_of_beam,amplitude2)
derivative_width_beam_sigma1 =derivative(width_of_beam,sigma1)
derivative_width_beam_sigma2 =derivative(width_of_beam,sigma2)
derivative_width_beam_mean1 =derivative(width_of_beam,mean1)
derivative_width_beam_mean2 =derivative(width_of_beam,mean2)
print derivative_width_beam_mean1,"derivative_width_beam_mean1"
print derivative_width_beam_mean2,"derivative_width_beam_mean2"
 




def set_derivative_value(amplitude1_fit, amplitude2_fit, sigma1_fit, sigma2_fit, mean1_fit, mean2_fit, position):

    amplitude1_value = derivative_width_beam_amplitude1.evalf(subs={amplitude1:amplitude1_fit, sigma1:sigma1_fit, mean1:mean1_fit, amplitude2:amplitude2_fit, sigma2:sigma2_fit, mean2:mean2_fit, x:position})

    amplitude2_value = derivative_width_beam_amplitude2.evalf(subs={amplitude1:amplitude1_fit, sigma1:sigma1_fit, mean1:mean1_fit ,amplitude2:amplitude2_fit, sigma2:sigma2_fit ,mean2:mean2_fit, x:position})


    sigma1_value = derivative_width_beam_sigma1.evalf(subs={amplitude1:amplitude1_fit, sigma1:sigma1_fit ,mean1:mean1_fit ,amplitude2:amplitude2_fit ,sigma2:sigma2_fit ,mean2:mean2_fit, x:position})

    sigma2_value = derivative_width_beam_sigma2.evalf(subs={amplitude1:amplitude1_fit, sigma1:sigma1_fit ,mean1:mean1_fit ,amplitude2:amplitude2_fit ,sigma2:sigma2_fit ,mean2:mean2_fit, x:position})

    mean1_value = derivative_width_beam_mean1.evalf(subs={amplitude1:amplitude1_fit, sigma1:sigma1_fit ,mean1:mean1_fit ,amplitude2:amplitude2_fit, sigma2:sigma2_fit ,mean2:mean2_fit, x:position})

    mean2_value = derivative_width_beam_mean2.evalf(subs={amplitude1:amplitude1_fit,sigma1:sigma1_fit ,mean1:mean1_fit ,amplitude2:amplitude2_fit, sigma2:sigma2_fit ,mean2:mean2_fit, x:position})


    return amplitude1_value, amplitude2_value, sigma1_value, sigma2_value, mean1_value, mean2_value




 # Value_der_Capsigma_amplitude1= derivative_value(derivative_width_beam_amplitude1,amplitude1_value,sigma1_value,mean1_value,amplitude2_value,sigma2_value,mean2_value,peak_position)


# #        print Value_der_Capsigma_amplitude1,"derivative_amp1"                                                                                                                                                    



#         Value_der_Capsigma_amplitude2= derivative_value(derivative_width_beam_amplitude2,amplitude1_value,sigma1_value,mean1_value,amplitude2_value,sigma2_value,mean2_value,peak_position)
# #        print derivative_amp2,"derivative_amp2"                                                                                                                                                                  



#         Value_der_Capsigma_sigma1= derivative_value(derivative_width_beam_sigma1,amplitude1_value,sigma1_value,mean1_value,amplitude2_value,sigma2_value,mean2_value,peak_position)




#         Value_der_Capsigma_sigma2= derivative_value(derivative_width_beam_sigma2,amplitude1_value,sigma1_value,mean1_value,amplitude2_value,sigma2_value,mean2_value,peak_position)




#         Value_der_Capsigma_mean1= derivative_value(derivative_width_beam_mean1,amplitude1_value,sigma1_value,mean1_value,amplitude2_value,sigma2_value,mean2_value,peak_position)

# #        print derivative_mean1,"derivative_mean1"                                                                                                                                                                


#         Value_der_Capsigma_mean2= derivative_value(derivative_width_beam_mean2,amplitude1_value,sigma1_value,mean1_value,amplitude2_value,sigma2_value,mean2_value,peak_position)






#     return required_value









