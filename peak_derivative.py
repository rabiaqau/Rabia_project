import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
import math
from sympy import exp





frequency = (1./8.89244e-05) #[s-1] 

def derivative(formula,parameter):


    derivative = formula.diff(parameter)
    

    return derivative











# variables define    
x,amplitude,sigma,mean = sym.symbols('x amplitude sigma mean',real=True)

# gaussian function

gauss= amplitude * (1.0 / ( math.sqrt(2 * math.pi) * sigma ) ) * exp( - (1.0 / 2.0) * ( (x - mean)**2 / ( sigma **2) ) )



# calling the derivative function for peak derivative

# amplitude
derivative_peak_amplitude =  derivative(gauss,amplitude)                                                                                                                      

# sigma                                                                                                                                                                 

derivative_peak_sigma= derivative(gauss,sigma)                                                                                                                                

# mean                                                                                                                                                                  

derivative_peak_mean= derivative(gauss,mean)                                                                                                                                  





def derivative_value_peak(derivative,amp,sig,mea,position):


    required_value = derivative.evalf(subs={amplitude:amp,sigma:sig,mean:mea, x:position})


    return required_value





def set_peak_derivative_value(amplitude1,amplitude2,mean1,mean2,sigma1,sigma2,position):
  # derivative value of amplitude1                                                                                                                                                                                                                                      
    der_peak_amplitude1= derivative_value_peak(derivative_peak_amplitude,amplitude1,sigma1,mean1,position)


    # derivative value of amplitude2                                                                                                                                                                                                                                      
    der_peak_amplitude2= derivative_value_peak(derivative_peak_amplitude,amplitude2,sigma2,mean2,position)
#        print derivative_amp2,"derivative_amp2"                                                                                                                                                                                                                              


    # derivative value of sigma1                                                                                                                                                                                                                                          
    der_peak_sigma1= derivative_value_peak(derivative_peak_sigma,amplitude1,sigma1,mean1,position)



    # derivative value of sigma2                                                                                                                                                                                                                                          
    der_peak_sigma2= derivative_value_peak(derivative_peak_sigma,amplitude2,sigma2,mean2,position)



    # derivative value of mean1                                                                                                                                                                              
    
    der_peak_mean1= derivative_value_peak(derivative_peak_mean,amplitude1,sigma1,mean1,position)

#        print derivative_mean1,"derivative_mean1"                                                                                                                                                                                                                            



    # derivative value of mean2                                                                                                                                                                                                                                           
    der_peak_mean2= derivative_value_peak(derivative_peak_mean,amplitude2,sigma2,mean2,position)
    


    return der_peak_amplitude1, der_peak_amplitude2, der_peak_sigma1, der_peak_sigma2, der_peak_mean1, der_peak_mean1













# function for computing the off diagnol elements

def offdiagnol(amplitude1,sigma1,mean1,amplitude2,sigma2,mean2, amp1_mean1,  amp1_sigma1, amp1_amp2, amp1_sigma2, amp1_mean2, mean1_sigma1, mean1_amp2, mean1_sigma2, mean1_mean2, sigma1_amp2, sigma1_sigma2, sigma1_mean2, amp2_sigma2, amp2_mean2, sigma2_mean2 ):



    # peak amplitude1
    off_diagnol_amplitude1= 2.0 * amplitude1 * mean1 * amp1_mean1 + 2.0 * amplitude1 * sigma1 * amp1_sigma1 + 2.0 * amplitude1 * amplitude2 * amp1_amp2 + 2.0 * amplitude1 * sigma2 * amp1_sigma2 + 2.0 * amplitude1 * mean2 * amp1_mean2




# off dignol elements w.r.t mean1                                                                                                                     
    off_diagnol_mean1 = 2.0 * mean1 * sigma1 * mean1_sigma1 + 2.0 * mean1 * amplitude2 * mean1_amp2 + 2.0 * mean1 * sigma2 * mean1_sigma2 + 2.0* mean1 * mean2 * mean1_mean2

# off diagnol elements w.r.t sigma1                                                                                                                   

    off_diagnol_sigma1 = 2.0 * sigma1 * amplitude2 * sigma1_amp2 + 2.0 * sigma1 * sigma2 * sigma1_sigma2 + 2.0 * sigma1 * mean2 * sigma1_mean2

        #3                                                                                                                                                    
    off_diagnol_amplitude2 = 2.0 * amplitude2 * sigma2 * amp2_sigma2 + 2.0 * amplitude2 * mean2 * amp2_mean2


    off_diagnol_sigma2 = 2.0 * sigma2 * mean2 * sigma2_mean2


    off_diagnol_elements= off_diagnol_amplitude1 + off_diagnol_mean1 + off_diagnol_sigma1 + off_diagnol_amplitude2 + off_diagnol_sigma2

    return off_diagnol_elements




# defining the sigma formulas


def diagnol_FOM(derivative_peak_parameter,intensities,luminosity,sigma_x,sigma_y,derivative_of_sigma_x_parameter,derivative_of_sigma_y_parameter):


    constant= (frequency * intensities)/ ( 2 * math.pi *( sigma_x * sigma_y)**2 )

    derivative_of_FOM_parameter = derivative_peak_parameter + constant * ( derivative_of_sigma_x_parameter * sigma_x +
 derivative_of_sigma_y_parameter * sigma_y )


    return derivative_of_FOM_parameter










