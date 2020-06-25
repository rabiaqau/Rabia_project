import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
import math
from sympy import exp





frequency = (1./8.89244e-05) #[s-1] 

def derivative(formula,parameter):


    derivative = formula.diff(parameter)
    

    return derivative




def derivative_value_peak(derivative,amp,sig,mea,position):


    required_value = derivative.evalf(subs={amplitude:amp,sigma:sig,mean:mea, x:position})


    return required_value







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



# variables define    
amplitude1, sigma1, mean1, amplitude2, sigma2, mean2 = sym.symbols(' amplitude1 sigma1 mean1 amplitude2 sigma2 mean2',real=True)



# sigma XY formula
sigma_XY = (amplitude1 + amplitude2) /( amplitude1 * (1.0 / ( math.sqrt(2 * math.pi) * sigma1 ) ) * exp( - (1.0 / 2.0) * ( (x - mean1)**2 / ( sigma1 **2) ) ) + amplitude2 * (1.0 / ( math.sqrt(2 * math.pi) * sigma2 ) ) * exp( - (1.0 / 2.0) * ( (x - mean2)**2 / ( sigma2 **2) ) ) )



# calling the derivative function

# amplitude
derivative_sigma_XY_amplitude1 =  derivative(sigma_XY,amplitude1)                             

#print derivative_sigma_XY_amplitude1,"derivative on other module"                                                                 

# # sigma                                                                                                                                                                 

# derivative_peak_sigma= derivative(gauss,sigma)                                                                                                                                

# # mean                                                                                                                                                                  

# derivative_peak_mean= derivative(gauss,mean)                                                                          








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










