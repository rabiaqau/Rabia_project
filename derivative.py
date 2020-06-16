import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
import math
from sympy import exp

    
x,amplitude,sigma,mean = sym.symbols('x amplitude sigma mean',real=True)
gauss= amplitude * (1.0 / ( math.sqrt(2 * math.pi) * sigma ) ) * exp( - (1.0 / 2.0) * ( (x - mean)**2 / ( sigma **2) ) )



def derivative(parameter):


    derivative_gauss=gauss.diff(parameter)
    

    return derivative_gauss



derivative_peak_amplitude =  derivative(amplitude)                                                                                                                      

#         # sigma                                                                                                                                                                 

derivative_peak_sigma= derivative(sigma)                                                                                                                                

#         # mean                                                                                                                                                                  

derivative_peak_mean= derivative(mean)                                                                                                                                  








 # value of amplitude 1                                                                                                                                
# Value_der_peak_amplitude1 =derivative_peak_amplitude.evalf(subs={amplitude:amplitude1_value,sigma:sigma1_value,mean:mean1_value,x:Xmax})

#         #  print Value_der_peak_amplitude1,"Value_der_peak_amplitude1"                                                                                       

#         # value of amplitude 2                                                                                                                                
# Value_der_peak_amplitude2 =derivative_peak_amplitude.evalf(subs={amplitude:amplitude2_value,sigma:sigma2_value,mean:mean2_value,x:Xmax})
# #        print Value_der_peak_amplitude2,"Value_der_peak_amplitude2"                                                                                          
#         # value of sigma1                                                                                                                                     
# Value_der_peak_sigma1 =derivative_peak_sigma.evalf(subs={ amplitude:amplitude1_value, sigma:sigma1_value, mean:mean1_value, x:Xmax})

# #        print Value_der_peak_sigma1,"Value_der_peak_sigma1"                                                                                                  

# Value_der_peak_sigma2 =derivative_peak_sigma.evalf(subs={ amplitude:amplitude2_value, sigma:sigma2_value, mean:mean2_value, x:Xmax})
# #        print Value_der_peak_sigma2,"Value_der_peak_sigma2"                                                                                                  


# Value_der_peak_mean1 =derivative_peak_mean.evalf(subs={ amplitude:amplitude1_value, sigma:sigma1_value, mean:mean1_value, x:Xmax})

# #        print Value_der_peak_mean1,"Value_der_peak_mean1"                                                                                                   


# Value_der_peak_mean2 =derivative_peak_mean.evalf(subs={ amplitude:amplitude2_value, sigma:sigma2_value, mean:mean2_value, x:Xmax})





#        print Value_der_peak_mean2,"Value_der_peak_mean2"                                                                                                    
def offdiagnol(amplitude1,sigma1,mean1,amplitude2,sigma2,mean2, amp1_mean1,  amp1_sigma1, amp1_amp2, amp1_sigma2, amp1_mean2, mean1_sigma1, mean1_amp2, mean1_sigma2, mean1_mean2, sigma1_amp2, sigma1_sigma2, sigma1_mean2, amp2_sigma2, amp2_mean2, sigma2_mean2 ):


    peak_amplitude1= 2.0 * amplitude1 * mean1 * amp1_mean1 + 2.0 * amplitude1 * sigma1 * amp1_sigma1 + 2.0 * amplitude1 * amplitude2 * amp1_amp2 + 2.0 * amplitude1 * sigma2 * amp1_sigma2 + 2.0 * amplitude1 * mean2 * amp1_mean2




# off dignol elements w.r.t mean1                                                                                                                     
    peak_mean1 = 2.0 * mean1 * sigma1 * mean1_sigma1 + 2.0 * mean1 * amplitude2 * mean1_amp2 + 2.0 * mean1 * sigma2 * mean1_sigma2 + 2.0* mean1 * mean2 * mean1_mean2

# off diagnol elements w.r.t sigma1                                                                                                                   

    peak_sigma1 = 2.0 * sigma1 * amplitude2 * sigma1_amp2 + 2.0 * sigma1 * sigma2 * sigma1_sigma2 + 2.0 * sigma1 * mean2 * sigma1_mean2

        #3                                                                                                                                                    
    peak_amplitude2 = 2.0 * amplitude2 * sigma2 * amp2_sigma2 + 2.0 * amplitude2 * mean2 * amp2_mean2


    peak_sigma2 = 2.0 * sigma2 * mean2 * sigma2_mean2


    off_diagnol_elements= peak_amplitude1 + peak_mean1 + peak_sigma1 + peak_amplitude2 + peak_sigma2

    return off_diagnol_elements












