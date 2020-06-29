import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
import math
from sympy import exp




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








def derivative_FOM(formula,parameter):


    derivative_formula = formula.diff(parameter)


    return derivative_formula

# constant is (2*pi / n1n2 * 1e22 * 8 * 1e-24)


y,n1n2, amplitude1_x, amplitude2_x, amplitude1_y, amplitude2_y, sigma1_y, mean1_y, amplitude2_y, sigma2_y,mean2_y = sym.symbols(' y n1n2 amplitude1_x amplitude2_x amplitude1_y amplitude2_y sigma1_y mean1_y amplitude2_y sigma2_y mean2_y',real=True)


FOM = ( 12.5 * (amplitude1_x + amplitude2_x) * (amplitude1_y + amplitude2_y)  ) / n1n2 / ( amplitude1_y * (1.0 / ( math.sqrt(2 * math.pi) * sigma1_y ) ) * exp( - (1.0 / 2.0) * ( (y - mean1_y)**2 / ( sigma1_y **2) ) ) +  amplitude2_y * (1.0 / (  math.sqrt(2 * math.pi) * sigma2_y ) ) * exp( - (1.0 / 2.0) * ( (y - mean2_y)**2 / ( sigma2_y **2 ) ) )  )

#  defining vaiables




derivative_FOM_amplitude1_x=derivative_FOM(FOM,amplitude1_x)
derivative_FOM_amplitude2_x=derivative_FOM(FOM,amplitude2_x)



#derivative_FOM_amplitude1_x=derivative_FOM(FOM,amplitude1_x)


def values_derivative_X(amp1,sigma1,mean1,amp2,sigma2,mean2,position,current,x_error_on_A1,x_error_on_A2,cov_A1_A2):


    derivative_FOM_amplitude1_x=derivative_FOM(FOM,amplitude1_x)
    derivative_FOM_amplitude2_x=derivative_FOM(FOM,amplitude2_x)


    x_derivtaive_FOM_A1 = derivative_FOM_amplitude1_x.evalf(subs={ amplitude1_y:amp1 ,sigma1_y:sigma1 , mean1_y:mean1 ,amplitude2_y:amp2 , sigma2_y: sigma2 ,mean2_y: mean2 , y:position ,n1n2:current})


    x_derivtaive_FOM_A2 = derivative_FOM_amplitude2_x.evalf(subs={ amplitude1_y:amp1 ,sigma1_y:sigma1 , mean1_y:mean1 ,amplitude2_y:amp2 , sigma2_y: sigma2 ,mean2_y: mean2 , y:position ,n1n2:current})



    dignol_terms =( x_derivtaive_FOM_A1 * x_error_on_A1) **2 + ( x_derivtaive_FOM_A2 * x_error_on_A2) **2
    

    off_dignol =2 * (x_derivtaive_FOM_A1)* (x_derivtaive_FOM_A2) *cov_A1_A2

    X_terms = dignol_terms + off_dignol


    return X_terms



# def Set_FOM_Values:

#     derivative_FOM_amplitude1_x=derivative_FOM(FOM,amplitude1_x)
#     derivative_FOM_amplitude2_x=derivative_FOM(FOM,amplitude2_x)




def values_derivative_Y(amp1_y,sig1_y,y_mea1,amp2_y,y_sig2,y_mea2,amp1_x,amp2_x,position,current,error_amp1_y, error_amp2_y, error_sig1_y, error_sig2_y, error_mea1_y, error_mea2_y, cov_amp1_mean1,  cov_amp1_sigma1, cov_amp1_amp2, cov_amp1_sigma2, cov_amp1_mean2, cov_mean1_sigma1, cov_mean1_amp2, cov_mean1_sigma2, cov_mean1_mean2, cov_sigma1_amp2, cov_sigma1_sigma2, cov_sigma1_mean2, cov_amp2_sigma2, cov_amp2_mean2, cov_sigma2_mean2):
  
    derivative_FOM_amplitude1_y=derivative_FOM(FOM,amplitude1_y)
    derivative_FOM_amplitude2_y=derivative_FOM(FOM,amplitude2_y)
    derivative_FOM_sigma1_y=derivative_FOM(FOM,sigma1_y)
    derivative_FOM_sigma2_y=derivative_FOM(FOM,sigma2_y)
    derivative_FOM_mean1_y=derivative_FOM(FOM,mean1_y)
    derivative_FOM_mean2_y=derivative_FOM(FOM,mean2_y)
    

    y_FOM_derivative_amplitude1= derivative_FOM_amplitude1_y.evalf(subs={ amplitude1_x:amp1_x, amplitude2_x:amp2_x ,amplitude1_y:amp1_y ,sigma1_y:sig1_y , mean1_y:y_mea1 ,amplitude2_y:amp2_y ,sigma2_y: y_sig2 ,mean2_y: y_mea2 , y:position ,n1n2:current})



    y_FOM_derivative_amplitude2=derivative_FOM_amplitude2_y.evalf(subs={ amplitude1_x:amp1_x, amplitude2_x:amp2_x ,amplitude1_y:amp1_y ,sigma1_y:sig1_y , mean1_y:y_mea1 ,amplitude2_y:amp2_y ,sigma2_y: y_sig2 ,mean2_y: y_mea2 , y:position ,n1n2:current})

    
    y_FOM_derivative_sigma1=derivative_FOM_sigma1_y.evalf(subs={ amplitude1_x:amp1_x, amplitude2_x:amp2_x ,amplitude1_y:amp1_y ,sigma1_y:sig1_y , mean1_y:y_mea1 ,amplitude2_y:amp2_y ,sigma2_y: y_sig2 ,mean2_y: y_mea2 , y:position ,n1n2:current})

    y_FOM_derivative_sigma2=derivative_FOM_sigma2_y.evalf(subs={ amplitude1_x:amp1_x, amplitude2_x:amp2_x ,amplitude1_y:amp1_y ,sigma1_y:sig1_y , mean1_y:y_mea1 ,amplitude2_y:amp2_y ,sigma2_y: y_sig2 ,mean2_y: y_mea2 , y:position ,n1n2:current})


    y_FOM_derivative_mean1=derivative_FOM_mean1_y.evalf(subs={ amplitude1_x:amp1_x, amplitude2_x:amp2_x ,amplitude1_y:amp1_y ,sigma1_y:sig1_y , mean1_y:y_mea1 ,amplitude2_y:amp2_y ,sigma2_y: y_sig2 ,mean2_y: y_mea2 , y:position ,n1n2:current})


    y_FOM_derivative_mean2=derivative_FOM_mean2_y.evalf(subs={ amplitude1_x:amp1_x, amplitude2_x:amp2_x ,amplitude1_y:amp1_y ,sigma1_y:sig1_y , mean1_y:y_mea1 ,amplitude2_y:amp2_y ,sigma2_y: y_sig2 ,mean2_y: y_mea2 , y:position ,n1n2:current})


    y_diagnol_terms=(y_FOM_derivative_amplitude1 *error_amp1_y)**2 + (y_FOM_derivative_amplitude2*error_amp2_y)**2 +(y_FOM_derivative_sigma1*error_sig1_y)**2+(y_FOM_derivative_sigma2 *error_sig2_y)**2 +(y_FOM_derivative_mean1 *error_mea1_y)**2 + (y_FOM_derivative_mean2*error_mea2_y)**2




    off_diagnol_y =offdiagnol(y_FOM_derivative_amplitude1,y_FOM_derivative_sigma1,y_FOM_derivative_mean1, y_FOM_derivative_amplitude2,y_FOM_derivative_sigma2,y_FOM_derivative_mean2,cov_amp1_mean1,  cov_amp1_sigma1, cov_amp1_amp2, cov_amp1_sigma2, cov_amp1_mean2, cov_mean1_sigma1, cov_mean1_amp2, cov_mean1_sigma2, cov_mean1_mean2, cov_sigma1_amp2, cov_sigma1_sigma2, cov_sigma1_mean2, cov_amp2_sigma2, cov_amp2_mean2, cov_sigma2_mean2)







    return y_diagnol_terms, off_diagnol_y












