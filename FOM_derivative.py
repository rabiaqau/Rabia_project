import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
import math
from sympy import exp
from sympy import Function, diff, Eq,Derivative


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






def FOM_x(sigma_y, n1n2,error_a1x, error_a2x,error_a1a2x):



    derivative_a1x = (31.3249 * sigma_y) / n1n2


    derivative_a2x = (31.3249 * sigma_y) / n1n2



    relative_error_a1x = (derivative_a1x * error_a1x)**2



    relative_error_a2x = (derivative_a2x * error_a2x)**2


    
    off_dignol_error_a1a2x = 2 * derivative_a1x * derivative_a2x * error_a1a2x



    x_terms = relative_error_a1x + relative_error_a2x + off_dignol_error_a1a2x


    
    return x_terms





def FOM_y(n1n2, a1x, a2x, der_sigma_y_a1y, der_sigma_y_a2y, der_sigma_y_s1y, der_sigma_y_s2y, der_sigma_y_m1y, der_sigma_y_m2y, error_a1y, error_a2y, error_s1y, error_s2y, error_m1y, error_m2y):


    derivative_a1y= (31.3249 * (a1x + a2x) * der_sigma_y_a1y ) / n1n2


    derivative_a2y= (31.3249 * (a1x + a2x) * der_sigma_y_a2y ) / n1n2


    derivative_s1y= (31.3249 * (a1x + a2x) * der_sigma_y_s1y ) / n1n2


    derivative_s2y= (31.3249 * (a1x + a2x) * der_sigma_y_s2y ) / n1n2


    derivative_m1y= (31.3249 * (a1x + a2x) * der_sigma_y_m1y ) / n1n2


    derivative_m2y= (31.3249 * (a1x + a2x) * der_sigma_y_m2y ) / n1n2



    diagnol_terms = ( derivative_a1y * error_a1y )**2 +( derivative_a2y * error_a2y )**2 + ( derivative_s1y * error_s1y )**2 + ( derivative_s2y * error_s2y )**2 + ( derivative_m1y * error_m1y )**2 + ( derivative_m2y * error_m2y )**2
    

    return derivative_a1y, derivative_a2y, derivative_s1y, derivative_s2y, derivative_m1y, derivative_m2y, diagnol_terms
    












# constant is (2*pi / n1n2 * 1e22 * 8 * 1e-24)


#y,n1n2, amplitude1_x, amplitude2_x, amplitude1_y, amplitude2_y, sigma1_y, mean1_y, amplitude2_y, sigma2_y,mean2_y = sym.symbols(' y n1n2 amplitude1_x amplitude2_x amplitude1_y amplitude2_y sigma1_y mean1_y amplitude2_y sigma2_y mean2_y',real=True)


#FOM = ( 12.5 * (amplitude1_x + amplitude2_x) * (amplitude1_y + amplitude2_y)  ) / n1n2 / ( amplitude1_y * (1.0 / ( math.sqrt(2 * math.pi) * sigma1_y ) ) * exp( - (1.0 / 2.0) * ( (y - mean1_y)**2 / ( sigma1_y **2) ) ) +  amplitude2_y * (1.0 / (  math.sqrt(2 * math.pi) * sigma2_y ) ) * exp( - (1.0 / 2.0) * ( (y - mean2_y)**2 / ( sigma2_y **2 ) ) )  )


#y,sigma_y, n1n2, amplitude1_x, amplitude2_x, amplitude1_y, sigma1_y, mean1_y, amplitude2_y, sigma2_y,mean2_y = sym.symbols(' y sigma_y n1n2 amplitude1_x amplitude2_x amplitude1_y sigma1_y mean1_y amplitude2_y sigma2_y mean2_y',real=True)
#f =Function('f')
#sigma_y = f(amplitude1_y, sigma1_y, mean1_y, amplitude2_y, sigma2_y,mean2_y)

# y,sigma_y, n1n2, amplitude1_x, amplitude2_x, amplitude1_y, sigma1_y, mean1_y, amplitude2_y, sigma2_y,mean2_y = sym.symbols(' y sigma_y n1n2 amplitude1_x amplitude2_x amplitude1_y sigma1_y mean1_y amplitude2_y sigma2_y mean2_y',real=True)

# def sigma_y(amplitude1_y, sigma1_y, mean1_y, amplitude2_y, sigma2_y, mean2_y):

#     sigma_y = (amplitude1_y + amplitude2_y)   /  amplitude1_y * (1.0 / ( math.sqrt(2 * math.pi) * sigma1_y ) ) * exp( - (1.0 / 2.0) * ( (y - mean1_y)**2 / ( sigma1_y **2) ) ) +  amplitude2_y * (1.0 / (  math.sqrt(2 * math.pi) * sigma2_y ) ) * exp( - (1.0 / 2.0) * ( (y - mean2_y)**2 / ( sigma2_y **2 ) ) ) 


    
#     return  31.3249 * sigma_y


# y, n1n2, amplitude1_x, amplitude2_x, amplitude1_y, sigma1_y, mean1_y, amplitude2_y, sigma2_y,mean2_y = sym.symbols(' y n1n2 amplitude1_x amplitude2_x amplitude1_y sigma1_y mean1_y amplitude2_y sigma2_y mean2_y',real=True)

# test = diff(sigma_y(amplitude1_y, sigma1_y, mean1_y, amplitude2_y, sigma2_y, mean2_y),amplitude1_y)



# print test.evalf(subs={ amplitude1_y:1 ,sigma1_y:1 , mean1_y:1 ,amplitude2_y:1 , sigma2_y:1 ,mean2_y: 1 , y:1}) 

#FOM = ( 31.3249 * (amplitude1_x + amplitude2_x) * sigma_y) / n1n2

#print derivative_FOM(sigma_y,amplitude1_y)
#print sigma_y.diff( amplitude1_y),"(amplitude1_y, sigma1_y, mean1_y, amplitude2_y, sigma2_y, mean2_y)"
   
#print Derivative(FOM,amplitude1_y)


#derivative_FOM_amplitude1_y=derivative_FOM(FOM,amplitude1_y)

#print derivative_FOM_amplitude1_y,"derivative_FOM_amplitude1_y"


#  defining vaiables




# derivative_FOM_amplitude1_x=derivative_FOM(FOM,amplitude1_x)
# derivative_FOM_amplitude2_x=derivative_FOM(FOM,amplitude2_x)



# #derivative_FOM_amplitude1_x=derivative_FOM(FOM,amplitude1_x)


# def values_derivative_X(amp1,sigma1,mean1,amp2,sigma2,mean2,position,current,x_error_on_A1,x_error_on_A2,cov_A1_A2):


#     derivative_FOM_amplitude1_x=derivative_FOM(FOM,amplitude1_x)
#     derivative_FOM_amplitude2_x=derivative_FOM(FOM,amplitude2_x)


#     x_derivtaive_FOM_A1 = derivative_FOM_amplitude1_x.evalf(subs={ amplitude1_y:amp1 ,sigma1_y:sigma1 , mean1_y:mean1 ,amplitude2_y:amp2 , sigma2_y: sigma2 ,mean2_y: mean2 , y:position ,n1n2:current})


#     x_derivtaive_FOM_A2 = derivative_FOM_amplitude2_x.evalf(subs={ amplitude1_y:amp1 ,sigma1_y:sigma1 , mean1_y:mean1 ,amplitude2_y:amp2 , sigma2_y: sigma2 ,mean2_y: mean2 , y:position ,n1n2:current})



#     dignol_terms =( x_derivtaive_FOM_A1 * x_error_on_A1) **2 + ( x_derivtaive_FOM_A2 * x_error_on_A2) **2
    

#     off_dignol =2 * (x_derivtaive_FOM_A1)* (x_derivtaive_FOM_A2) *cov_A1_A2

#     X_terms = dignol_terms + off_dignol


#     return X_terms



# # def Set_FOM_Values:

# #     derivative_FOM_amplitude1_x=derivative_FOM(FOM,amplitude1_x)
# #     derivative_FOM_amplitude2_x=derivative_FOM(FOM,amplitude2_x)




# def values_derivative_Y(amp1_y,sig1_y,y_mea1,amp2_y,y_sig2,y_mea2,amp1_x,amp2_x,position,current,error_amp1_y, error_amp2_y, error_sig1_y, error_sig2_y, error_mea1_y, error_mea2_y, cov_amp1_mean1,  cov_amp1_sigma1, cov_amp1_amp2, cov_amp1_sigma2, cov_amp1_mean2, cov_mean1_sigma1, cov_mean1_amp2, cov_mean1_sigma2, cov_mean1_mean2, cov_sigma1_amp2, cov_sigma1_sigma2, cov_sigma1_mean2, cov_amp2_sigma2, cov_amp2_mean2, cov_sigma2_mean2):
  
#     derivative_FOM_amplitude1_y=derivative_FOM(FOM,amplitude1_y)
#     derivative_FOM_amplitude2_y=derivative_FOM(FOM,amplitude2_y)
#     derivative_FOM_sigma1_y=derivative_FOM(FOM,sigma1_y)
#     derivative_FOM_sigma2_y=derivative_FOM(FOM,sigma2_y)
#     derivative_FOM_mean1_y=derivative_FOM(FOM,mean1_y)
#     derivative_FOM_mean2_y=derivative_FOM(FOM,mean2_y)
    

#     y_FOM_derivative_amplitude1= derivative_FOM_amplitude1_y.evalf(subs={ amplitude1_x:amp1_x, amplitude2_x:amp2_x ,amplitude1_y:amp1_y ,sigma1_y:sig1_y , mean1_y:y_mea1 ,amplitude2_y:amp2_y ,sigma2_y: y_sig2 ,mean2_y: y_mea2 , y:position ,n1n2:current})



#     y_FOM_derivative_amplitude2=derivative_FOM_amplitude2_y.evalf(subs={ amplitude1_x:amp1_x, amplitude2_x:amp2_x ,amplitude1_y:amp1_y ,sigma1_y:sig1_y , mean1_y:y_mea1 ,amplitude2_y:amp2_y ,sigma2_y: y_sig2 ,mean2_y: y_mea2 , y:position ,n1n2:current})

    
#     y_FOM_derivative_sigma1=derivative_FOM_sigma1_y.evalf(subs={ amplitude1_x:amp1_x, amplitude2_x:amp2_x ,amplitude1_y:amp1_y ,sigma1_y:sig1_y , mean1_y:y_mea1 ,amplitude2_y:amp2_y ,sigma2_y: y_sig2 ,mean2_y: y_mea2 , y:position ,n1n2:current})

#     y_FOM_derivative_sigma2=derivative_FOM_sigma2_y.evalf(subs={ amplitude1_x:amp1_x, amplitude2_x:amp2_x ,amplitude1_y:amp1_y ,sigma1_y:sig1_y , mean1_y:y_mea1 ,amplitude2_y:amp2_y ,sigma2_y: y_sig2 ,mean2_y: y_mea2 , y:position ,n1n2:current})


#     y_FOM_derivative_mean1=derivative_FOM_mean1_y.evalf(subs={ amplitude1_x:amp1_x, amplitude2_x:amp2_x ,amplitude1_y:amp1_y ,sigma1_y:sig1_y , mean1_y:y_mea1 ,amplitude2_y:amp2_y ,sigma2_y: y_sig2 ,mean2_y: y_mea2 , y:position ,n1n2:current})


#     y_FOM_derivative_mean2=derivative_FOM_mean2_y.evalf(subs={ amplitude1_x:amp1_x, amplitude2_x:amp2_x ,amplitude1_y:amp1_y ,sigma1_y:sig1_y , mean1_y:y_mea1 ,amplitude2_y:amp2_y ,sigma2_y: y_sig2 ,mean2_y: y_mea2 , y:position ,n1n2:current})


#     y_diagnol_terms=(y_FOM_derivative_amplitude1 *error_amp1_y)**2 + (y_FOM_derivative_amplitude2*error_amp2_y)**2 +(y_FOM_derivative_sigma1*error_sig1_y)**2+(y_FOM_derivative_sigma2 *error_sig2_y)**2 +(y_FOM_derivative_mean1 *error_mea1_y)**2 + (y_FOM_derivative_mean2*error_mea2_y)**2




#     off_diagnol_y =offdiagnol(y_FOM_derivative_amplitude1,y_FOM_derivative_sigma1,y_FOM_derivative_mean1, y_FOM_derivative_amplitude2,y_FOM_derivative_sigma2,y_FOM_derivative_mean2,cov_amp1_mean1,  cov_amp1_sigma1, cov_amp1_amp2, cov_amp1_sigma2, cov_amp1_mean2, cov_mean1_sigma1, cov_mean1_amp2, cov_mean1_sigma2, cov_mean1_mean2, cov_sigma1_amp2, cov_sigma1_sigma2, cov_sigma1_mean2, cov_amp2_sigma2, cov_amp2_mean2, cov_sigma2_mean2)







#     return y_diagnol_terms, off_diagnol_y












