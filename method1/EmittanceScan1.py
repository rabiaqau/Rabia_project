

# #Rabia shaheen
# #16_nov_2019

# here i have all the initial paramters values giving to fit function


# Parts of script

# import libararies
from ROOT import *
import AtlasStyle
import AtlasUtils
import ROOT
import math
import numpy as np
import os.path
import pickle
import os.path
import pickle
import sys
import glob



<<<<<<< HEAD

=======
# def GetInitialParameters(RunNumber, Xscan, Early):
# if run_number==350013 and Xscan and Early:
#     A1 = 3.22
#     .
#     .
#     .
# return [A1, A2, s1, s2, m1, m2]
>>>>>>> be0a288319fb9e829fc479f331c3f7ca319a5d0a


# initial paramter function
def initialparameter(fit, Amp1, Amp2, Mean1, Mean2 ,Sigma1, Sigma2 ):

    fit.SetParameter(0,Amp1)#amplitude                                                                                                            
    fit.SetParameter(3,Amp2)#amp                                                                                                                  
    fit.SetParameter(4,Sigma2)#sigma2                                                                                                             
    fit.SetParameter(1,Mean1)#mean1                                                                                                               
    fit.SetParameter(2,Sigma1)#sigma1                                                                                                             
    fit.SetParameter(5,Mean1)#mean 




    

    # calling function Initial_par (for those runs which have different values of initial paramters)
    # < 2000
    if run_number==350013:
        if early and not Xscan:
            initialparameter(double_gauss_fit, 0.3*amp1, 0.9*amp1, m1, m1, s1, 1.3*s1)
        if not early and Xscan:
            initialparameter(double_gauss_fit, 0.711, 0.24, -m1*0.1, m1*0.5, 0.0015, 0.5)
   # <2000         
    
    if run_number==359286:
        if Xscan:
            initialparameter(double_gauss_fit, 0.6*amp1, 0.4*amp1, m1, m1, 0.9*s1, 1.3*s1)
        if not Xscan:
            initialparameter(double_gauss_fit, 0.5*amp1, 0.5*amp1, m1, m1, s1, 2*s1)
    #<2000
            
    if run_number==350880:
        if early and Xscan:
            initialparameter(double_gauss_fit, 0.5*amp1, 0.5*amp1, m1, m1, s1, 2*s1)

        if early and not Xscan:
            initialparameter(double_gauss_fit, 0.9*amp1, 0.1*amp1, m1, m1, s1, 2*s1)

        if not early and not Xscan:
            initialparameter(double_gauss_fit, amp1, amp1, m1, m1, s1, 2*s1)

    #<2000
    if run_number==350144:
        if Xscan and early:
            initialparameter(double_gauss_fit, amp1, amp1, m1, m1, s1, 2*s1)


    #
    if run_number==350220:
        if Xscan :
            initialparameter(double_gauss_fit, amp1*0.52, amp1*0.49, m1, m1, s1*1.2, 2*s1)
        if not Xscan:
            initialparameter(double_gauss_fit,2.7, 0.105, -0.0012,0.0198 ,0.0207, 0.0239)
    #

    if run_number==349327:

        if Xscan :
            initialparameter(double_gauss_fit, amp1*0.43, amp1*0.54, m1, m1*1.08, s1*0.8, 1.2*s1)
        if not Xscan:
            initialparameter(double_gauss_fit, amp1*0.43, amp1*0.54, m1, m1*1.08, s1*0.8, 1.2*s1)


    #
    if run_number==348618:
        if not Xscan:
            initialparameter(double_gauss_fit,amp1*0.13, amp1*0.86, -m1*5.4, m1*2.3 ,s1*0.8, s1*1.1)

    #
    if run_number==349111:
        if Xscan:
            initialparameter(double_gauss_fit,amp1*0.46, amp1*0.52, m1*1.1, m1*0.99 ,s1*1.23, s1*0.9)

    #
    if run_number==357451:
        if Xscan:
            initialparameter(double_gauss_fit,0.00588, 0.499, 0.00306,0.00032 ,0.007930, 0.0146)
        if not Xscan:
            initialparameter(double_gauss_fit,amp1*0.6, amp1*0.37, m1, m1 ,s1, s1*2)
    #

    if run_number==349169:
        if early and Xscan:
            initialparameter(double_gauss_fit,amp1*0.27, amp1*0.75, m1, m1 ,s1, s1*2)

        if not early and Xscan:
            initialparameter(double_gauss_fit,1.05, 0.00638,-0.00028, 0.015 ,0.0145, 0.0073)

        if early and not Xscan:
            initialparameter(double_gauss_fit,amp1*0.6, amp1*0.4, m1, m1 ,s1, s1*2)
        if not early and not Xscan:
            initialparameter(double_gauss_fit,amp1, amp1, m1, m1 ,s1, s1*2)
    #


    if run_number==349335:
        if not Xscan:
            initialparameter(double_gauss_fit,2.8, 0.175,-0.00145, -0.00665 ,0.0216, 0.015)


# 
    if run_number== 349481:
        if early and Xscan:
            
            initialparameter(double_gauss_fit,amp1*0.5, amp1*0.5, m1*0.9, m1*1.2 ,s1, s1*2)
        if not early and not Xscan:
            initialparameter(double_gauss_fit,0.28, 0.8,-0.0018, 0.00125 ,0.0168, 0.0206)

#
    if run_number==349533:
        if Xscan:
            initialparameter(double_gauss_fit,amp1*0.39, amp1*0.63, m1*1.15, m1*0.9 ,s1*1.32, s1*0.93)
        if not Xscan:
            initialparameter(double_gauss_fit,amp1, amp1, m1, m1 ,s1, s1*2)

#
    if run_number==349637:
        if early and Xscan:
            
            initialparameter(double_gauss_fit,amp1*0.5, amp1*0.4, -m1, m1 ,s1, s1*1.8)
        if not early and Xscan:
            initialparameter(double_gauss_fit,amp1*0.5, amp1*0.5, -m1, m1 ,s1, s1*2)
        if not early and not Xscan:
            initialparameter(double_gauss_fit,amp1*0.76, amp1*0.24, m1*2.9, -m1*4.5 ,s1, s1*2)

#

    if run_number==349646:
        if Xscan:
            initialparameter(double_gauss_fit,amp1*0.54, amp1*0.46, m1, m1*1.1 ,s1*0.9, s1*1.23)
        if not Xscan:
            initialparameter(double_gauss_fit,amp1*0.7, amp1*0.3, m1, m1 ,s1, s1*1.2)

#
    if run_number==349693:
        if early and Xscan:
            initialparameter(double_gauss_fit,0.907, 0.822, -0.0033, 0.0036 ,0.01165, 0.01619)
        if early and  not Xscan:
            initialparameter(double_gauss_fit,0.086, 2.4, -0.00546,-0.000138 ,0.0147, 0.0210)
        if not early and Xscan:
            initialparameter(double_gauss_fit,amp1, amp1, m1, m1 ,s1, s1)
        
#
    if run_number==349944:
        if early and Xscan:
            initialparameter(double_gauss_fit,amp1*0.83, amp1*0.17, m1, m1*1.1 ,s1*1.1, s1*0.75)
        if early and not Xscan:
            initialparameter(double_gauss_fit,amp1*0.9, amp1*0.1, m1, m1 ,s1, s1*2)
        if not early and Xscan:
            initialparameter(double_gauss_fit,amp1*0.9, amp1*0.1, m1, m1 ,s1, s1*2)
        if not early and not Xscan:
            initialparameter(double_gauss_fit,amp1*0.3, amp1*0.8, m1, m1 ,s1, s1*1.1)

# 
    if run_number==349977:
        if Xscan:
            initialparameter(double_gauss_fit,amp1*0.51, amp1*0.5, m1*0.96, m1*1.1 ,s1*0.9, s1*1.3)


#
    if run_number==350067:
        if Xscan:
            initialparameter(double_gauss_fit,amp1*0.48, amp1*0.53, m1*0.79, m1*1.3 ,s1*0.87, s1*1.12)
        if not Xscan:
            initialparameter(double_gauss_fit,amp1*0.6, amp1*0.4, m1, m1 ,s1, s1*2)

# check this scan early y and late xscan these are not responiding correctly
    if run_number==350121:
        if Xscan and early:
            initialparameter(double_gauss_fit,amp1*0.52, amp1*0.49, m1, m1 ,s1, s1*2)
        if not Xscan and early:
            initialparameter(double_gauss_fit,amp1, amp1, m1, m1 ,s1, s1*1.8)
        # if not early and Xscan:
        #     initialparameter(double_gauss_fit,amp1*0.5, amp1*0.5, m1, m1 ,s1, s1*1.28)
        if not early and not Xscan:
            initialparameter(double_gauss_fit,amp1, amp1, m1, m1 ,s1, s1*2)
##    
    if run_number==350184:
        if early and Xscan:
            initialparameter(double_gauss_fit,amp1*0.53, amp1*0.49, m1*1.25, m1*0.9 ,s1*1.2, s1*0.8)
        if early and not Xscan:
            initialparameter(double_gauss_fit,amp1, amp1, m1, m1 ,s1, s1*2)

        if not early and Xscan:
            initialparameter(double_gauss_fit,amp1*0.9, amp1*0.9, m1, m1 ,s1, s1*1.5)
##
    if run_number==350440:
        if Xscan and early:
            initialparameter(double_gauss_fit,amp1*0.7, amp1, m1, m1 ,s1, s1*1.1)
        else:
            initialparameter(double_gauss_fit,amp1, amp1, m1, m1 ,s1, s1*2)

    if run_number==350479:
        if Xscan and early:
            initialparameter(double_gauss_fit,0.473, 1.33,-0.00120, -0.00224 ,0.00966, 0.014)
        if not early and Xscan:
            initialparameter(double_gauss_fit,1.00, 0.336, 0.00016,0.29 ,0.0140, 0.047)
        else:
            initialparameter(double_gauss_fit,amp1, amp1, m1, m1 ,s1, s1*2)

    if run_number==350531:
        if early and Xscan:
            initialparameter(double_gauss_fit,amp1*0.58, amp1*0.42, m1*1.3, m1*0.6 ,s1*1.2, s1*0.9)
        if not early and Xscan:
            initialparameter(double_gauss_fit,amp1*0.5, amp1*0.08, m1, m1*8 ,s1, s1*4.7)
#
    if run_number==350682:
        if early and Xscan:
            initialparameter(double_gauss_fit,amp1*0.9, amp1*0.3, m1, m1 ,s1, s1*1.3)
        if early and not Xscan:
            initialparameter(double_gauss_fit,2.5, 0.0184, 0.000104, 0.036 ,0.0214, 0.0137)

        if not early and Xscan:
            initialparameter(double_gauss_fit,0.056, 0.65,-0.0031, 0.00103 ,0.022, 0.0138)
        if not early and not Xscan:
            initialparameter(double_gauss_fit,amp1*0.5, amp1*0.5, -0.000506, 0.312 ,0.018, 0.312)

#
    if run_number==350842:
        if early and Xscan:
            initialparameter(double_gauss_fit,amp1*0.44, amp1*0.57, m1*1.09, m1 ,s1*0.9, s1*1.3)


#
    if run_number==351062:
        if Xscan:
            initialparameter(double_gauss_fit,amp1*0.44, amp1*0.57, m1*1.09, m1 ,s1*0.9, s1*1.3)
    #
    if run_number==351223:
        if early and Xscan:

            initialparameter(double_gauss_fit,amp1*0.53, amp1*0.48, m1*0.89, m1*1.2 ,s1*0.9, s1*1.2)
        if not early and not Xscan:
            initialparameter(double_gauss_fit,0.440, 0.63,-0.000312, 0.00139 ,0.01611, 0.020)

    if run_number==351364:
        if early and Xscan:
            initialparameter(double_gauss_fit,amp1*0.41, amp1*0.6, m1*0.99, m1*1.1 ,s1*1.2, s1*0.9)
        if not early and Xscan:
            initialparameter(double_gauss_fit,amp1*0.41, amp1*0.6, m1*0.99, m1*1.1 ,s1*1.2, s1*0.9)
        if not early and not Xscan:
            initialparameter(double_gauss_fit,amp1*0.76, amp1*0.23, -m1*0.4, m1*7 ,s1*0.9, s1*1.23)

    if run_number==351550:
        if Xscan:
            initialparameter(double_gauss_fit,amp1*0.64, amp1*0.37, m1*9.9, m1*12 ,s1*0.9, s1*1.23)


    if run_number==351671:
        if Xscan and early:
            initialparameter(double_gauss_fit,amp1*0.29, amp1*0.7, m1*1.2, m1*0.99 ,s1*1.2, s1*0.9)
        if  early and not Xscan:
            initialparameter(double_gauss_fit,0.0767, 2.5, 0.0263, 0.000862 ,0.0169, 0.0206)

    if run_number==351969:
        if early and Xscan:
            initialparameter(double_gauss_fit,0.667, 0.913, -0.0024, -0.00210 ,0.0154, 0.011)

        if not early and Xscan:
            initialparameter(double_gauss_fit,0.033,0.626,-0.00294, 0.000610 ,0.0188, 0.0140)

        if not early and not Xscan:
            initialparameter(double_gauss_fit,0.155,0.725,-0.00157, 0.000150 ,0.015, 0.0199)



    if run_number==352107:
        if Xscan:
            initialparameter(double_gauss_fit,0.636, 0.937, -0.000579,-0.000184  ,0.0104, 0.01513)

        if not Xscan:
            initialparameter(double_gauss_fit,amp1*0.4, amp1*0.6, m1*1.7, 0.5*m1 ,s1*0.8, s1*1.2)
        
    if run_number==352448:
        if Xscan:
            initialparameter(double_gauss_fit,0.155,0.725,-0.00157, 0.000150 ,0.015, 0.0199)

        if not Xscan:
            initialparameter(double_gauss_fit,amp1*0.8, amp1*0.2, m1*3.2, m1*2.4 ,s1*1.1, s1*0.8)


    ### this run is very problemtic
    # if run_number==352448:
    #     if Xscan:
    #         initialparameter(double_gauss_fit,0.635, 0.00069, -0.0, 0.0 ,0.01431, 0.0069)
            

    if run_number==354396:
        if Xscan:
            initialparameter(double_gauss_fit,amp1*0.73, amp1*0.26, m1*0.3, m1*3 ,s1*0.9, s1*1.23)

        if not Xscan:
            initialparameter(double_gauss_fit,1.21, 0.35, -0.000407, 0.0027 ,0.019, 0.023)

    if run_number==355261:
        if early and Xscan:
            initialparameter(double_gauss_fit,0.626,1.16,-0.00067, -0.000948 ,0.0107, 0.0152)
        if early and not Xscan:
            initialparameter(double_gauss_fit,amp1, amp1, m1, m1 ,s1, s1*2)
        if not early and not Xscan:
            initialparameter(double_gauss_fit,amp1*0.3, amp1*0.7, m1, m1 ,s1, s1*2)

        if not early and Xscan:
            initialparameter(double_gauss_fit,amp1*0.3, amp1*0.7, m1, m1 ,s1, s1*2)
            
    if run_number==355416:
        if Xscan:
            initialparameter(double_gauss_fit,0.802,0.86,-0.00198, -0.00125 ,0.0169, 0.013)
        if not Xscan:
            initialparameter(double_gauss_fit,amp1*0.1, amp1*0.9, -m1*14, m1*0.7 ,s1*0.8, s1*1.05)

    if run_number==355468:
        if early and Xscan:
            initialparameter(double_gauss_fit,0.732, 0.92, -0.00247,-0.0026 ,0.0108, 0.0162)
        if not early and Xscan:
            initialparameter(double_gauss_fit,amp1*0.1, amp1*0.9, -m1*14, m1*0.7 ,s1*0.8, s1*1.05)

    if run_number==355529:
        if not Xscan:
            initialparameter(double_gauss_fit,0.188,2.37, -0.00235,0.001028 ,0.0160, 0.0209)
            

    if run_number==355544:
        if Xscan:
            initialparameter(double_gauss_fit,amp1*0.66, amp1*0.37, m1*1.08, m1*0.99 ,s1, s1*2)
    

    if run_number==355754:
        if early and Xscan:
            initialparameter(double_gauss_fit,amp1*0.64, amp1*0.37, m1, m1 ,s1, s1*2)

    

    if run_number==356177:
        if not Xscan:

            initialparameter(double_gauss_fit,amp1*0.43, amp1*0.58, m1, m1 ,s1*1.2, s1*0.8)


    if run_number==356250:
        if Xscan:
            initialparameter(double_gauss_fit,amp1*0.43, amp1*0.58, m1, m1 ,s1*1.2, s1*0.8)
        
    if run_number==357283:
        if early and Xscan:
            initialparameter(double_gauss_fit,amp1*0.77, amp1*0.23, m1*1.57, m1*0.099 ,s1*0.98, s1*1.4)

    if run_number==357409:
        if early and Xscan:
            initialparameter(double_gauss_fit,amp1*0.77, amp1*0.23, m1*1.57, m1*0.099 ,s1*0.98, s1*1.4)
        if early and not Xscan:
            initialparameter(double_gauss_fit,0.286, 2.039,-0.00142, 0.000374 ,0.0163, 0.0215)

        else:
            initialparameter(double_gauss_fit,amp1, amp1, m1, m1 ,s1, s1*2)

    if run_number==357500:
        if Xscan:
            initialparameter(double_gauss_fit,amp1*0.67, amp1*0.35, m1*1.12, m1*0.99 ,s1*0.8, s1*1.2)

    if run_number==357539:
        if not Xscan:
            initialparameter(double_gauss_fit,amp1, amp1, m1, m1 ,s1, s1*2)


    if run_number==357620:
        if early and Xscan:
            initialparameter(double_gauss_fit,amp1*0.38, amp1*0.6, m1*2.4, m1*0.02 ,s1*0.8, s1*1.2)
        if early and not Xscan:
            initialparameter(double_gauss_fit,amp1*0.38, amp1*0.6, m1*2.4, m1*0.02 ,s1*0.8, s1*1.2)

        if not early and Xscan:
            initialparameter(double_gauss_fit,amp1*0.79, amp1*0.2, m1, m1 ,s1, s1*2)

        if not early and not Xscan:
           
            initialparameter(double_gauss_fit,amp1*0.17, amp1*0.83, -m1*1.8, m1*0.90 ,s1*0.8, s1*2)
            

    if run_number==357713:
        if early and Xscan:
            initialparameter(double_gauss_fit,amp1*0.95, amp1*0.52, m1*0.09, m1*3 ,s1, s1*0.5)
        if early and not Xscan:
            initialparameter(double_gauss_fit,0.9*amp1,0.1*amp1,m1,m1 ,s1,2*s1)

        if not early and Xscan:
            initialparameter(double_gauss_fit,0.787,0.00567, -0.00234, 0.0205 ,0.0137, 0.0205)

        if not early and not Xscan:
            initialparameter(double_gauss_fit,amp1*0.95, amp1*0.52, m1*1.2, m1*1.8 ,s1*1.09, s1*0.5)
        
    if run_number==357887:
        if early and Xscan:
            initialparameter(double_gauss_fit,amp1*0.52, amp1*0.49, m1*0.98, m1 ,s1*0.9, s1*1.2)

        if early and not Xscan:
            initialparameter(double_gauss_fit,amp1*0.8, amp1*0.2, m1, m1 ,s1, s1*2)


        if not early and Xscan:
            initialparameter(double_gauss_fit,amp1, amp1, m1, m1 ,s1, s1*2)
        
    if run_number==358031:
        if early and Xscan:
            initialparameter(double_gauss_fit,amp1*0.84, amp1*0.18, m1, m1*10 ,s1*0.94, s1*1.5)

        if not early and Xscan:
            initialparameter(double_gauss_fit,amp1*0.84, amp1*0.18, m1, m1*10 ,s1*0.94, s1*1.5)

        if not early and not Xscan:
            initialparameter(double_gauss_fit,amp1*0.84, amp1*0.18, m1, m1*10 ,s1*0.94, s1*1.5)

    if run_number==358300:
        if early and Xscan:
            initialparameter(double_gauss_fit,0.873, 0.82,-0.00098, -0.00170 ,0.0112, 0.0156)

    if run_number==358395:
        if early and Xscan:
            initialparameter(double_gauss_fit,amp1*0.49, amp1*0.52, m1, m1 ,s1*1.2, s1*0.8)
        if not early and Xscan:
            initialparameter(double_gauss_fit,amp1, amp1*0.1, m1, m1 ,s1, s1*2)


    if run_number==358541:
        if not early and Xscan:
            initialparameter(double_gauss_fit,amp1*0.3, amp1*0.7, m1, m1 ,s1, s1*2)
        if not early and not Xscan:
            initialparameter(double_gauss_fit,amp1*0.35, amp1*0.67, -m1*0.07, m1 ,s1*0.87, s1*1.05)

    if run_number==358615:

        if early and not Xscan:
            
            initialparameter(double_gauss_fit,amp1*0.3, amp1*0.75, m1, m1 ,s1, s1*2)

        if not early and Xscan:
            initialparameter(double_gauss_fit,0.003, amp1*0.75, 0.0075, m1 ,0.003, s1*2)

        if not early and not Xscan:
            initialparameter(double_gauss_fit,amp1*0.3, amp1*0.75, m1, m1 ,s1, s1*2)


    if run_number==359058:
        if Xscan:
            initialparameter(double_gauss_fit,amp1*0.7, amp1*0.36, m1, m1*1.1 ,s1*1.3, s1*0.8)

    if run_number==359170:
        if Xscan:
            initialparameter(double_gauss_fit,amp1*0.8, amp1*0.18, m1*0.98, m1*1.4 ,s1*0.98, s1*1.6)
    

    if run_number==359191:
        if early and Xscan:

            initialparameter(double_gauss_fit,0.5, 1.01, m1, m1 ,s1, s1*2)

        if early and not Xscan:
            initialparameter(double_gauss_fit,0.5*amp1, 0.08, m1, m1 ,s1, s1*2)


        if not early and Xscan:
            initialparameter(double_gauss_fit,amp1*0.17, amp1*0.83, -m1*1.8, m1*0.90 ,s1*0.8, s1*2)

        if not early and not Xscan:
            initialparameter(double_gauss_fit,amp1*0.37, amp1*0.67, m1, m1 ,s1, s1*2)

    if run_number==359355:
        if Xscan:
            initialparameter(double_gauss_fit,0.580,1.07,-0.00187, -0.00244 ,0.0106, 0.01497)
 

    if run_number==359472:

        if Xscan:
            initialparameter(double_gauss_fit,0.580,1.07,-0.00187, -0.00244 ,0.0106, 0.01497)
        if not Xscan:
            initialparameter(double_gauss_fit,0.15,0.621,-0.00128, 0.00135 ,0.01417, 0.0197)
 

    if run_number==363400:
        if Xscan:
            initialparameter(double_gauss_fit,amp1*0.97, amp1*0.48, m1, m1*1.6 ,s1*0.8, s1*2.1)
        if not Xscan:
            initialparameter(double_gauss_fit,0.62, 1.58,-0.0024, 0.00207 ,0.013, 0.0198)

            # 1_x scan is problemeatic
    if run_number==363710:
        if early and Xscan:
            initialparameter(double_gauss_fit,amp1*0.56, amp1*0.46, m1*0.98, m1 ,s1*0.9, s1*1.3)
        if early and not Xscan:
            initialparameter(double_gauss_fit,amp1*0.9, amp1*0.1, m1, m1 ,s1, s1*2)

        if not early and not Xscan:
            initialparameter(double_gauss_fit,amp1*0.5, amp1*0.67, m1, m1 ,s1, s1*2.0)


    if run_number==363830:
        if early and Xscan:
            initialparameter(double_gauss_fit,1.5, 0.8, -0.00179,-0.15 ,0.0128, 0.176)
        if not early and not Xscan:
            initialparameter(double_gauss_fit,0.552, 0.0738,-0.00000122, 0.004905 ,0.016, 0.029)

    if run_number==363910:
        if Xscan:
            initialparameter(double_gauss_fit,0.552, 0.0738,-0.00000122, 0.004905 ,0.016, 0.029)
        if not Xscan:
            initialparameter(double_gauss_fit,0.552, 0.0738,-0.00000122, 0.004905 ,0.016, 0.029)


    if run_number==363947:
        if Xscan:
            initialparameter(double_gauss_fit,amp1*0.4, amp1*0.6, m1, m1 ,s1*1.3, s1*0.9)


    if run_number==364030:
        if early and Xscan:
            initialparameter(double_gauss_fit,0.64, 0.98,-0.000119, 0.000312 ,0.016, 0.0116)
        if early and not Xscan:
            initialparameter(double_gauss_fit,amp1*0.3, amp1*0.7, m1, m1 ,s1, s1*1.2)
        if not early and Xscan:
            initialparameter(double_gauss_fit,amp1, amp1, m1, m1 ,s1, s1*2.0)
        if not Xscan and not Xscan:
            initialparameter(double_gauss_fit,amp1*0.4, amp1*0.7, -m1, m1 ,s1, s1*2.0)


    if run_number==364076:
        if not Xscan:
            initialparameter(double_gauss_fit,amp1*0.42, amp1*0.6, m1, m1 ,s1, s1*2.0)


    if run_number==364098:
        if early and Xscan:
            initialparameter(double_gauss_fit,0.9,0.7, m1, m1*1.06 ,s1*0.9, s1*1.29)
        if early and not Xscan:
            initialparameter(double_gauss_fit,amp1*0.9, amp1*0.1, m1, m1 ,s1, s1*2.0)


        if not early and Xscan:

            initialparameter(double_gauss_fit,0.058,0.767,-0.00618 , 0.00104 ,0.0148, 0.0141)

        if not early and not Xscan:
            initialparameter(double_gauss_fit,amp1, amp1*0.6, m1, m1 ,s1, s1*2.0)


    if run_number==364292:
        if early and Xscan:
            initialparameter(double_gauss_fit,amp1*0.35, amp1*0.69, m1, m1 ,s1, s1*1.1)   
    if run_number==348354:
        if early and Xscan:
            initialparameter(double_gauss_fit,amp1*0.19, amp1*0.82, m1, m1 ,s1*1.9, s1*0.7)   


    if run_number==348251:
        if early and Xscan:
            initialparameter(double_gauss_fit,amp1*0.37, amp1*0.7, m1, -m1 ,s1*0.78, s1*0.8)
        if early and not Xscan:
            initialparameter(double_gauss_fit,amp1*0.5, amp1*0.5, -m1, m1 ,s1, s1*2.0)
        else:
            initialparameter(double_gauss_fit,amp1, amp1, m1, m1 ,s1, s1*2.0)
    if run_number==352107:
        if Xscan:
            initialparameter(double_gauss_fit,amp1*0.4, amp1*0.6, m1*1.7, m1*0.5 ,s1*0.8, s1*1.2)


    if run_number==352274:
        if Xscan:
            initialparameter(double_gauss_fit,amp1*0.41, amp1*0.56, m1*0.91, m1*1.1 ,s1*0.9, s1*1.2)   
        if not Xscan:
            initialparameter(double_gauss_fit,amp1*0.4, amp1*0.6, m1*1.7, m1*0.5 ,s1*0.8, s1*1.2)   
    if run_number==354124:
        if Xscan:
            initialparameter(double_gauss_fit,amp1*0.5, 0.82, m1, m1 ,s1*0.8, s1*2)   
        if not Xscan:
            initialparameter(double_gauss_fit,2.19, 0.198, 0.000944, -0.0056 ,0.0213, 0.0168)   
    if run_number ==355389:
        if Xscan:
            initialparameter(double_gauss_fit,amp1*0.42,amp1*0.59, m1*4.5, m1*2,s1*0.8, s1*1.2)
    

    if run_number==355529:
        if not Xscan:
            initialparameter(double_gauss_fit,amp1*0.1,amp1*0.9, m1, m1,s1, s1*2)


    if run_number==357821:
        if not Xscan:
            initialparameter(double_gauss_fit,amp1,amp1, m1, m1,s1, s1*2)

        if not Xscan and not early:
            initialparameter(double_gauss_fit,amp1*0.18,amp1*0.83, m1, m1,s1, s1*2)
    if run_number==359735:
        if early and not Xscan:
            initialparameter(double_gauss_fit,1.95,0.735,-0.00138, 0.00810,0.01969, 0.0205)

        if not early and Xscan:
            initialparameter(double_gauss_fit,amp1*0.18,0.27, m1, m1,s1, s1*2)

        if not early and not Xscan:
            initialparameter(double_gauss_fit,amp1*0.18,amp1*0.83, m1, m1,s1, s1*2)

    if run_number==355995:
        if early and Xscan:
            initialparameter(double_gauss_fit,amp1*0.5,amp1*0.48, m1, m1, s1, s1*2)

        if early and not Xscan:
            initialparameter(double_gauss_fit,amp1*0.5,amp1*0.48, m1, m1, s1, s1*2)

        if not early and not Xscan:
            initialparameter(double_gauss_fit,amp1*0.4,amp1*0.6, m1, m1,s1, s1*2)


    if run_number==349011:
        if early and Xscan:
            initialparameter(double_gauss_fit,amp1*0.6,amp1*0.4, m1, m1,s1, s1*1.3)
        if not early:
            initialparameter(double_gauss_fit,amp1*0.9,amp1*0.1, m1, m1,s1, s1*2)


    if run_number==349051:
        if early and Xscan:
            initialparameter(double_gauss_fit,amp1,amp1, m1, m1,s1, s1*2)
        if early and not Xscan:
            initialparameter(double_gauss_fit,amp1*0.2,amp1*0.8, m1, m1,s1, s1*2)
        if not early and Xscan:

            initialparameter(double_gauss_fit,0.112,1.60, 0.051,-0.000196,0.075, 0.0150)
        if not early and not Xscan:
            initialparameter(double_gauss_fit,amp1,amp1, m1, m1,s1, s1*2)

    if run_number==349268:
        if not early and not Xscan:
            initialparameter(double_gauss_fit,amp1*0.3,amp1*0.7, m1, m1,s1, s1*1.0011)


    if run_number==350160:
        if not Xscan:
            initialparameter(double_gauss_fit,amp1*0.2,amp1*0.9, -m1, m1,s1, s1*1.2)

    if run_number==358516:
        if not Xscan:
            initialparameter(double_gauss_fit,amp1,amp1, m1, m1,s1, s1*2)

    if run_number==358656:
        if not Xscan:
            initialparameter(double_gauss_fit,amp1,amp1, m1, m1,s1, s1*2)
            

    if run_number==359124:
        if Xscan:
            initialparameter(double_gauss_fit,amp1*0.3,amp1*0.7, m1, m1,s1, s1*2)
        if not Xscan:
            initialparameter(double_gauss_fit,amp1,amp1, m1, m1,s1, s1*2)


    if run_number==359310:
        if not Xscan:
            initialparameter(double_gauss_fit,0.4,0.29, -0.000523, 0.001521,0.0156, 0.020)
    if run_number==359717:# problematic scan
        if not Xscan:
            initialparameter(double_gauss_fit,amp1*0.4,amp1*0.7, m1, m1,s1, s1*2)
    if run_number==359823:
        if early and not Xscan:
            initialparameter(double_gauss_fit,amp1*0.4,amp1*0.7, m1, m1,s1, s1*2)

        if not early and Xscan:
            
            initialparameter(double_gauss_fit,amp1*0.3,amp1*0.7, m1, m1,2*s1, s1)

        if not early and not Xscan:
            initialparameter(double_gauss_fit,amp1*0.3,amp1*0.7, m1, m1,2*s1, s1)

    if run_number==359918:
        if Xscan:
            initialparameter(double_gauss_fit,0.9,0.7, -0.0021,-0.00189,0.0159, 0.0110)
        if not Xscan:
            initialparameter(double_gauss_fit,2.5,0.306,0.000826,0.00624,0.0255,0.0206 )


    if run_number==360244:
        if Xscan:
            
            initialparameter(double_gauss_fit,amp1*0.9,amp1*0.1, m1, m1,s1, s1*2)

        if not Xscan:
            initialparameter(double_gauss_fit,amp1*0.1,amp1, m1, m1,s1*2, s1)
        

    if run_number==360414:
        if Xscan:
            initialparameter(double_gauss_fit,amp1*0.6,amp1*0.4, m1, m1,s1*2, s1)
        if not Xscan:
            initialparameter(double_gauss_fit,amp1*0.6,amp1*0.4, m1, m1,s1*2, s1)

    if run_number==361738:

        if early and not Xscan:
            initialparameter(double_gauss_fit,2.3,0.01, m1, m1,s1, s1*2)


        if not early and Xscan:
            initialparameter(double_gauss_fit,amp1*0.6,amp1*0.6, m1, m1,s1, s1*1.2)

        if not early and not Xscan:
            initialparameter(double_gauss_fit,amp1*0.6,amp1*0.6, m1, m1,s1, s1*1.2)


    if run_number==362445:

        if not early and Xscan:
            initialparameter(double_gauss_fit,amp1*0.5,amp1*0.5, m1, m1,s1, s1*1.2)

        if early and not Xscan:
            initialparameter(double_gauss_fit,amp1*0.6,amp1*0.5, m1, m1,s1, s1*2)


        if not early and not Xscan:
            initialparameter(double_gauss_fit,amp1*0.6,amp1*0.0, m1, m1,s1, s1*1.2)


    if run_number==362661:
        if early and Xscan:
            initialparameter(double_gauss_fit,amp1*0.5,amp1*0.5, m1, m1,s1, s1*1.2)
        if early and not Xscan:
            initialparameter(double_gauss_fit,amp1*0.5,amp1*0.5, m1, m1,s1, s1*1.2)

        if not early and Xscan:
            initialparameter(double_gauss_fit,0.934,0.036,-0.001311,0.00604,0.0146, 0.0276)

        if not early and not Xscan:
            initialparameter(double_gauss_fit,0.934,0.036,-0.001311,0.00604,0.0146, 0.0276)


    if run_number==363033:

        initialparameter(double_gauss_fit,amp1,amp1, m1, m1,s1, s1*2)
            
    if run_number==363198:
        if early and Xscan:
            initialparameter(double_gauss_fit,0.85,0.74, m1, m1,s1, s1*1.2)
        if early and not Xscan:
            initialparameter(double_gauss_fit,0.85,0.74, m1, m1,s1, s1*1.2)

        if not early and Xscan:
            initialparameter(double_gauss_fit,0.47,0.0136, m1, m1,s1, s1*2)

        if not Xscan and not early:
            initialparameter(double_gauss_fit,0.213,0.381, m1, m1,s1, s1*2)


    if run_number==363979:#problematic
        if Xscan:
            initialparameter(double_gauss_fit,0.26,0.37, 0.00076,0.04,0.0210, 0.176)

        if not Xscan:
            initialparameter(double_gauss_fit,0.26,0.37, 0.00076,0.04,0.0210, 0.176)



    if run_number==364214:
#        if not early:
        if Xscan:
            initialparameter(double_gauss_fit,0.26,0.37, 0.00076,0.04,0.0210, 0.176)

        if not Xscan:
            initialparameter(double_gauss_fit,amp1*0.6,amp1*0.4, m1, m1,s1, s1*1.8)


    if run_number==357772:
        if Xscan:
            initialparameter(double_gauss_fit,0.46,0.869,-0.00361,-0.00331,0.0157, 0.0125)
        if not Xscan:
            initialparameter(double_gauss_fit,0.46,0.869,-0.00361,-0.00331,0.0157, 0.0125)
           
    if run_number==357679:
        if Xscan:
            initialparameter(double_gauss_fit,0.26,0.37, 0.00076,0.04,0.0210, 0.176)
        if not Xscan:
            initialparameter(double_gauss_fit,0.68,3.5, m1, m1,s1, s1*2)


    if run_number==351628:
        if Xscan:
            initialparameter(double_gauss_fit,0.68,3.5, m1, m1,s1, s1*2)
        if not Xscan:
            initialparameter(double_gauss_fit,2.0,0.0152, m1, m1,s1, s1*2)
        
    if run_number==350361:
        if early and Xscan:
            initialparameter(double_gauss_fit,amp1*0.6,amp1*0.4, m1, m1,s1, s1*1.2)
        if early and not Xscan:
            initialparameter(double_gauss_fit,0.48,2.0, m1, m1,s1, s1*1.2)

            
    if run_number==354359:
        if Xscan:
            initialparameter(double_gauss_fit,0.865,0.869,-0.00361,-0.00331,0.0157, 0.0125)
        if not Xscan:
            initialparameter(double_gauss_fit,amp1,amp1, m1, m1,s1, s1*1.2)

    if run_number==349592:
        if Xscan:
            initialparameter(double_gauss_fit,1.09,0.765, m1, m1,s1, s1)
        if not Xscan:
            initialparameter(double_gauss_fit,0.26,0.37, 0.00076,0.04,0.0210, 0.176)

    if run_number==349114:

        if not early and Xscan:
            initialparameter(double_gauss_fit,0.566,0.59,-0.00133,0.000124,0.0158, 0.152)

    if run_number==349309:
        if Xscan:
            initialparameter(double_gauss_fit,0.6,2.5, m1, m1,s1, s1*2)

        if not Xscan:
            initialparameter(double_gauss_fit,0.6,2.5, m1, m1,s1, s1*2)


    if run_number==349033:
        if not early and Xscan:
            initialparameter(double_gauss_fit,1.16,0.1317, m1, m1,s1, s1*2)
        if not early and not Xscan:
            initialparameter(double_gauss_fit,amp1,amp1, m1, m1,s1, s1*2)

    if run_number==349014:
        if not Xscan:
            initialparameter(double_gauss_fit,amp1,amp1, m1, m1,s1, s1*2)



















