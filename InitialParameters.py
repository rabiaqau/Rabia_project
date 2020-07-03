

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




def initialparameter(fit, Amp1, Amp2, Mean1, Mean2 ,Sigma1, Sigma2 ):

    fit.SetParameter(0,Amp1)#amplitude                                                                                            
    fit.SetParameter(3,Amp2)#amp                                                                                                  
    fit.SetParameter(4,Sigma2)#sigma2                                                                                             
    fit.SetParameter(1,Mean1)#mean1                                                                                               
    fit.SetParameter(2,Sigma1)#sigma1                                                                                             
    fit.SetParameter(5,Mean1)#mean                                                                                                


def SetInitialParameters(Fit,RunNumber,xscan,Early,amplitude,sigma,mean):

    if RunNumber==350531:


        if Early and xscan:
            initialparameter(Fit,amplitude*0.58, amplitude*0.42, mean*1.3, mean*0.6 ,sigma*1.2, sigma*0.9)

        if not Early and xscan:
            initialparameter(Fit,0.94,0.0796, mean, mean*10 ,0.014,0.066)





    if RunNumber==350479:
    
        if xscan and Early:
            initialparameter(Fit,0.473, 1.33,-0.00120, -0.00224 ,0.00966, 0.014)
        if not Early and xscan:
            initialparameter(Fit,amplitude*0.98,amplitude*0.043,mean,mean ,0.0140, 0.047)
        else:
            initialparameter(Fit,amplitude, amplitude, mean, mean ,sigma, sigma*2)




    if RunNumber==350184:

        if Early and xscan:
            initialparameter(Fit,amplitude*0.53, amplitude*0.49, mean*1.25, mean*0.9 ,sigma*1.2, sigma*0.8)

        if Early and not xscan:
            initialparameter(Fit,amplitude, amplitude, mean, mean ,sigma, sigma*2)


        if not Early and xscan:
            initialparameter(Fit,amplitude*0.93, amplitude*0.1 , mean, -mean ,sigma, sigma*1.6)




    if RunNumber==350220:
        if xscan :
            initialparameter(Fit, amplitude*0.52, amplitude*0.49, mean, mean, sigma*1.2, 2*sigma)
        if not xscan:

            initialparameter(Fit, amplitude*0.52, 0.105, mean, 0.198, sigma*1.2, 0.0239)










    if RunNumber==350013:
        if Early and not xscan:
            initialparameter(Fit, 0.3*amplitude, 0.9*amplitude, mean, mean, sigma, 1.3*sigma)
        if not Early and xscan:
            initialparameter(Fit, 0.711, 0.24, -mean*0.1, mean*0.5, 0.0015, 0.5)



#    # <2000         
    
    if RunNumber==359286:
        if xscan:
            initialparameter(Fit, 0.6*amplitude, 0.4*amplitude, mean, mean, 0.9*sigma, 1.3*sigma)
        if not xscan:
            initialparameter(Fit, 0.5*amplitude, 0.5*amplitude, mean, mean, sigma, 2*sigma)
    #<2000
            
    if RunNumber==350880:
        if Early and xscan:
            initialparameter(Fit, 0.5*amplitude, 0.5*amplitude, mean, mean, sigma, 2*sigma)

        if Early and not xscan:
            initialparameter(Fit, 0.9*amplitude, 0.1*amplitude, mean, mean, sigma, 2*sigma)


        if not Early and xscan:

            initialparameter(Fit,amplitude ,amplitude*0.3, mean, -mean, sigma, 2*sigma)

        if not Early and not xscan:
            initialparameter(Fit, amplitude, amplitude, mean, mean, sigma, 2*sigma)

    #<2000
    if RunNumber==350144:
        if xscan and Early:
            initialparameter(Fit, amplitude, amplitude, mean, mean, sigma, 2*sigma)


    #
            #initialparameter(Fit,2.7, 0.105, -0.0012,0.0198 ,0.0207, 0.0239)
            #initialparameter(Fit,2.7, 0.105, -0.0012,0.0198 ,0.0207, 0.0239)
    #

    if RunNumber==349327:

        if xscan :
            initialparameter(Fit, amplitude*0.43, amplitude*0.54, mean, mean*1.08, sigma*0.8, 1.2*sigma)
        if not xscan:
            initialparameter(Fit, amplitude*0.43, amplitude*0.54, mean, mean*1.08, sigma*0.8, 1.2*sigma)


    #
    if RunNumber==348618:
        if not xscan:
            initialparameter(Fit,amplitude*0.13, amplitude*0.86, -mean*5.4, mean*2.3 ,sigma*0.8, sigma*1.1)

    #
    if RunNumber==349111:
        if xscan:
            initialparameter(Fit,amplitude*0.46, amplitude*0.52, mean*1.1, mean*0.99 ,sigma*1.23, sigma*0.9)

    #
    if RunNumber==357451:
        if xscan:
            initialparameter(Fit,0.00588, 0.499, 0.00306,0.00032 ,0.007930, 0.0146)
        if not xscan:
            initialparameter(Fit,amplitude*0.6, amplitude*0.37, mean, mean ,sigma, sigma*2)
    #

    if RunNumber==349169:
        if Early and xscan:
            initialparameter(Fit,amplitude*0.27, amplitude*0.75, mean, mean ,sigma, sigma*2)

        if not Early and xscan:
            initialparameter(Fit,1.05, 0.00638,-0.00028, 0.015 ,0.0145, 0.0073)

        if Early and not xscan:
            initialparameter(Fit,amplitude*0.6, amplitude*0.4, mean, mean ,sigma, sigma*2)
        if not Early and not xscan:
            initialparameter(Fit,amplitude, amplitude, mean, mean ,sigma, sigma*2)
    #


    if RunNumber==349335:
        if not xscan:
            initialparameter(Fit,2.8, 0.175,-0.00145, -0.00665 ,0.0216, 0.015)


# 
    if RunNumber== 349481:
        if Early and xscan:
            
            initialparameter(Fit,amplitude*0.5, amplitude*0.5, mean*0.9, mean*1.2 ,sigma, sigma*2)
        if not Early and not xscan:
            initialparameter(Fit,0.28, 0.8,-0.0018, 0.00125 ,0.0168, 0.0206)

#
    if RunNumber==349533:
        if xscan:
            initialparameter(Fit,amplitude*0.39, amplitude*0.63, mean*1.15, mean*0.9 ,sigma*1.32, sigma*0.93)
        if not xscan:
            initialparameter(Fit,amplitude, amplitude, mean, mean ,sigma, sigma*2)

#
    if RunNumber==349637:
        if Early and xscan:
            
            initialparameter(Fit,amplitude*0.5, amplitude*0.4, -mean, mean ,sigma, sigma*1.8)
        if not Early and xscan:
            initialparameter(Fit,amplitude*0.5, amplitude*0.5, -mean, mean ,sigma, sigma*2)
        if not Early and not xscan:
            initialparameter(Fit,amplitude*0.76, amplitude*0.24, mean*2.9, -mean*4.5 ,sigma, sigma*2)

#

    if RunNumber==349646:
        if xscan:
            initialparameter(Fit,amplitude*0.54, amplitude*0.46, mean, mean*1.1 ,sigma*0.9, sigma*1.23)
        if not xscan:
            initialparameter(Fit,amplitude*0.7, amplitude*0.3, mean, mean ,sigma, sigma*1.2)

#
    if RunNumber==349693:
        if Early and xscan:
            initialparameter(Fit,0.907, 0.822, -0.0033, 0.0036 ,0.01165, 0.01619)
        if Early and  not xscan:
            initialparameter(Fit,0.086, 2.4, -0.00546,-0.000138 ,0.0147, 0.0210)
        if not Early and xscan:
            initialparameter(Fit,amplitude, amplitude, mean, mean ,sigma, sigma)
        
#
    if RunNumber==349944:
        if Early and xscan:
            initialparameter(Fit,amplitude*0.83, amplitude*0.17, mean, mean*1.1 ,sigma*1.1, sigma*0.75)
        if Early and not xscan:
            initialparameter(Fit,amplitude*0.9, amplitude*0.1, mean, mean ,sigma, sigma*2)
        if not Early and xscan:
            initialparameter(Fit,amplitude*0.9, amplitude*0.1, mean, mean ,sigma, sigma*2)
        if not Early and not xscan:
            initialparameter(Fit,amplitude*0.3, amplitude*0.8, mean, mean ,sigma, sigma*1.1)

# 
    if RunNumber==349977:
        if xscan:
            initialparameter(Fit,amplitude*0.51, amplitude*0.5, mean*0.96, mean*1.1 ,sigma*0.9, sigma*1.3)


#
    if RunNumber==350067:
        if xscan:
            initialparameter(Fit,amplitude*0.48, amplitude*0.53, mean*0.79, mean*1.3 ,sigma*0.87, sigma*1.12)
        if not xscan:
            initialparameter(Fit,amplitude*0.6, amplitude*0.4, mean, mean ,sigma, sigma*2)

# check this scan Early y and late xscan these are not responiding correctly
    if RunNumber==350121:


        if xscan and Early:
            initialparameter(Fit,amplitude*0.52, amplitude*0.49, mean, mean ,sigma, sigma*2)
        if not xscan and Early:
            initialparameter(Fit,amplitude, amplitude, mean, mean ,sigma, sigma*1.8)


        # if not Early and xscan:
        #     initialparameter(Fit,amplitude*0.5, amplitude*0.5, mean, mean ,sigma, sigma*1.28)


        if not Early and not xscan:
            initialparameter(Fit,amplitude, amplitude, mean, mean ,sigma, sigma*2)
   





##
    if RunNumber==350440:
        if xscan and Early:
            initialparameter(Fit,amplitude*0.7, amplitude, mean, mean ,sigma, sigma*1.1)
        else:
            initialparameter(Fit,amplitude, amplitude, mean, mean ,sigma, sigma*2)







#
    if RunNumber==350682:
        if Early and xscan:
            initialparameter(Fit,amplitude*0.9, amplitude*0.3, mean, mean ,sigma, sigma*1.3)
        if Early and not xscan:
            initialparameter(Fit,2.5, 0.0184, 0.000104, 0.036 ,0.0214, 0.0137)

        if not Early and xscan:
            initialparameter(Fit,0.056, 0.65,-0.0031, 0.00103 ,0.022, 0.0138)
        if not Early and not xscan:
            initialparameter(Fit,amplitude*0.5, amplitude*0.5, -0.000506, 0.312 ,0.018, 0.312)

#
    if RunNumber==350842:
        if Early and xscan:
            initialparameter(Fit,amplitude*0.44, amplitude*0.57, mean*1.09, mean ,sigma*0.9, sigma*1.3)


#
    if RunNumber==351062:
        if xscan:
            initialparameter(Fit,amplitude*0.44, amplitude*0.57, mean*1.09, mean ,sigma*0.9, sigma*1.3)
    #
    if RunNumber==351223:
        if Early and xscan:

            initialparameter(Fit,amplitude*0.53, amplitude*0.48, mean*0.89, mean*1.2 ,sigma*0.9, sigma*1.2)
        if not Early and not xscan:
            initialparameter(Fit,0.440, 0.63,-0.000312, 0.00139 ,0.01611, 0.020)

    if RunNumber==351364:
        if Early and xscan:
            initialparameter(Fit,amplitude*0.41, amplitude*0.6, mean*0.99, mean*1.1 ,sigma*1.2, sigma*0.9)
        if not Early and xscan:
            initialparameter(Fit,amplitude*0.41, amplitude*0.6, mean*0.99, mean*1.1 ,sigma*1.2, sigma*0.9)
        if not Early and not xscan:
            initialparameter(Fit,amplitude*0.76, amplitude*0.23, -mean*0.4, mean*7 ,sigma*0.9, sigma*1.23)

    if RunNumber==351550:
        if xscan:
            initialparameter(Fit,amplitude*0.64, amplitude*0.37, mean*9.9, mean*12 ,sigma*0.9, sigma*1.23)


    if RunNumber==351671:
        if xscan and Early:
            initialparameter(Fit,amplitude*0.29, amplitude*0.7, mean*1.2, mean*0.99 ,sigma*1.2, sigma*0.9)
        if  Early and not xscan:
            initialparameter(Fit,0.0767, 2.5, 0.0263, 0.000862 ,0.0169, 0.0206)

    if RunNumber==351969:
        if Early and xscan:
            initialparameter(Fit,0.667, 0.913, -0.0024, -0.00210 ,0.0154, 0.011)

        if not Early and xscan:
            initialparameter(Fit,0.033,0.626,-0.00294, 0.000610 ,0.0188, 0.0140)

        if not Early and not xscan:
            initialparameter(Fit,0.155,0.725,-0.00157, 0.000150 ,0.015, 0.0199)



    if RunNumber==352107:
        if xscan:
            initialparameter(Fit,0.636, 0.937, -0.000579,-0.000184  ,0.0104, 0.01513)

        if not xscan:
            initialparameter(Fit,amplitude*0.4, amplitude*0.6, mean*1.7, 0.5*mean ,sigma*0.8, sigma*1.2)
        
    if RunNumber==352448:
        if xscan:
            initialparameter(Fit,0.155,0.725,-0.00157, 0.000150 ,0.015, 0.0199)

        if not xscan:
            initialparameter(Fit,amplitude*0.8, amplitude*0.2, mean*3.2, mean*2.4 ,sigma*1.1, sigma*0.8)


    ### this run is very problemtic
    # if RunNumber==352448:
    #     if xscan:
    #         initialparameter(Fit,0.635, 0.00069, -0.0, 0.0 ,0.01431, 0.0069)
            

    if RunNumber==354396:
        if xscan:
            initialparameter(Fit,amplitude*0.73, amplitude*0.26, mean*0.3, mean*3 ,sigma*0.9, sigma*1.23)

        if not xscan:
            initialparameter(Fit,1.21, 0.35, -0.000407, 0.0027 ,0.019, 0.023)

    if RunNumber==355261:
        if Early and xscan:
            initialparameter(Fit,0.626,1.16,-0.00067, -0.000948 ,0.0107, 0.0152)
        if Early and not xscan:
            initialparameter(Fit,amplitude, amplitude, mean, mean ,sigma, sigma*2)
        if not Early and not xscan:
            initialparameter(Fit,amplitude*0.3, amplitude*0.7, mean, mean ,sigma, sigma*2)

        if not Early and xscan:
            initialparameter(Fit,amplitude*0.3, amplitude*0.7, mean, mean ,sigma, sigma*2)
            
    if RunNumber==355416:
        if xscan:
            initialparameter(Fit,0.802,0.86,-0.00198, -0.00125 ,0.0169, 0.013)
        if not xscan:
            initialparameter(Fit,amplitude*0.1, amplitude*0.9, -mean*14, mean*0.7 ,sigma*0.8, sigma*1.05)

    if RunNumber==355468:
        if Early and xscan:
            initialparameter(Fit,0.732, 0.92, -0.00247,-0.0026 ,0.0108, 0.0162)
        if not Early and xscan:
            initialparameter(Fit,amplitude*0.1, amplitude*0.9, -mean*14, mean*0.7 ,sigma*0.8, sigma*1.05)

    if RunNumber==355529:
        if not xscan:
            initialparameter(Fit,0.188,2.37, -0.00235,0.001028 ,0.0160, 0.0209)
            

    if RunNumber==355544:
        if xscan:
            initialparameter(Fit,amplitude*0.66, amplitude*0.37, mean*1.08, mean*0.99 ,sigma, sigma*2)
    

    if RunNumber==355754:
        if Early and xscan:
            initialparameter(Fit,amplitude*0.64, amplitude*0.37, mean, mean ,sigma, sigma*2)

    

    if RunNumber==356177:
        if not xscan:

            initialparameter(Fit,amplitude*0.43, amplitude*0.58, mean, mean ,sigma*1.2, sigma*0.8)


    if RunNumber==356250:
        if xscan:
            initialparameter(Fit,amplitude*0.43, amplitude*0.58, mean, mean ,sigma*1.2, sigma*0.8)
        
    if RunNumber==357283:
        if Early and xscan:
            initialparameter(Fit,amplitude*0.77, amplitude*0.23, mean*1.57, mean*0.099 ,sigma*0.98, sigma*1.4)

    if RunNumber==357409:
        if Early and xscan:
            initialparameter(Fit,amplitude*0.77, amplitude*0.23, mean*1.57, mean*0.099 ,sigma*0.98, sigma*1.4)
        if Early and not xscan:
            initialparameter(Fit,0.286, 2.039,-0.00142, 0.000374 ,0.0163, 0.0215)

        else:
            initialparameter(Fit,amplitude, amplitude, mean, mean ,sigma, sigma*2)

    if RunNumber==357500:
        if xscan:
            initialparameter(Fit,amplitude*0.67, amplitude*0.35, mean*1.12, mean*0.99 ,sigma*0.8, sigma*1.2)

    if RunNumber==357539:
        if not xscan:
            initialparameter(Fit,amplitude, amplitude, mean, mean ,sigma, sigma*2)


    if RunNumber==357620:
        if Early and xscan:
            initialparameter(Fit,amplitude*0.38, amplitude*0.6, mean*2.4, mean*0.02 ,sigma*0.8, sigma*1.2)
        if Early and not xscan:
            initialparameter(Fit,amplitude*0.38, amplitude*0.6, mean*2.4, mean*0.02 ,sigma*0.8, sigma*1.2)

        if not Early and xscan:
            initialparameter(Fit,amplitude*0.79, amplitude*0.2, mean, mean ,sigma, sigma*2)

        if not Early and not xscan:
           
            initialparameter(Fit,amplitude*0.17, amplitude*0.83, -mean*1.8, mean*0.90 ,sigma*0.8, sigma*2)
            

    if RunNumber==357713:
        if Early and xscan:
            initialparameter(Fit,amplitude*0.95, amplitude*0.52, mean*0.09, mean*3 ,sigma, sigma*0.5)
        if Early and not xscan:
            initialparameter(Fit,0.9*amplitude,0.1*amplitude,mean,mean ,sigma,2*sigma)

        if not Early and xscan:
            initialparameter(Fit,0.787,0.00567, -0.00234, 0.0205 ,0.0137, 0.0205)

        if not Early and not xscan:
            initialparameter(Fit,amplitude*0.95, amplitude*0.52, mean*1.2, mean*1.8 ,sigma*1.09, sigma*0.5)
        
    if RunNumber==357887:
        if Early and xscan:
            initialparameter(Fit,amplitude*0.52, amplitude*0.49, mean*0.98, mean ,sigma*0.9, sigma*1.2)

        if Early and not xscan:
            initialparameter(Fit,amplitude*0.8, amplitude*0.2, mean, mean ,sigma, sigma*2)


        if not Early and xscan:
            initialparameter(Fit,amplitude, amplitude, mean, mean ,sigma, sigma*2)
        
    if RunNumber==358031:
        if Early and xscan:
            initialparameter(Fit,amplitude*0.84, amplitude*0.18, mean, mean*10 ,sigma*0.94, sigma*1.5)

        if not Early and xscan:
            initialparameter(Fit,amplitude*0.84, amplitude*0.18, mean, mean*10 ,sigma*0.94, sigma*1.5)

        if not Early and not xscan:
            initialparameter(Fit,amplitude*0.84, amplitude*0.18, mean, mean*10 ,sigma*0.94, sigma*1.5)

    if RunNumber==358300:
        if Early and xscan:
            initialparameter(Fit,0.873, 0.82,-0.00098, -0.00170 ,0.0112, 0.0156)

    if RunNumber==358395:
        if Early and xscan:
            initialparameter(Fit,amplitude*0.49, amplitude*0.52, mean, mean ,sigma*1.2, sigma*0.8)
        if not Early and xscan:
            initialparameter(Fit,amplitude, amplitude*0.1, mean, mean ,sigma, sigma*2)


    if RunNumber==358541:
        if not Early and xscan:
            initialparameter(Fit,amplitude*0.3, amplitude*0.7, mean, mean ,sigma, sigma*2)
        if not Early and not xscan:
            initialparameter(Fit,amplitude*0.35, amplitude*0.67, -mean*0.07, mean ,sigma*0.87, sigma*1.05)

    if RunNumber==358615:

        if Early and not xscan:
            
            initialparameter(Fit,amplitude*0.3, amplitude*0.75, mean, mean ,sigma, sigma*2)

        if not Early and xscan:
            initialparameter(Fit,0.003, amplitude*0.75, 0.0075, mean ,0.003, sigma*2)

        if not Early and not xscan:
            initialparameter(Fit,amplitude*0.3, amplitude*0.75, mean, mean ,sigma, sigma*2)


    if RunNumber==359058:
        if xscan:
            initialparameter(Fit,amplitude*0.7, amplitude*0.36, mean, mean*1.1 ,sigma*1.3, sigma*0.8)

    if RunNumber==359170:
        if xscan:
            initialparameter(Fit,amplitude*0.8, amplitude*0.18, mean*0.98, mean*1.4 ,sigma*0.98, sigma*1.6)
    

    if RunNumber==359191:
        if Early and xscan:

            initialparameter(Fit,0.5, 1.01, mean, mean ,sigma, sigma*2)

        if Early and not xscan:
            initialparameter(Fit,0.5*amplitude, 0.08, mean, mean ,sigma, sigma*2)


        if not Early and xscan:
            initialparameter(Fit,amplitude*0.17, amplitude*0.83, -mean*1.8, mean*0.90 ,sigma*0.8, sigma*2)

        if not Early and not xscan:
            initialparameter(Fit,amplitude*0.37, amplitude*0.67, mean, mean ,sigma, sigma*2)

    if RunNumber==359355:
        if xscan:
            initialparameter(Fit,0.580,1.07,-0.00187, -0.00244 ,0.0106, 0.01497)
 

    if RunNumber==359472:

        if xscan:
            initialparameter(Fit,0.580,1.07,-0.00187, -0.00244 ,0.0106, 0.01497)
        if not xscan:
            initialparameter(Fit,0.15,0.621,-0.00128, 0.00135 ,0.01417, 0.0197)
 

    if RunNumber==363400:
        if xscan:
            initialparameter(Fit,amplitude*0.97, amplitude*0.48, mean, mean*1.6 ,sigma*0.8, sigma*2.1)
        if not xscan:
            initialparameter(Fit,0.62, 1.58,-0.0024, 0.00207 ,0.013, 0.0198)

            # 1_x scan is problemeatic
    if RunNumber==363710:
        if Early and xscan:
            initialparameter(Fit,amplitude*0.56, amplitude*0.46, mean*0.98, mean ,sigma*0.9, sigma*1.3)
        if Early and not xscan:
            initialparameter(Fit,amplitude*0.9, amplitude*0.1, mean, mean ,sigma, sigma*2)

        if not Early and not xscan:
            initialparameter(Fit,amplitude*0.5, amplitude*0.67, mean, mean ,sigma, sigma*2.0)


    if RunNumber==363830:
        if Early and xscan:
            initialparameter(Fit,1.5, 0.8, -0.00179,-0.15 ,0.0128, 0.176)
        if not Early and not xscan:
            initialparameter(Fit,0.552, 0.0738,-0.00000122, 0.004905 ,0.016, 0.029)

    if RunNumber==363910:
        if xscan:
            initialparameter(Fit,0.552, 0.0738,-0.00000122, 0.004905 ,0.016, 0.029)
        if not xscan:
            initialparameter(Fit,0.552, 0.0738,-0.00000122, 0.004905 ,0.016, 0.029)


    if RunNumber==363947:
        if xscan:
            initialparameter(Fit,amplitude*0.4, amplitude*0.6, mean, mean ,sigma*1.3, sigma*0.9)


    if RunNumber==364030:
        if Early and xscan:
            initialparameter(Fit,0.64, 0.98,-0.000119, 0.000312 ,0.016, 0.0116)
        if Early and not xscan:
            initialparameter(Fit,amplitude*0.3, amplitude*0.7, mean, mean ,sigma, sigma*1.2)
        if not Early and xscan:
            initialparameter(Fit,amplitude, amplitude, mean, mean ,sigma, sigma*2.0)
        if not xscan and not xscan:
            initialparameter(Fit,amplitude*0.4, amplitude*0.7, -mean, mean ,sigma, sigma*2.0)


    if RunNumber==364076:
        if not xscan:
            initialparameter(Fit,amplitude*0.42, amplitude*0.6, mean, mean ,sigma, sigma*2.0)


    if RunNumber==364098:
        if Early and xscan:
            initialparameter(Fit,0.9,0.7, mean, mean*1.06 ,sigma*0.9, sigma*1.29)
        if Early and not xscan:
            initialparameter(Fit,amplitude*0.9, amplitude*0.1, mean, mean ,sigma, sigma*2.0)


        if not Early and xscan:

            initialparameter(Fit,0.058,0.767,-0.00618 , 0.00104 ,0.0148, 0.0141)

        if not Early and not xscan:
            initialparameter(Fit,amplitude, amplitude*0.6, mean, mean ,sigma, sigma*2.0)


    if RunNumber==364292:
        if Early and xscan:
            initialparameter(Fit,amplitude*0.35, amplitude*0.69, mean, mean ,sigma, sigma*1.1)   
    if RunNumber==348354:
        if Early and xscan:
            initialparameter(Fit,amplitude*0.19, amplitude*0.82, mean, mean ,sigma*1.9, sigma*0.7)   


    if RunNumber==348251:
        if Early and xscan:
            initialparameter(Fit,amplitude*0.37, amplitude*0.7, mean, -mean ,sigma*0.78, sigma*0.8)
        if Early and not xscan:
            initialparameter(Fit,amplitude*0.5, amplitude*0.5, -mean, mean ,sigma, sigma*2.0)
        else:
            initialparameter(Fit,amplitude, amplitude, mean, mean ,sigma, sigma*2.0)
    if RunNumber==352107:
        if xscan:
            initialparameter(Fit,amplitude*0.4, amplitude*0.6, mean*1.7, mean*0.5 ,sigma*0.8, sigma*1.2)


    if RunNumber==352274:
        if xscan:
            initialparameter(Fit,amplitude*0.41, amplitude*0.56, mean*0.91, mean*1.1 ,sigma*0.9, sigma*1.2)   
        if not xscan:
            initialparameter(Fit,amplitude*0.4, amplitude*0.6, mean*1.7, mean*0.5 ,sigma*0.8, sigma*1.2)   
    if RunNumber==354124:
        if xscan:
            initialparameter(Fit,amplitude*0.5, 0.82, mean, mean ,sigma*0.8, sigma*2)   
        if not xscan:
            initialparameter(Fit,2.19, 0.198, 0.000944, -0.0056 ,0.0213, 0.0168)   
    if RunNumber ==355389:
        if xscan:
            initialparameter(Fit,amplitude*0.42,amplitude*0.59, mean*4.5, mean*2,sigma*0.8, sigma*1.2)
    

    if RunNumber==355529:
        if not xscan:
            initialparameter(Fit,amplitude*0.1,amplitude*0.9, mean, mean,sigma, sigma*2)


    if RunNumber==357821:
        if not xscan:
            initialparameter(Fit,amplitude,amplitude, mean, mean,sigma, sigma*2)

        if not xscan and not Early:
            initialparameter(Fit,amplitude*0.18,amplitude*0.83, mean, mean,sigma, sigma*2)
    if RunNumber==359735:
        if Early and not xscan:
            initialparameter(Fit,1.95,0.735,-0.00138, 0.00810,0.01969, 0.0205)

        if not Early and xscan:
            initialparameter(Fit,amplitude*0.18,0.27, mean, mean,sigma, sigma*2)

        if not Early and not xscan:
            initialparameter(Fit,amplitude*0.18,amplitude*0.83, mean, mean,sigma, sigma*2)

    if RunNumber==355995:
        if Early and xscan:
            initialparameter(Fit,amplitude*0.5,amplitude*0.48, mean, mean, sigma, sigma*2)

        if Early and not xscan:
            initialparameter(Fit,amplitude*0.5,amplitude*0.48, mean, mean, sigma, sigma*2)

        if not Early and not xscan:
            initialparameter(Fit,amplitude*0.4,amplitude*0.6, mean, mean,sigma, sigma*2)


    if RunNumber==349011:
        if Early and xscan:
            initialparameter(Fit,amplitude*0.6,amplitude*0.4, mean, mean,sigma, sigma*1.3)
        if not Early:
            initialparameter(Fit,amplitude*0.9,amplitude*0.1, mean, mean,sigma, sigma*2)


    if RunNumber==349051:
        if Early and xscan:
            initialparameter(Fit,amplitude,amplitude, mean, mean,sigma, sigma*2)
        if Early and not xscan:
            initialparameter(Fit,amplitude*0.2,amplitude*0.8, mean, mean,sigma, sigma*2)
        if not Early and xscan:

            initialparameter(Fit,0.112,1.60, 0.051,-0.000196,0.075, 0.0150)
        if not Early and not xscan:
            initialparameter(Fit,amplitude,amplitude, mean, mean,sigma, sigma*2)

    if RunNumber==349268:
        if not Early and not xscan:
            initialparameter(Fit,amplitude*0.3,amplitude*0.7, mean, mean,sigma, sigma*1.0011)


    if RunNumber==350160:
        if not xscan:
            initialparameter(Fit,amplitude*0.2,amplitude*0.9, -mean, mean,sigma, sigma*1.2)

    if RunNumber==358516:
        if not xscan:
            initialparameter(Fit,amplitude,amplitude, mean, mean,sigma, sigma*2)

    if RunNumber==358656:
        if not xscan:
            initialparameter(Fit,amplitude,amplitude, mean, mean,sigma, sigma*2)
            

    if RunNumber==359124:
        if xscan:
            initialparameter(Fit,amplitude*0.3,amplitude*0.7, mean, mean,sigma, sigma*2)
        if not xscan:
            initialparameter(Fit,amplitude,amplitude, mean, mean,sigma, sigma*2)


    if RunNumber==359310:
        if not xscan:
            initialparameter(Fit,0.4,0.29, -0.000523, 0.001521,0.0156, 0.020)
    if RunNumber==359717:# problematic scan
        if not xscan:
            initialparameter(Fit,amplitude*0.4,amplitude*0.7, mean, mean,sigma, sigma*2)
    if RunNumber==359823:
        if Early and not xscan:
            initialparameter(Fit,amplitude*0.4,amplitude*0.7, mean, mean,sigma, sigma*2)

        if not Early and xscan:
            
            initialparameter(Fit,amplitude*0.3,amplitude*0.7, mean, mean,2*sigma, sigma)

        if not Early and not xscan:
            initialparameter(Fit,amplitude*0.3,amplitude*0.7, mean, mean,2*sigma, sigma)

    if RunNumber==359918:
        if xscan:
            initialparameter(Fit,0.9,0.7, -0.0021,-0.00189,0.0159, 0.0110)
        if not xscan:
            initialparameter(Fit,2.5,0.306,0.000826,0.00624,0.0255,0.0206 )


    if RunNumber==360244:
        if xscan:
            
            initialparameter(Fit,amplitude*0.9,amplitude*0.1, mean, mean,sigma, sigma*2)

        if not xscan:
            initialparameter(Fit,amplitude*0.1,amplitude, mean, mean,sigma*2, sigma)
        

    if RunNumber==360414:
        if xscan:
            initialparameter(Fit,amplitude*0.6,amplitude*0.4, mean, mean,sigma*2, sigma)
        if not xscan:
            initialparameter(Fit,amplitude*0.6,amplitude*0.4, mean, mean,sigma*2, sigma)

    if RunNumber==361738:

        if Early and not xscan:
            initialparameter(Fit,2.3,0.01, mean, mean,sigma, sigma*2)


        if not Early and xscan:
            initialparameter(Fit,amplitude*0.6,amplitude*0.6, mean, mean,sigma, sigma*1.2)

        if not Early and not xscan:
            initialparameter(Fit,amplitude*0.6,amplitude*0.6, mean, mean,sigma, sigma*1.2)


    if RunNumber==362445:

        if not Early and xscan:
            initialparameter(Fit,amplitude*0.5,amplitude*0.5, mean, mean,sigma, sigma*1.2)

        if Early and not xscan:
            initialparameter(Fit,amplitude*0.6,amplitude*0.5, mean, mean,sigma, sigma*2)


        if not Early and not xscan:
            initialparameter(Fit,amplitude*0.6,amplitude*0.0, mean, mean,sigma, sigma*1.2)


    if RunNumber==362661:
        if Early and xscan:
            initialparameter(Fit,amplitude*0.5,amplitude*0.5, mean, mean,sigma, sigma*1.2)
        if Early and not xscan:
            initialparameter(Fit,amplitude*0.5,amplitude*0.5, mean, mean,sigma, sigma*1.2)

        if not Early and xscan:
            initialparameter(Fit,0.934,0.036,-0.001311,0.00604,0.0146, 0.0276)

        if not Early and not xscan:
            initialparameter(Fit,0.934,0.036,-0.001311,0.00604,0.0146, 0.0276)


    if RunNumber==363033:

        initialparameter(Fit,amplitude,amplitude, mean, mean,sigma, sigma*2)
            
    if RunNumber==363198:
        if Early and xscan:
            initialparameter(Fit,0.85,0.74, mean, mean,sigma, sigma*1.2)
        if Early and not xscan:
            initialparameter(Fit,0.85,0.74, mean, mean,sigma, sigma*1.2)

        if not Early and xscan:
            initialparameter(Fit,0.47,0.0136, mean, mean,sigma, sigma*2)

        if not xscan and not Early:
            initialparameter(Fit,0.213,0.381, mean, mean,sigma, sigma*2)


    if RunNumber==363979:#problematic
        if xscan:
            initialparameter(Fit,0.26,0.37, 0.00076,0.04,0.0210, 0.176)

        if not xscan:
            initialparameter(Fit,0.26,0.37, 0.00076,0.04,0.0210, 0.176)



    if RunNumber==364214:
#        if not Early:
        if xscan:
            initialparameter(Fit,0.26,0.37, 0.00076,0.04,0.0210, 0.176)

        if not xscan:
            initialparameter(Fit,amplitude*0.6,amplitude*0.4, mean, mean,sigma, sigma*1.8)


    if RunNumber==357772:
        if xscan:
            initialparameter(Fit,0.46,0.869,-0.00361,-0.00331,0.0157, 0.0125)
        if not xscan:
            initialparameter(Fit,0.46,0.869,-0.00361,-0.00331,0.0157, 0.0125)
           
    if RunNumber==357679:
        if xscan:
            initialparameter(Fit,0.26,0.37, 0.00076,0.04,0.0210, 0.176)
        if not xscan:
            initialparameter(Fit,0.68,3.5, mean, mean,sigma, sigma*2)


    if RunNumber==351628:
        if xscan:
            initialparameter(Fit,0.68,3.5, mean, mean,sigma, sigma*2)
        if not xscan:
            initialparameter(Fit,2.0,0.0152, mean, mean,sigma, sigma*2)
        
    if RunNumber==350361:
        if Early and xscan:
            initialparameter(Fit,amplitude*0.6,amplitude*0.4, mean, mean,sigma, sigma*1.2)
        if Early and not xscan:
            initialparameter(Fit,0.48,2.0, mean, mean,sigma, sigma*1.2)

            
    if RunNumber==354359:
        if xscan:
            initialparameter(Fit,0.865,0.869,-0.00361,-0.00331,0.0157, 0.0125)
        if not xscan:
            initialparameter(Fit,amplitude,amplitude, mean, mean,sigma, sigma*1.2)

    if RunNumber==349592:
        if xscan:
            initialparameter(Fit,1.09,0.765, mean, mean,sigma, sigma)
        if not xscan:
            initialparameter(Fit,0.26,0.37, 0.00076,0.04,0.0210, 0.176)

    if RunNumber==349114:

        if not Early and xscan:
            initialparameter(Fit,0.566,0.59,-0.00133,0.000124,0.0158, 0.152)

    if RunNumber==349309:
        if xscan:
            initialparameter(Fit,0.6,2.5, mean, mean,sigma, sigma*2)

        if not xscan:
            initialparameter(Fit,0.6,2.5, mean, mean,sigma, sigma*2)


    if RunNumber==349033:
        if not Early and xscan:
            initialparameter(Fit,1.16,0.1317, mean, mean,sigma, sigma*2)
        if not Early and not xscan:
            initialparameter(Fit,amplitude,amplitude, mean, mean,sigma, sigma*2)

    if RunNumber==349014:
        if not xscan:
            initialparameter(Fit,amplitude,amplitude, mean, mean,sigma, sigma*2)



















