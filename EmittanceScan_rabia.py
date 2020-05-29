

#Rabia shaheen
#16_nov_2019

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
from initialPar import *




# importing pickle files

#fit_pickle = sorted(glob.glob("fits_4WP/integrated/fit_integrated_trk_run_350880.pickle"))
fit_pickle = sorted(glob.glob("fits_4WP/integrated/*trk*"))
#fit_pickle = sorted(glob.glob("test/*trk*"))


# 1d double gauusian fit function
def gaussD(x, p):                               

    return p[0] * (1.0 / ( math.sqrt(2 * math.pi) * p[2] ) ) * np.exp( - (1.0 / 2.0) * ( (x[0] - p[1])**2 / ( p[2] **2) ) ) + p[3] * (1.0 / (  math.sqrt(2 * math.pi) * p[4] ) ) * np.exp( - (1.0 / 2.0) * ( (x[0] - p[5])**2 / ( p[4] **2) ) )                         


# expected _luminosity fucntion      
freq = (1./8.89244e-05) #[s-1]
inel = 8. * 1e-24

def expected_luminosity( sigma_x, sigma_y, current ):
    return freq * current * 1e22 / ( sigma_x * sigma_y  * 2.0 * math.pi) / (freq / inel )



    # 3.main function (analyze function)
def analyse_scan(scan,separation, luminosity, error, path, run_number, bunches, current, Xscan,early,date):

    # separation --> separation of points(but not corrected one)
    # luminosity 
    # error
    # canvas_name 
    # run_number
    # bunches
    # current

    separation_points = len(separation) # number of points (8,9)depends on scan

    # TGraph(A Graph is a graphics object made of two arrays X and Y with npoints each.)
    graph = ROOT.TGraph(separation_points,separation,luminosity)

    #TGraphErrors(A TGraphErrors is a TGraph with error bars.)
    graph = ROOT.TGraphErrors(separation_points, separation, luminosity, np.zeros(separation_points), error ) 
    #TCanvas
    canvas_name = ROOT.TCanvas( 'scan_name', 'scan_name',800, 600 )

    # Canvas Title
    graph.GetXaxis().SetTitle("separation[mm]")
    graph.GetYaxis().SetTitle("mu")
    graph.GetXaxis().SetTitleSize(0.064)
    graph.GetYaxis().SetTitleSize(0.064)

    #Single gauss Fit 
    single_gauss_fit = ROOT.TF1( "single_gauss_fit", 'gaus', -0.05, 0.05 )# single gauss (3 parameters)
    #Double gauss Fit
    double_gauss_fit = ROOT.TF1("double_gauss_fit", gaussD, -0.05, 0.05, 6)#2d gauss (5 parameters)

    # fitting the graph with single gauss fit
    graph.Fit('single_gauss_fit','Q')
    # title
    graph.GetXaxis().SetTitle("separation")
    graph.GetYaxis().SetTitle("mu")


# for my fucntions

    amp1 = single_gauss_fit.GetParameter(0)
    #print amp1,"amp1"

    m1 = single_gauss_fit.GetParameter(1)
#    print m1,"m1"
    s1 = single_gauss_fit.GetParameter(2)
#    print s1,"s1"

    # getting value for the double gauss fit before we have only constant term

    amp1 = amp1 * s1 * math.sqrt(2*math.pi)


    ######## this is general initial values for fit paramters

    double_gauss_fit.SetParameter(0,0.5*amp1)#amplitude          
    double_gauss_fit.SetParameter(3,0.5*amp1)#amplitude2

    double_gauss_fit.SetParameter(4,2*s1)#sigma2
    double_gauss_fit.SetParameter(2,s1)#sigma1 

    double_gauss_fit.SetParameter(1,m1)#mean1 
    double_gauss_fit.SetParameter(5,m1)#mean
#    print bunches,"bunches"
    
    # Do a fit for a particular run number, X or Y, early or late
    # RunNumber = 123456
    # Xscan = True
    # Early = False
    # [A1, A2, s1, s2, m1, m2] = GetInitialParameters(RunNumber, Xscan, Early)
    # fit.SetParameter(0,A1) #amplitude                                                                                                            
    # fit.SetParameter(3,A2) #amp                                                                                                                  
    # fit.SetParameter(4,s2) #sigma2                                                                                                             
    # fit.SetParameter(1,m1) #mean1                                                                                                               
    # fit.SetParameter(2,s1) #sigma1                                                                                                             
    # fit.SetParameter(5,m2) #mean 

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






    # Fit double gauss                                                        
 
    legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
    legend.SetHeader("Double Gauss Fit")
    legend.AddEntry(graph,"data","P")
    legend.AddEntry(double_gauss_fit,"fit","L")
 

    
    status = graph.Fit('double_gauss_fit','SNQ')#fitting of double gaussain



    # this step is for those scans where fit fails

    if int(status) ==4:
        ROOT.TVirtualFitter.Fitter(graph).SetMaxIterations(10000)# increase iterations

        status = graph.Fit('double_gauss_fit','SNQ')#fitting of double gaussain

        # again fitting if fit fails
        if int(status)==4:
            status = graph.Fit('double_gauss_fit','SNQ')#fitting of double gaussain



    status_of_fit = int(status)
#    print status_of_fit,"status_of_fit"
    # chi2 from fit before correction
    chi2_double_gauss_fit = double_gauss_fit.GetChisquare()

    # NDF

    NDF_double_gauss_fit =  double_gauss_fit.GetNDF()

    # probabiltiy 
    Prob_double_gauss_fit = double_gauss_fit.GetProb()
    # chi/NDF

    chi_NDF_double_gauss_fit= chi2_double_gauss_fit/NDF_double_gauss_fit#ratio chi/NDF 
    # round off chi2

    round_chi_double_gauss_fit = round(chi2_double_gauss_fit,2)#round off the value    
    # round off chi2
    round_chi_NDF_double_gauss_fit = round(chi_NDF_double_gauss_fit,2)# 

#    print chi2_double_gauss_fit
 

    #cosmetics

    graph.SetLineColor( 38 )
    graph.SetMarkerColor( 4 )
    graph.SetMarkerStyle( 20 )
    graph.SetMarkerSize( 1.1 )
  
    # colors                            
    canvas_name.SetFillColor(18)
    canvas_name.SetGrid()
    canvas_name.SetGridx()
    canvas_name.SetGridy()

    # legend
    legend.SetTextSize(0.03)
    legend.Draw()

    graph.Draw("AP")

    # ****** without corrected separation *****    

    # chi2 value
    # luminsoity from fit (each point)
  
    chi2_uncorrected = np.array([])
    lumi_fit = gaussD([separation], double_gauss_fit.GetParameters())


    # looping on all separation points
    # to get lumionsity of each separation points
    # lastly calculate the chi2
    subtraction_lumi = np.array([])# array for subration_lumi

    for sep in range(separation_points):  # looping 
        subtraction_luminosity = (luminosity[sep]- lumi_fit[sep])
        subtraction_lumi= np.append(subtraction_lumi,subtraction_luminosity)
        chi2 = ((luminosity[sep]- lumi_fit[sep])/error[sep])**2# chi2 formula
        chi2_uncorrected = np.append(chi2_uncorrected,chi2)
        
    
    #canvas
    error_tgraph = ROOT.TGraphErrors(len(separation), separation,subtraction_lumi,np.zeros(separation_points), error)
    
    c = ROOT.TCanvas('sub_lumi', 'sub_lumi')
    c.Divide(1,2) # 1 row, 2 columns
    c.cd(1)
    graph.GetXaxis().SetTitle("separation[mm]")
    graph.GetYaxis().SetTitle("mu")
    graph.GetXaxis().SetTitleSize(0.064)
    graph.GetYaxis().SetTitleSize(0.064)


    graph.Draw("AP")

    legend = ROOT.TLegend(0.75, 0.6, 0.9, 0.86)
    legend.AddEntry(graph,"data","P")# points                                      
    legend.AddEntry(double_gauss_fit," Fit","L")# 
    legend.AddEntry(0,'#chi^{2}='+str(round_chi_double_gauss_fit), '')# chi2
    legend.AddEntry(0,'#chi^{2}/ndf='+str(round_chi_NDF_double_gauss_fit), '')# chi
    legend.SetTextSize(0.055)
   

    # latex
    label= ROOT.TLatex()
    label.SetTextSize(0.058)
    label.DrawLatexNDC(0.75,0.90,"#bf{Double Gauss Fit}")# horizontal, vertical
    label.SetTextColor(4)

    legend.Draw()


    c.cd(2)
    error_tgraph.Draw("AP")


    # line
    error_tgraph.SetLineColor( 38 )
    error_tgraph.SetMarkerColor( 4 )
    error_tgraph.SetMarkerStyle( 20 )
    error_tgraph.SetMarkerSize( 1.3 )




    #cosmetics
    error_tgraph.GetXaxis().SetTitle("Separation[mm]");
    error_tgraph.GetYaxis().SetTitle("lumi[data]-lumi[fit]");
    error_tgraph.GetXaxis().SetTitleSize(0.064)
    error_tgraph.GetYaxis().SetTitleSize(0.064)

    # line                                                                                                                                     
    line = TLine(-0.037,0,0.037,0.0)
    line.SetLineColor(2)
    line.SetLineWidth(1)
    line.Draw()





    # path
    fileName = path + str(run_number) +"_"+str(scan)+"_"
    if Xscan== True:
        fileName += "x_scan_comp.png"
    else:
        fileName += "y_scan_comp.png"
    c.Print(fileName)




       
     # path
    graph.Draw("AP")
    fileName = path + str(run_number) +"_"+str(scan)+"_"
    if Xscan == True:
        fileName = fileName +"x_scan.png"
    else:
        fileName = fileName+ "y_scan.png"
        
    canvas_name.Print(fileName)
    





    # covarinace matrix

    cov= status.GetCovarianceMatrix ()                                                                               
  

    # correlaation matrix
    corr= status.GetCorrelationMatrix ()                                                                             

    ### elements of Covariance matrix                                                                                
    e01 = cov(0,1)#1                                                                                                 
    e02 = cov(0,2)#2                                                                                                 
    e03 = cov(0,3)#3                                                                                                 
    e04 = cov(0,4)#4                                                                                                 
    e05 = cov(0,5)#5                                                                                                 
    e12 = cov(1,2)#6                                                                                                 
    e13 = cov(1,3)#7                                                                                                 
    e14 = cov(1,4)#8                                                                                                 
    e15 = cov(1,5)#9                                                                                                 
    e23 = cov(2,3)#10                                                                                                
    e24 = cov(2,4)#11                                                                                                
    e25 = cov(2,5)#12                                                                                                
    e34 = cov(3,4)#13                                                                                                
    e35 = cov(3,5)#14                                                                                                
    e45 = cov(4,5)#15                     
                                                                           
    # gives me status of covarinace matrix
    staus_cov_matrix = status.CovMatrixStatus()
#    print staus_cov_matrix,"staus_cov_matrix"
    # paramters from fit function
    # A1

    A1= double_gauss_fit.GetParameter(0)                                                                             

    error_a1=double_gauss_fit.GetParError(0)                                                                         
    relative_a1 = (error_a1/A1)* 100                                                                                 
    # A2


    A2= double_gauss_fit.GetParameter(3)                                                                   
    error_a2=double_gauss_fit.GetParError(3)                                  
    relative_a2 = (error_a2/A2)* 100                                          

#     #M1                                                                       
    M1= double_gauss_fit.GetParameter(1)                                      
    error_m1= double_gauss_fit.GetParError(1)                                 
    absolute_m1= error_m1                                                     



#     # M2                                                                      
    M2=double_gauss_fit.GetParameter(5)                                       
    error_m2= double_gauss_fit.GetParError(5)                                 
    absolute_m2= error_m2  



    # S1
    S1= double_gauss_fit.GetParameter(2)                                      
    error_s1= double_gauss_fit.GetParError(2)                                 
    relative_s1= (error_s1/S1)*100                                            

    # S2                                                                      

    S2= double_gauss_fit.GetParameter(4)                                      
    error_s2= double_gauss_fit.GetParError(4)                                 
    relative_s2= (error_s2/S2)*100                                            
#     print relative_s2,"s2 relative "  




    

    #  cap sigma function

    def cap_sigma():

        # maximum of the fit fucntion
        peak = double_gauss_fit.GetMaximum()
        
        # position of the peak
        Xmax=  double_gauss_fit.GetMaximumX()

        # limits for the integration
        
        intlimit = 0.1        
        
        # integral of the fit function( it gives me area under the curve)
        
        integral = double_gauss_fit.Integral(-intlimit,intlimit)
        
        # formula for finding the cap sigma
        sigma = (1 / math.sqrt (2 * math.pi)) * integral /peak 
        
        # derivatives of the paramters 

        # a1                                                                                                         
        der_a1= (1.0 / (math.sqrt(2 * math.pi) * S1) ) * np.exp( - (1.0 / 2.0) * ( ( Xmax - M1)**2 / (S1 **2) ) )# checked  


                 

         #a2                                                                                                          

        der_a2= (1.0 / ( math.sqrt(2 * math.pi) * S2) ) * np.exp( - (1.0 / 2.0) * ( ( Xmax - M2)**2 / (S2 **2) ) )# \checked
        
        
        # s1
        der_s1= - A1 * (1.0/ ( math.sqrt(2 * math.pi) * S1 **2) ) * np.exp( - (1.0 / 2.0) * ( ( Xmax - M1)**2 / (S1 **2) ) ) + A1 * (1.0/ ( math.sqrt(2 * math.pi) * S1 **4 ) )* (Xmax - M1)**2 * np.exp( - (1.0 / 2.0) * ( ( Xmax - M1)**2 / (S1 **2) ) )# checke                                                                                              


         #s2 
        der_s2= - A2 * (1.0/ ( math.sqrt(2 * math.pi) * S2 **2 ) ) * np.exp( - (1.0 / 2.0) * ( ( Xmax - M2)**2 / (S2 **2) ) ) + A2* (1.0/ ( math.sqrt(2 * math.pi) * S2 **4 ) )  * (Xmax - M2)**2 * np.exp( - (1.0 / 2.0) * ( ( Xmax - M2)**2 / (S2 **2) ) )# checked 

         #m1
        
        der_m1=  A1 * (Xmax - M1) * ( 1.0/ ( math.sqrt(2 * math.pi) * S1 **3) ) * np.exp( - (1.0 / 2.0) * ( ( Xmax - M1) **2 / (S1 **2) ) ) # checked



         #m2                                                                                                          
        der_m2=  A2 * (Xmax - M2) * ( 1.0/ ( math.sqrt(2 * math.pi) * S2 **3) ) * np.exp( - (1.0 / 2.0) * ( ( Xmax - M2) **2 / (S2 **2) ) ) #checked                                                                                      



        ### off diagnol terms
        
         #0                                                                                                           
        off_diagnol_0= 2.0 * der_a1 * der_m1 * e01 + 2.0 * der_a1 * der_s1 * e02 + 2.0 * der_a1 * der_a2 * e03 + 2.0 * der_a1 * der_s2 * e04 + 2.0 * der_a1 * der_m2 * e05



        off_diagnol_1 = 2.0 * der_m1 * der_s1 * e12 + 2.0 * der_m1 * der_a2 * e13 + 2.0 * der_m1 * der_s2 * e14 + 2.0* der_m1 * der_m2 * e15                                                                              



         #2  
        off_diagnol_2 = 2.0 * der_s1 * der_a2 * e23 + 2.0 * der_s1 * der_s2 * e24 + 2.0 * der_s1 * der_m2 * e25 


         #3                                                                                                           
        off_diagnol_3 = 2.0 * der_a2 * der_s2 * e34 + 2.0 * der_a2 * der_m2 *e35                              



         #4                                                                                                           
        off_diagnol_4 = 2.0 * der_s2 * der_m2 * e45                                                              




         #unceratnity in peak 
        peak_error = math.sqrt( (error_a1 * der_a1)**2 + (error_a2 * der_a2)**2 + (error_s1 * der_s1)**2 + (error_s2 * der_s2)**2 + (error_m1 * der_m1)**2 + (error_m2* der_m2)**2 + off_diagnol_0 + off_diagnol_1 + off_diagnol_2 + off_diagnol_3 + off_diagnol_4 )  



#        print peak_error,"peak error"
#        print peak,"peak"




             
        return sigma, peak, peak_error

    # calling the cap_sigma fucntion


    sigma_fit = cap_sigma()
    sigma = sigma_fit[0]
    peak =  sigma_fit[1]
    error_peak= sigma_fit[2]
    





    # 0  --> sigma
    # 1  --> chi/NDF
    # 2  --> peak
    # 3  --> chi2
    # 4  --> status of fit
    # 5  --> status of covarinace matrix
    # 6  --> peak error
    # 7  --> relative_a1
    # 8  --> relative_a2
    # 9  --> relative_s1
    # 10 --> relative_s2
    # 11 --> absolute_m1
    # 12 --> absolute_m2 
    # 13 --> A1
    # 14 --> A2
    # 15 --> S1
    # 16 --> S2
    # 17 --> M1
    # 18 --> M2


    return sigma, chi_NDF_double_gauss_fit, peak, chi2_double_gauss_fit, status_of_fit, staus_cov_matrix, error_peak,relative_a1,relative_a2, relative_s1, relative_s2, absolute_m1, absolute_m2, A1, A2, S1, S2, M1, M2


#### empty array portion


chi2_NDF_x=np.array([])
chi2_NDF_y=np.array([])


# peak (relative error)
rel_error_on_peak_x = np.array([])
rel_error_on_peak_y = np.array([])

# a1_x(relative error)
rel_error_on_a1_x  = np.array([])
rel_error_on_a2_x = np.array([])

# a2_y(relative error)
rel_error_on_a1_y  = np.array([])
rel_error_on_a2_y = np.array([])

#  relatuve error on sigma or width of the fit fucntion
# x scan

rel_error_on_s1_x  = np.array([])
rel_error_on_s2_x = np.array([])

# (relative error)
# y scan

rel_error_on_s1_y  = np.array([])
rel_error_on_s2_y = np.array([])

### absolute errors on m1
abs_error_on_m1_x  = np.array([])
abs_error_on_m2_x = np.array([])

# (absolute error) m2
# y scan

abs_error_on_m1_y  = np.array([])
abs_error_on_m2_y = np.array([])

### fit parameters values


# x scan
fit_A1_x =np.array([])
fit_A2_x =np.array([])
fit_S1_x =np.array([]) 
fit_S2_x =np.array([])
fit_M1_x =np.array([])
fit_M2_x =np.array([])

# y scan

fit_A1_y =np.array([])
fit_A2_y =np.array([])
fit_S1_y =np.array([]) 
fit_S2_y =np.array([])
fit_M1_y =np.array([])
fit_M2_y =np.array([])


step_1=np.array([])
step_2=np.array([])
step_3=np.array([])

step_4=np.array([])
step_5=np.array([])

FOM=np.array([])








# chi_NDF_early_x

for filename in fit_pickle:
    with open(filename,'rb') as fit_pickle:
        scan = pickle.load(fit_pickle)
        scans = scan['TightModLumi']

        for i_scan in range(len(scans)):
            output = scans[i_scan]
            #print i_scan,"scans"
            sep_x = output['sep_x']
            exp_lumi_x = output['exp_lumi_x']
            err_x = output['err_x']
            curr_0 = output['curr'][0]
            date = output['date']
            
            sep_y = output['sep_y']
            exp_lumi_y = output['exp_lumi_y']
            err_y = output['err_y']
            curr_1 = output['curr'][1]

            run = output['run']
            bunches= output['bunches']
            early = output['early']


            result_x = analyse_scan(i_scan, sep_x, exp_lumi_x ,err_x ,"Emittance_scan/", run, bunches, curr_0, True, early,date )
            result_y = analyse_scan(i_scan, sep_y, exp_lumi_y ,err_y ,"Emittance_scan/", run, bunches, curr_1, False, early,date)
            
        
            # intensities
            n1= np.mean(output['curr'][0])/output['bunches']
            n2= np.mean(output['curr'][1])/output['bunches']
            
            # product of intensities
            product = n1*n2
            



            # selection of the scans
            
            #xscans 

            if output['bunches'] >2000: # bunches cut
                step_1=np.append(step_1,output['run'])

                print len(step_1),"step1 bunches"

                if result_x[4]==0 and result_y[4]==0 :# converge fit

                    if result_x[4]>0 :
                        print result_x[4],"X fit is not coverging"
                        print output['run']
                        

                    if result_y[4]>0:
                        print result_y[4],"Y fit is not converging"
                        print output['run']

                    step_2= np.append(step_2, result_x[4])
                    step_3=np.append(step_3, result_y[4])

                    print len(step_2),"step_2 converging fit"
                    print len(step_3),"step_3 converging fit"



                    if result_x[5] == 3 and result_y[5]==3: ## covariance matrix cut


                        if result_x[5]<3 :
                            print result_x[5],"X covarinace matrix is not converging"
                            print output['run']

                        if result_y[5]<3:
                            print result_y[5],"Y covarinace matrix is not converging"
                            print output['run']


                        step_4=np.append(step_4,result_x[5])
                        step_5=np.append(step_5,result_y[5])


                        print len(step_4),"step_4 covarinace matrix"
                        print len(step_5),"step_5 covarinace matrix"

                        if result_x[1]<7 and result_y[1]<7:# chi/ndf cut


                            print len(chi2_NDF_x),"chi_NDF_x"
                            print len(chi2_NDF_y),"chi_NDF_y"
                            print output['run'],"run"
                            print result_x[1],"Xscan"
                            print result_y[1],"y scan"
                            



                            #chi2/ndf X scan
                            chi2_NDF_x= np.append(chi2_NDF_x,result_x[1])
                            histogram(100,0,4,chi2_NDF_x,"Chi_NDF_x","plots/Chi2_NDF_x.png","chi_NDF_x")
                            #chi2/ndf Y scan
                            chi2_NDF_y=np.append(chi2_NDF_y,result_y[1])
                            histogram(100,0,4,chi2_NDF_y,"Chi_NDF_y","plots/Chi2_NDF_y.png","chi_NDF_y")

                            # peak x relative error
                            relative_error_peak_x = (result_x[6]/result_x[2])*100
                            # appending peak x array
                            rel_error_on_peak_x = np.append(rel_error_on_peak_x,relative_error_peak_x)
                             # histogram
                            histogram(100,0,0.5,rel_error_on_peak_x,"rel_error_on_peak_x","plots/rel_error_on_peak_x.png","rel_error_on_peak_x[%]")
                            
                            # now expected luminosity
                            

                            expected_lumi =expected_luminosity(result_x[0],result_y[0],product)

                            # figure of merit

                            Figure_of_merit= result_x[2]/expected_lumi

                            # appending
                            FOM=np.append(FOM,Figure_of_merit)
                            # histogram
                            histogram(100,0.9,1.05,FOM,"FOM","plots/FOM.png","Figure of merit(peak/lumi)")

                        # here i will plot relative errors on all paramters


                        # relative error on a1   
                            rel_error_on_a1_x = np.append(rel_error_on_a1_x, result_x[7])
                        # histogram
                            histogram(100,0,2,rel_error_on_a1_x,"rel_error_on_a1_x","plots/rel_error_on_a1_x.png","rel_error_on_a1_x[%]")
                        
                        
                        # relative error on a2
                            rel_error_on_a2_x = np.append(rel_error_on_a2_x, result_x[8])
                        # histogram
                            histogram(100,0,3,rel_error_on_a2_x,"rel_error_on_a2_x","plots/rel_error_on_a2_x.png","rel_error_on_a2_x[%]")



                        # relative error on s1    
                            rel_error_on_s1_x = np.append(rel_error_on_s1_x, result_x[9])
                        # histogram
                            histogram(100,0,4,rel_error_on_s1_x,"rel_error_on_s1_x","plots/rel_error_on_s1_x.png","rel_error_on_s1_x[%]")
                        
                        
                        # relative error on s2
                            rel_error_on_s2_x = np.append(rel_error_on_s2_x, result_x[10])
                        # histogram
                            histogram(100,0,4,rel_error_on_s2_x,"rel_error_on_s2_x","plots/rel_error_on_s2_x.png","rel_error_on_s2_x[%]")
                             
                        # absolute error on the m1

                            abs_error_on_m1_x = np.append(abs_error_on_m1_x, result_x[11])
                        # histogram
                            histogram(100,-0.1,0.1,abs_error_on_m1_x,"abs_error_on_m1_x","plots/abs_error_on_m1_x.png","abs_error_on_m1_x")
                        
                        
                        # absolute error on m2
                            abs_error_on_m2_x = np.append(abs_error_on_m2_x, result_x[12])
                        # histogram
                            histogram(100,-0.1,0.1,abs_error_on_m2_x,"abs_error_on_m2_x","plots/abs_error_on_m2_x.png","abs_error_on_m2_x")


                        # here i am plotting the values for the paramters


                            # A1
                            fit_A1_x =np.append(fit_A1_x,result_x[13])
                            histogram(100,0.0,5,fit_A1_x,"fit_A1_x","plots/fit_A1_x.png","fit_A1_x[value]")



                            #A2
                            fit_A2_x =np.append(fit_A2_x,result_x[14])
                            histogram(100,0.0,5,fit_A2_x,"fit_A2_x","plots/fit_A2_x.png","fit_A2_x[value]")
                            
                            comp_histogram(100,0.0,10,fit_A1_x,fit_A2_x,"plots/fit_A1A2_x.png")
                            
                            


                            # S1
                            fit_S1_x =np.append(fit_S1_x,result_x[15])
                            histogram(100,0.01,0.09,fit_S1_x,"fit_S1_x","plots/fit_S1_x.png","fit_S1_x[value]")
                            
                            # S1
                            fit_S2_x =np.append(fit_S2_x,result_x[16])

                            histogram(100,0.01,0.09,fit_S2_x,"fit_S2_x","plots/fit_S2_x.png","fit_S2_x[value]")

                            comp_histogram(100,0.0,10,fit_S1_x,fit_S2_x,"plots/fit_S1S2_x.png")


                            # M1
                            fit_M1_x =np.append(fit_M1_x,result_x[17])

                            histogram(100,-0.1,0.1,fit_M1_x,"fit_M1_x","plots/fit_M1_x.png","fit_M1_x[value]")

                            #M2
                            fit_M2_x =np.append(fit_M2_x,result_x[18])

                            histogram(100,-0.1,0.1,fit_M2_x,"fit_M2_x","plots/fit_M2_x.png","fit_M2_x[value]")





                
                    


















    

            #             # relative error on a1

            #             # relative error on the sigma or width








            #             # fit parameters value
            #             # x scan


                        






            # if result_y[4]==0:# converge fit

            #     if output['bunches'] >2000: # bunches cut

            #         if result_y[5] == 3: ## covariance matriy cut



            #             if result_y[1]<4:

            #             # relative erron on peak of y scan
            #                 relative_error_peak_y = (result_y[6]/result_y[2])*100
                        
            #             # appending array
            #                 rel_error_on_peak_y = np.append(rel_error_on_peak_y,relative_error_peak_y)
                        
            #             # histogram for the peak y
            #                 histogram(100,0,1,rel_error_on_peak_y,"rel_error_on_peak_y","plots/rel_error_on_peak_y.png","rel_error_on_peak_y[%]")

            #             # relative error on a1
            #                 rel_error_on_a1_y = np.append(rel_error_on_a1_y, result_y[7])
                       
            #             # histogram
            #                 histogram(100,0,50,rel_error_on_a1_y,"rel_error_on_a1_y","plots/rel_error_on_a1_y.png","rel_error_on_a1_y[%]")
            #             # relative error on a2
                        
            #                 rel_error_on_a2_y = np.append(rel_error_on_a2_y, result_y[8])
                            
            #             # histogram
            #                 histogram(100,0,50,rel_error_on_a2_y,"rel_error_on_a2_y","plots/rel_error_on_a2_y.png","rel_error_on_a2_y[%]")

            #             # relative error on the sigma or width

            #                 rel_error_on_s1_y = np.append(rel_error_on_s1_y, result_y[9])

            #             # histogram
            #                 histogram(100,0,20,rel_error_on_s1_y,"rel_error_on_s1_y","plots/rel_error_on_s1_y.png","rel_error_on_s1_y[%]")
            #             # relative error on a2
            #                 rel_error_on_s2_y = np.append(rel_error_on_s2_y, result_y[10])
 
            #            # histogram
            #                 histogram(100,0,20,rel_error_on_s2_y,"rel_error_on_s2_y","plots/rel_error_on_s2_y.png","rel_error_on_s2_y[%]")

            #             # absolute error on the m1

            #                 abs_error_on_m1_y = np.append(abs_error_on_m1_y, result_y[11])
            #             # histogram
            #                 histogram(100,-0.0001,0.1,abs_error_on_m1_y,"abs_error_on_m1_y","plots/abs_error_on_m1_y.png","abs_error_on_m1_y[%]")
                        
                        
            #             # absolute error on m2
            #                 abs_error_on_m2_y = np.append(abs_error_on_m2_y, result_y[12])
            #             # histogram
            #                 histogram(100,-0.0001,0.1,abs_error_on_m2_y,"abs_error_on_m2_y","plots/abs_error_on_m2_y.png","abs_error_on_m2_y[%]")




            #             # fit paramters y scan

            #                 fit_A1_y =np.append(fit_A1_y,result_y[13])

            #                 fit_A2_y =np.append(fit_A2_y,result_y[14])
                            
            #                 fit_S1_y =np.append(fit_S1_y,result_y[15])
                            
            #                 fit_S2_y =np.append(fit_S2_y,result_y[16])

            #                 fit_M1_y =np.append(fit_M1_y,result_y[17])

            #                 fit_M2_y =np.append(fit_M2_y,result_y[18])


            #                 histogram(100,0,2,fit_A1_y,"fit_A1_y","plots/fit_A1_y.png","fit_A1_y[value]")
                            
            #                 histogram(100,0,2,fit_A2_y,"fit_A2_y","plots/fit_A2_y.png","fit_A2_y[value]")
                        
            #                 histogram(100,0.01,0.03,fit_S1_y,"fit_S1_y","plots/fit_S1_y.png","fit_S1_y[value]")
                        
            #                 histogram(100,0.01,0.04,fit_S2_y,"fit_S2_y","plots/fit_S2_y.png","fit_S2_y[value]")
                        
            #                 histogram(100,-0.1,0.1,fit_M1_y,"fit_M1_y","plots/fit_M1_y.png","fit_M1_y[value]")

            #                 histogram(100,-0.1,0.1,fit_M2_y,"fit_M2_y","plots/fit_M2_y.png","fit_M2_y[value]")
                        
