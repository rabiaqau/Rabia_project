
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
from functions import*
from ROOT import TCanvas, TPad, TFile, TPaveLabel, TPaveText


#
fit_pickle_trk = sorted(glob.glob("fits_4WP/integrated/fit_integrated_trk_run_349111.pickle"))
#fit_pickle_lcd = sorted(glob.glob("fits_4WP/integrated/*lcd*"))
#
#fit_pickle_trk = sorted(glob.glob("fits_4WP/integrated/*trk*"))




# def gaussD(x, p):                               

#     return p[0] *( 1.0 / (math.sqrt(2 * math.pi) * p[2]) *( p[3] * np.exp( - (1.0 / 2.0) *  ( (x[0] - p[1])**2 / p[2] **2 )  ) ) + ( 1.0  / (math.sqrt(2 * math.pi) * p[4])  )*((1- p[3]) * np.exp( - (1.0 / 2.0) * ( (x[0] - p[5])**2 / p[4] **2  ) ) ) )


def gaussD(x, p):

    return p[0] * (1.0 / ( math.sqrt(2 * math.pi) * p[2] ) ) * np.exp( - (1.0 /2.0) * ( (x[0] - p[1])**2 / ( p[2] **2) ) ) + p[3] * (1.0 / (  math.sqrt(2 * math.pi) * p[4] ) ) * np.exp( - (1.0 / 2.0) * ( (x[0] - p[5])**2 / ( p[4] **2) ))


# expected _luminosity         
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
    graph.Fit('single_gauss_fit','SNQ')



# for my fucntions

    amp1 = single_gauss_fit.GetParameter(0)

    error_amp1= double_gauss_fit.GetParError(0)
    #print (error_amp1/amp1)*100,"consant error"
    m1 = single_gauss_fit.GetParameter(1)

    error_m1= double_gauss_fit.GetParError(1)
   # print (error_m1/m1),"single mean error"
#    print bunches,"bunches"

#    print m1,"m1"
    s1 = single_gauss_fit.GetParameter(2)
    error_s1= double_gauss_fit.GetParError(2)
#    print (error_s1/s1)*100 ,"sigma error"


 #   print s1,"s1"

    # getting value for the double gauss fit before we have only constant term

    amplitude = amp1 * s1 * math.sqrt(2*math.pi)
  #  print amplitude,"amplitude"
#    print bunches,"bunches"

    m1_initial= m1
    s1_initial = s1
    s2_initial= 2*s1


    double_gauss_fit.SetParameter(4,s1*2)
    double_gauss_fit.SetParameter(1,m1)
    double_gauss_fit.SetParameter(2,s1)
    double_gauss_fit.SetParameter(5,m1)
    double_gauss_fit.SetParameter(0,amplitude*0.5)
    double_gauss_fit.SetParameter(3,amplitude*0.5)


    if run_number==359286:
        if Xscan:
            double_gauss_fit.SetParameter(4,s1*1.3)
            double_gauss_fit.SetParameter(2,s1*0.9)
            double_gauss_fit.SetParameter(0,amplitude*0.6)
            double_gauss_fit.SetParameter(3,amplitude*0.4)

        if not Xscan:
            double_gauss_fit.SetParameter(0,amplitude*0.5)
            double_gauss_fit.SetParameter(3,amplitude*0.5)







    if run_number==350013:
        if early and not Xscan:
            double_gauss_fit.SetParameter(4,s1*1.3)
            double_gauss_fit.SetParameter(0,amplitude*0.3)
            double_gauss_fit.SetParameter(3,amplitude*0.9)

        if not early and Xscan:
            double_gauss_fit.SetParameter(4,0.5)
            double_gauss_fit.SetParameter(2,0.0015)
            double_gauss_fit.SetParameter(0,0.711)
            double_gauss_fit.SetParameter(3,0.24)
            double_gauss_fit.SetParameter(1,-m1*0.1)

            double_gauss_fit.SetParameter(5,m1*0.5)

    # if run_number==350013:
    #     if not Xscan and not early:
    #         double_gauss_fit.SetParameter(4,s1*1.2)
    #         double_gauss_fit.SetParameter(0,amplitude*0.7)
    #         double_gauss_fit.SetParameter(3,amplitude*0.3)




    if run_number==350880:

        if early and Xscan:
            double_gauss_fit.SetParameter(0,amplitude*0.5)
            double_gauss_fit.SetParameter(3,amplitude*0.5)

        if not early and not Xscan:
            double_gauss_fit.SetParameter(0,amplitude)
            double_gauss_fit.SetParameter(3,amplitude)




        if early and not Xscan:
            double_gauss_fit.SetParameter(0,amplitude*0.9)
            double_gauss_fit.SetParameter(3,amplitude*0.1)
        if not early and not Xscan:
            double_gauss_fit.SetParameter(0,amplitude*0.9)
            double_gauss_fit.SetParameter(3,amplitude*0.1)

        # else:
        #     double_gauss_fit.SetParameter(0,amplitude)
        #     double_gauss_fit.SetParameter(3,amplitude)






    if run_number==350220:
        if Xscan:
            double_gauss_fit.SetParameter(2,s1*1.2)

            double_gauss_fit.SetParameter(0,amplitude*0.52)
            double_gauss_fit.SetParameter(3,amplitude*0.49)
        if not Xscan:

            double_gauss_fit.SetParameter(4,s1*1.1)
            double_gauss_fit.SetParameter(5,-m1*31)
            double_gauss_fit.SetParameter(3,amplitude*0.35)



    if run_number == 350144:
        if early and Xscan:
            
            double_gauss_fit.SetParameter(0,amplitude)
            double_gauss_fit.SetParameter(3,amplitude)
        if not early and Xscan:
            double_gauss_fit.SetParameter(5,m1*0.00098)



    if run_number==349327:
        if Xscan:

            double_gauss_fit.SetParameter(4,s1*1.2)
            double_gauss_fit.SetParameter(1,m1)
            double_gauss_fit.SetParameter(2,s1*0.8)
            double_gauss_fit.SetParameter(5,m1*1.08)
            double_gauss_fit.SetParameter(0,amplitude*0.43)
            double_gauss_fit.SetParameter(3,amplitude*0.54)





        if not Xscan:
            double_gauss_fit.SetParameter(0,amplitude)
            double_gauss_fit.SetParameter(3,amplitude)





    if run_number==348618:
        if not Xscan:
            double_gauss_fit.SetParameter(4,s1*1.1)
            double_gauss_fit.SetParameter(1,-m1*5.4)
            double_gauss_fit.SetParameter(2,s1*0.8)
            double_gauss_fit.SetParameter(5,m1*2.3)
            double_gauss_fit.SetParameter(0,amplitude*0.13)
            double_gauss_fit.SetParameter(3,amplitude*0.86)



    if run_number==349111:
        if Xscan:
            double_gauss_fit.SetParameter(4,s1*0.9)
            double_gauss_fit.SetParameter(2,s1*1.23)
            double_gauss_fit.SetParameter(1,m1*1.1)
            double_gauss_fit.SetParameter(5,m1*0.99)

            double_gauss_fit.SetParameter(0,amplitude*0.46)
            double_gauss_fit.SetParameter(3,amplitude*0.52)

    if run_number==357451:
        if not Xscan:

            double_gauss_fit.SetParameter(0,amplitude*0.6)
            double_gauss_fit.SetParameter(3,amplitude*0.37)




    if run_number== 349169:#### x early scan is still problem

        if early and Xscan:
            # double_gauss_fit.SetParameter(4,1.2*s1)
            # double_gauss_fit.SetParameter(5,1.1*s1)
            # double_gauss_fit.SetParameter(1,0.88*m1)
            # double_gauss_fit.SetParameter(2,0.8*s1)
            double_gauss_fit.SetParameter(0,amplitude*0.27)
            double_gauss_fit.SetParameter(3,amplitude*0.75)



        if early and not Xscan:
            double_gauss_fit.SetParameter(0,0.6*amplitude)
            double_gauss_fit.SetParameter(3,0.4*amplitude)
        if not early and not Xscan:
            double_gauss_fit.SetParameter(0,amplitude)
            double_gauss_fit.SetParameter(3,amplitude)






    if run_number==349335:
        if not Xscan:
#            double_gauss_fit.SetParameter(4,s1*0.72)
            double_gauss_fit.SetParameter(1,m1*0.99)
            double_gauss_fit.SetParameter(5,m1*3.7)
#            double_gauss_fit.SetParameter(0,amplitude*0.94)
            double_gauss_fit.SetParameter(3,amplitude*0.06)


    
    if run_number==349481:
        if early and Xscan:
            double_gauss_fit.SetParameter(5,m1*1.2)
            double_gauss_fit.SetParameter(1,m1*0.9)

        if  not early and  not Xscan:
            double_gauss_fit.SetParameter(4,s1)
            double_gauss_fit.SetParameter(0,amplitude*0.2)
            double_gauss_fit.SetParameter(3,amplitude*0.8)





    if run_number== 349533:

        if Xscan:
            double_gauss_fit.SetParameter(5,0.9*m1)
            double_gauss_fit.SetParameter(1,1.15*m1)
            double_gauss_fit.SetParameter(2,s1*1.32)
            double_gauss_fit.SetParameter(4,s1*0.93)
            double_gauss_fit.SetParameter(0,amplitude*0.39)
            double_gauss_fit.SetParameter(3,amplitude*0.63)


        if not Xscan:

            double_gauss_fit.SetParameter(0,amplitude)
            double_gauss_fit.SetParameter(3,amplitude)





    if run_number==349637:

        if not early and Xscan:
            double_gauss_fit.SetParameter(1,-m1)
            double_gauss_fit.SetParameter(0,amplitude*0.5)
            double_gauss_fit.SetParameter(3,amplitude*0.5)

        if not early and not Xscan:
#            double_gauss_fit.SetParameter(4,s1*0.8)
            double_gauss_fit.SetParameter(1,m1*2.9)
#            double_gauss_fit.SetParameter(2,s1*1.1)
            double_gauss_fit.SetParameter(5,-m1*4.5)
            double_gauss_fit.SetParameter(0,amplitude*0.76)
            double_gauss_fit.SetParameter(3,amplitude*0.24)



    if run_number== 349646:
        if Xscan:
            double_gauss_fit.SetParameter(5,m1*1.1)
            double_gauss_fit.SetParameter(0,amplitude*0.54)
            double_gauss_fit.SetParameter(3,amplitude*0.46)
            double_gauss_fit.SetParameter(2,s1*0.9)
            double_gauss_fit.SetParameter(4,s1*1.23)
        if not Xscan:
            double_gauss_fit.SetParameter(0,amplitude*0.7)
            double_gauss_fit.SetParameter(3,amplitude*0.3)
            double_gauss_fit.SetParameter(4,s1*1.2)





    if run_number==349693:
        if early and Xscan:
            double_gauss_fit.SetParameter(4,s1*1.3)
            double_gauss_fit.SetParameter(1,m1*0.99)
            double_gauss_fit.SetParameter(2,s1*0.9)
            double_gauss_fit.SetParameter(5,m1*1.08)
            double_gauss_fit.SetParameter(0,amplitude*0.53)
            double_gauss_fit.SetParameter(3,amplitude*0.48)


    if run_number==349841:
        if early and Xscan:
            double_gauss_fit.SetParameter(5,m1*2.5)
            # double_gauss_fit.SetParameter(1,m1*0.73)
            # double_gauss_fit.SetParameter(4,s1*1.4)
            # double_gauss_fit.SetParameter(2,s1*0.99)

            double_gauss_fit.SetParameter(0,amplitude*0.68)
            double_gauss_fit.SetParameter(3,amplitude*0.33)

        if not Xscan:
            double_gauss_fit.SetParameter(4,s1*2)
            double_gauss_fit.SetParameter(0,amplitude*0.9)
            double_gauss_fit.SetParameter(3,amplitude*0.1111)



    if run_number==349944:
        if early and Xscan:


            double_gauss_fit.SetParameter(4,s1*0.75)
            double_gauss_fit.SetParameter(1,m1)
            double_gauss_fit.SetParameter(2,s1*1.1)
            double_gauss_fit.SetParameter(5,m1*1.1)
            double_gauss_fit.SetParameter(0,amplitude*0.83)
            double_gauss_fit.SetParameter(3,amplitude*0.17)




        if early and not Xscan:
            double_gauss_fit.SetParameter(0,amplitude*0.9)
            double_gauss_fit.SetParameter(3,amplitude*0.1)
    
        if not early and Xscan:
            double_gauss_fit.SetParameter(0,amplitude*0.9)
            double_gauss_fit.SetParameter(3,amplitude*0.1)
    
        if not early and not Xscan:
            double_gauss_fit.SetParameter(0,amplitude*0.3)
            double_gauss_fit.SetParameter(3,amplitude*0.8)
            double_gauss_fit.SetParameter(4,s1*1.1)



    if run_number==349977:
        if Xscan:
            double_gauss_fit.SetParameter(4,s1*1.3)
            double_gauss_fit.SetParameter(1,m1*0.96)
            double_gauss_fit.SetParameter(2,s1*0.9)
            double_gauss_fit.SetParameter(5,m1*1.1)
            double_gauss_fit.SetParameter(0,amplitude*0.51)
            double_gauss_fit.SetParameter(3,amplitude*0.5)



    if run_number==350067:

        if Xscan:
            double_gauss_fit.SetParameter(4,s1*1.12)
            double_gauss_fit.SetParameter(1,m1*0.79)
            double_gauss_fit.SetParameter(2,s1*0.87)
            double_gauss_fit.SetParameter(5,m1*1.3)

            double_gauss_fit.SetParameter(0,amplitude*0.48)
            double_gauss_fit.SetParameter(3,amplitude*0.53)


        if not Xscan:
            double_gauss_fit.SetParameter(0,amplitude*0.6)
            double_gauss_fit.SetParameter(3,amplitude*0.4)


    if run_number==350121:


        if Xscan and early:
            # double_gauss_fit.SetParameter(4,s1*0.9)
            # double_gauss_fit.SetParameter(1,m1*0.99)
            # double_gauss_fit.SetParameter(2,s1*1.3)
            # double_gauss_fit.SetParameter(5,m1*1.1)
            double_gauss_fit.SetParameter(0,amplitude*0.52)
            double_gauss_fit.SetParameter(3,amplitude*0.49)



        if not early and Xscan:

            double_gauss_fit.SetParameter(4,1.28*s1)
            double_gauss_fit.SetParameter(0,amplitude*0.5)
            double_gauss_fit.SetParameter(3,amplitude*0.5)


        if not early and not Xscan:

            double_gauss_fit.SetParameter(0,amplitude)
            double_gauss_fit.SetParameter(3,amplitude)

    if run_number==350184:
        if early and Xscan:
            double_gauss_fit.SetParameter(1,m1*1.25)
            double_gauss_fit.SetParameter(5,m1*0.9)
            double_gauss_fit.SetParameter(2,s1*1.2)
            double_gauss_fit.SetParameter(4,s1*0.8)
            double_gauss_fit.SetParameter(0,amplitude*0.53)
            double_gauss_fit.SetParameter(3,amplitude*0.49)
        
        if early and not Xscan:
            double_gauss_fit.SetParameter(0,amplitude)
            double_gauss_fit.SetParameter(3,amplitude)
        if not early and Xscan:
            double_gauss_fit.SetParameter(2,s1*1.5)
            double_gauss_fit.SetParameter(4,s1)
            double_gauss_fit.SetParameter(0,amplitude*0.9)
            double_gauss_fit.SetParameter(3,amplitude*0.9)



    if run_number==350440:

        if Xscan and early:
            double_gauss_fit.SetParameter(4,s1*1.1)
            double_gauss_fit.SetParameter(0,amplitude*0.7)
            double_gauss_fit.SetParameter(3,amplitude)




        else:
            double_gauss_fit.SetParameter(0,amplitude)
            double_gauss_fit.SetParameter(3,amplitude)



    if run_number==350479:
        if Xscan and early:
            double_gauss_fit.SetParameter(4,s1*1.1)
            double_gauss_fit.SetParameter(1,m1*0.64)
            double_gauss_fit.SetParameter(2,s1*0.74)
            double_gauss_fit.SetParameter(5,m1*1.27)
            double_gauss_fit.SetParameter(0,amplitude*0.27)
            double_gauss_fit.SetParameter(3,amplitude*0.74)


        if Xscan and not early:
            double_gauss_fit.SetParameter(4,s1*1.41)
            double_gauss_fit.SetParameter(1,m1*0.99)
            double_gauss_fit.SetParameter(2,s1*0.99)
            double_gauss_fit.SetParameter(5,m1*100)
            double_gauss_fit.SetParameter(0,amplitude*0.98)
            double_gauss_fit.SetParameter(3,amplitude*0.7)



        else:    
            double_gauss_fit.SetParameter(0,amplitude)
            double_gauss_fit.SetParameter(3,amplitude)



    if run_number==350531:
        if early and Xscan:
            double_gauss_fit.SetParameter(4,s1*0.9)
            double_gauss_fit.SetParameter(1,m1*1.3)
            double_gauss_fit.SetParameter(2,s1*1.2)
            double_gauss_fit.SetParameter(5,m1*0.6)
            double_gauss_fit.SetParameter(0,amplitude*0.58)
            double_gauss_fit.SetParameter(3,amplitude*0.42)
        if Xscan and not early:

            double_gauss_fit.SetParameter(4,s1*4.7)
#            double_gauss_fit.SetParameter(1,m1*0.99)
#            double_gauss_fit.SetParameter(2,s1)
            double_gauss_fit.SetParameter(5,m1*8)
#            double_gauss_fit.SetParameter(0,amplitude*0.96)
            double_gauss_fit.SetParameter(3,amplitude*0.08)





    if run_number==350682:
        if early and Xscan:
            double_gauss_fit.SetParameter(0,amplitude*0.9)
            double_gauss_fit.SetParameter(4,s1*1.3)
            double_gauss_fit.SetParameter(3,amplitude*0.3)

        if early and not Xscan:
            
            double_gauss_fit.SetParameter(0,amplitude)

            double_gauss_fit.SetParameter(3,amplitude)
        if not early:

            double_gauss_fit.SetParameter(0,amplitude)
                        
            double_gauss_fit.SetParameter(3,amplitude)



    if run_number==350842:
        if early and Xscan:
            double_gauss_fit.SetParameter(4,s1*1.3)
            double_gauss_fit.SetParameter(1,m1*1.09)
            double_gauss_fit.SetParameter(2,s1*0.9)
            double_gauss_fit.SetParameter(5,m1)
            double_gauss_fit.SetParameter(0,amplitude*0.44)
            double_gauss_fit.SetParameter(3,amplitude*0.57)



    if run_number==351062:
        if Xscan:
            double_gauss_fit.SetParameter(4,s1*0.9)
            double_gauss_fit.SetParameter(1,m1*1.07)
            double_gauss_fit.SetParameter(2,s1*1.3)
            double_gauss_fit.SetParameter(5,m1*0.98)
            double_gauss_fit.SetParameter(0,amplitude*0.44)
            double_gauss_fit.SetParameter(3,amplitude*0.57)
        if not Xscan:
            double_gauss_fit.SetParameter(4,s1)
            double_gauss_fit.SetParameter(1,-m1*2.4)
            double_gauss_fit.SetParameter(2,s1*0.99)
            double_gauss_fit.SetParameter(5,m1*7)
            double_gauss_fit.SetParameter(0,amplitude*0.89)
            double_gauss_fit.SetParameter(3,amplitude*0.17)




    if run_number==351223:
        if early and Xscan:
            double_gauss_fit.SetParameter(4,s1*1.2)
            double_gauss_fit.SetParameter(1,m1*0.89)
            double_gauss_fit.SetParameter(2,s1*0.9)
            double_gauss_fit.SetParameter(5,m1*1.2)
            double_gauss_fit.SetParameter(0,amplitude*0.53)
            double_gauss_fit.SetParameter(3,amplitude*0.48)


        if not Xscan and not early:


            double_gauss_fit.SetParameter(0,amplitude)
            double_gauss_fit.SetParameter(3,amplitude)







    if run_number==351364:

        if early and Xscan:
            double_gauss_fit.SetParameter(4,s1*0.9)
            double_gauss_fit.SetParameter(1,m1*0.99)
            double_gauss_fit.SetParameter(2,s1*1.2)
            double_gauss_fit.SetParameter(5,m1*1.1)
            double_gauss_fit.SetParameter(0,amplitude*0.41)
            double_gauss_fit.SetParameter(3,amplitude*0.6)


        if not early and Xscan:
            double_gauss_fit.SetParameter(0,amplitude*0.9)
            double_gauss_fit.SetParameter(3,amplitude*0.1)
        if not early and not Xscan:
            double_gauss_fit.SetParameter(4,s1*1.23)
            double_gauss_fit.SetParameter(1,-m1*0.4)
            double_gauss_fit.SetParameter(2,s1*0.9)
            double_gauss_fit.SetParameter(5,m1*7)
            double_gauss_fit.SetParameter(0,amplitude*0.76)
            double_gauss_fit.SetParameter(3,amplitude*0.23)









    if run_number==351550:
        if Xscan:
            double_gauss_fit.SetParameter(4,s1*1.23)
            double_gauss_fit.SetParameter(1,m1*9.9)
            double_gauss_fit.SetParameter(2,s1*0.9)
            double_gauss_fit.SetParameter(5,m1*12)
            double_gauss_fit.SetParameter(0,amplitude*0.64)
            double_gauss_fit.SetParameter(3,amplitude*0.37)






    if run_number==351671:
        if Xscan and early:
            double_gauss_fit.SetParameter(4,s1*0.9)
            double_gauss_fit.SetParameter(1,m1*1.2)
            double_gauss_fit.SetParameter(2,s1*1.3)
            double_gauss_fit.SetParameter(5,m1*0.99)
            double_gauss_fit.SetParameter(0,amplitude*0.29)
            double_gauss_fit.SetParameter(3,amplitude*0.7)



        if Xscan and not early:
            # double_gauss_fit.SetParameter(4,s1*2)
            # double_gauss_fit.SetParameter(1,m1)
            # double_gauss_fit.SetParameter(2,s1)
            # double_gauss_fit.SetParameter(5,m1)
            double_gauss_fit.SetParameter(0,amplitude)
            double_gauss_fit.SetParameter(3,amplitude)




    if run_number==351969:
        if Xscan and early:
            double_gauss_fit.SetParameter(4,s1*0.9)
            double_gauss_fit.SetParameter(1,m1*1.1)
            double_gauss_fit.SetParameter(2,s1*1.2)
            double_gauss_fit.SetParameter(5,m1*0.99)
            double_gauss_fit.SetParameter(0,amplitude*0.4)
            double_gauss_fit.SetParameter(3,amplitude*0.58)

            
        if not  Xscan and not early:
            double_gauss_fit.SetParameter(4,s1)
            double_gauss_fit.SetParameter(1,-m1*1.8)
            double_gauss_fit.SetParameter(2,s1*0.8)
            double_gauss_fit.SetParameter(5,m1*1.71)
            double_gauss_fit.SetParameter(0,amplitude*0.18)
            double_gauss_fit.SetParameter(3,amplitude*0.83)



    if run_number==352107:
        if Xscan:
            double_gauss_fit.SetParameter(4,s1*1.2)
            double_gauss_fit.SetParameter(1,m1*1.7)
            double_gauss_fit.SetParameter(2,s1*0.8)
            double_gauss_fit.SetParameter(5,m1*0.5)
            double_gauss_fit.SetParameter(0,amplitude*0.4)
            double_gauss_fit.SetParameter(3,amplitude*0.6)


    if run_number==352274:
        if Xscan:
            double_gauss_fit.SetParameter(4,s1*1.2)
            double_gauss_fit.SetParameter(1,m1*0.91)
            double_gauss_fit.SetParameter(2,s1*0.9)
            double_gauss_fit.SetParameter(5,m1*1.1)
            double_gauss_fit.SetParameter(0,amplitude*0.41)
            double_gauss_fit.SetParameter(3,amplitude*0.56)
        if not Xscan:
            double_gauss_fit.SetParameter(4,s1*1.2)
            double_gauss_fit.SetParameter(1,m1*1.7)
            double_gauss_fit.SetParameter(2,s1*0.8)
            double_gauss_fit.SetParameter(5,m1*0.5)
            double_gauss_fit.SetParameter(0,amplitude*0.4)
            double_gauss_fit.SetParameter(3,amplitude*0.6)



    if run_number==352448:
        if not Xscan:
            double_gauss_fit.SetParameter(4,s1*0.8)
            double_gauss_fit.SetParameter(1,m1*3.2)
            double_gauss_fit.SetParameter(2,s1*1.1)
            double_gauss_fit.SetParameter(5,m1*2.4)
            double_gauss_fit.SetParameter(0,amplitude*0.8)
            double_gauss_fit.SetParameter(3,amplitude*0.2)





    if run_number==354124:
        if Xscan:
            double_gauss_fit.SetParameter(4,s1*2)
            double_gauss_fit.SetParameter(1,m1)
            double_gauss_fit.SetParameter(2,s1)
            double_gauss_fit.SetParameter(5,m1)
            double_gauss_fit.SetParameter(0,amplitude*0.5)
            double_gauss_fit.SetParameter(3,amplitude*0.5)

        if not Xscan:
            double_gauss_fit.SetParameter(4,s1*0.8)
            double_gauss_fit.SetParameter(1,m1*2.9)
            double_gauss_fit.SetParameter(2,s1*1.04)
            double_gauss_fit.SetParameter(5,-m1*15)
            double_gauss_fit.SetParameter(0,amplitude*0.9)
            double_gauss_fit.SetParameter(3,amplitude*0.08)



    if run_number==354396:
        if Xscan:

            double_gauss_fit.SetParameter(4,s1*1.23)
            double_gauss_fit.SetParameter(1,m1*0.3)
            double_gauss_fit.SetParameter(2,s1*0.9)
            double_gauss_fit.SetParameter(5,m1*3)
            double_gauss_fit.SetParameter(0,amplitude*0.73)
            double_gauss_fit.SetParameter(3,amplitude*0.26)


        if not Xscan:
            double_gauss_fit.SetParameter(4,s1*1.15)
            double_gauss_fit.SetParameter(1,-m1*1.5)
            double_gauss_fit.SetParameter(2,s1*0.99)
            double_gauss_fit.SetParameter(5,m1*10)
            double_gauss_fit.SetParameter(0,amplitude*0.8)
            double_gauss_fit.SetParameter(3,amplitude*0.22)




    if run_number==355261:
        if early:
            if Xscan:

                double_gauss_fit.SetParameter(4,s1*1.1)
                double_gauss_fit.SetParameter(1,m1*0.8)
                double_gauss_fit.SetParameter(2,s1*0.8)
                double_gauss_fit.SetParameter(5,m1*1.15)
                double_gauss_fit.SetParameter(0,amplitude*0.36)
                double_gauss_fit.SetParameter(3,amplitude*0.63)




            if not Xscan:    
                double_gauss_fit.SetParameter(0,amplitude)
                double_gauss_fit.SetParameter(3,amplitude)

        if not early:
            if not Xscan:
                double_gauss_fit.SetParameter(0,amplitude*0.3)
                double_gauss_fit.SetParameter(3,amplitude*0.7)





    if run_number==355389:
        if Xscan:
            double_gauss_fit.SetParameter(4,s1*1.2)
            double_gauss_fit.SetParameter(1,m1*4.5)
            double_gauss_fit.SetParameter(2,s1*0.8)
            double_gauss_fit.SetParameter(5,m1*2)
            double_gauss_fit.SetParameter(0,amplitude*0.42)
            double_gauss_fit.SetParameter(3,amplitude*0.59)
    



    if run_number== 355416:
        if Xscan:
            double_gauss_fit.SetParameter(4,s1*0.8)
            double_gauss_fit.SetParameter(1,m1*4.5)
            double_gauss_fit.SetParameter(2,s1*1.2)
            double_gauss_fit.SetParameter(5,m1*2)
            double_gauss_fit.SetParameter(0,amplitude*0.42)
            double_gauss_fit.SetParameter(3,amplitude*0.59)

        if not Xscan:
            double_gauss_fit.SetParameter(4,s1*1.05)
            double_gauss_fit.SetParameter(1,-m1*14)
            double_gauss_fit.SetParameter(2,s1*0.8)
            double_gauss_fit.SetParameter(5,m1*0.7)
            double_gauss_fit.SetParameter(0,amplitude*0.1)
            double_gauss_fit.SetParameter(3,amplitude*0.9)






    if run_number==355468:
        if early and Xscan:


            double_gauss_fit.SetParameter(4,s1*0.8)
            double_gauss_fit.SetParameter(1,m1*1.03)
            double_gauss_fit.SetParameter(2,s1*1.3)
            double_gauss_fit.SetParameter(5,m1*0.98)

            double_gauss_fit.SetParameter(0,amplitude*0.39)
            double_gauss_fit.SetParameter(3,amplitude*0.45)





    if run_number==355529:
        if not Xscan:
#            double_gauss_fit.SetParameter(4,s1)
#            double_gauss_fit.SetParameter(1,-m1*3.2)
#            double_gauss_fit.SetParameter(2,s1*0.8)
#            double_gauss_fit.SetParameter(5,m1*1.4)
            double_gauss_fit.SetParameter(0,amplitude*0.1)
            double_gauss_fit.SetParameter(3,amplitude*0.9)





    if run_number== 355529:
        if not Xscan:


            double_gauss_fit.SetParameter(0,amplitude*0.9)
            double_gauss_fit.SetParameter(3,amplitude*0.1)





    if run_number==355544:
        if Xscan:
#            double_gauss_fit.SetParameter(4,s1*0.8)
            double_gauss_fit.SetParameter(1,m1*1.08)
#            double_gauss_fit.SetParameter(2,s1*1.2)
            double_gauss_fit.SetParameter(5,m1*0.99)
            double_gauss_fit.SetParameter(0,amplitude*0.66)
            double_gauss_fit.SetParameter(3,amplitude*0.37)




    if run_number==355754:
        if Xscan and early:
            double_gauss_fit.SetParameter(0,amplitude*0.64)
            double_gauss_fit.SetParameter(3,amplitude*0.37)

        if not Xscan and not early:

            double_gauss_fit.SetParameter(4,s1*0.99)
#            double_gauss_fit.SetParameter(1,m1*13)
#            double_gauss_fit.SetParameter(2,s1*1.8)
            double_gauss_fit.SetParameter(5,m1*0.45)
            double_gauss_fit.SetParameter(0,amplitude*0.089)
            double_gauss_fit.SetParameter(3,amplitude*0.9)




    if run_number==356177:
        double_gauss_fit.SetParameter(0,amplitude)
        double_gauss_fit.SetParameter(3,amplitude)



    if run_number==356250:
        if Xscan:

            double_gauss_fit.SetParameter(4,s1*1.2)
            double_gauss_fit.SetParameter(2,s1*0.8)
            double_gauss_fit.SetParameter(0,amplitude*0.43)
            double_gauss_fit.SetParameter(3,amplitude*0.58)


    if run_number==357283:
        if early and Xscan:
            double_gauss_fit.SetParameter(4,s1*1.4)
            double_gauss_fit.SetParameter(1,m1*1.57)
            double_gauss_fit.SetParameter(2,s1*0.98)
            double_gauss_fit.SetParameter(5,m1*0.099)
            double_gauss_fit.SetParameter(0,amplitude*0.77)
            double_gauss_fit.SetParameter(3,amplitude*0.23)



    if run_number==357409:
        if early and Xscan:


            double_gauss_fit.SetParameter(4,s1*1.2)
            double_gauss_fit.SetParameter(1,m1*0.99)
            double_gauss_fit.SetParameter(2,s1*0.8)
            double_gauss_fit.SetParameter(5,m1*1.33)
            double_gauss_fit.SetParameter(0,amplitude*0.43)
            double_gauss_fit.SetParameter(3,amplitude*0.59)


        else:


            double_gauss_fit.SetParameter(0,amplitude)
            double_gauss_fit.SetParameter(3,amplitude)



    if run_number==357500:
        if Xscan:
            double_gauss_fit.SetParameter(3,amplitude*0.67)
            double_gauss_fit.SetParameter(0,amplitude*0.35)
            double_gauss_fit.SetParameter(4,s1*1.2)
            double_gauss_fit.SetParameter(1,m1*1.12)
            double_gauss_fit.SetParameter(2,s1*0.8)
            double_gauss_fit.SetParameter(5,m1*0.99)




    if run_number==357539:

        if Xscan:

            double_gauss_fit.SetParameter(4,s1*1.02)
            double_gauss_fit.SetParameter(1,m1*2.8)
            double_gauss_fit.SetParameter(2,s1*0.78)
            double_gauss_fit.SetParameter(5,-m1*0.099)
            double_gauss_fit.SetParameter(3,amplitude*0.96)
            double_gauss_fit.SetParameter(0,amplitude*0.32)



        if not Xscan:

            double_gauss_fit.SetParameter(0,amplitude)
            double_gauss_fit.SetParameter(3,amplitude)



    if run_number==357620:

        if early and Xscan:
            double_gauss_fit.SetParameter(4,s1*1.2)
            double_gauss_fit.SetParameter(1,m1*2.4)
            double_gauss_fit.SetParameter(2,s1*0.8)
            double_gauss_fit.SetParameter(5,m1*0.02)
            double_gauss_fit.SetParameter(0,amplitude*0.38)
            double_gauss_fit.SetParameter(3,amplitude*0.6)


        if not early and Xscan:
            double_gauss_fit.SetParameter(0,amplitude*0.79)
            double_gauss_fit.SetParameter(3,amplitude*0.2)
        if not early and not Xscan:
            double_gauss_fit.SetParameter(4,s1)
            double_gauss_fit.SetParameter(1,-m1*1.8)
            double_gauss_fit.SetParameter(2,s1*0.8)
            double_gauss_fit.SetParameter(5,m1*0.90)
            double_gauss_fit.SetParameter(0,amplitude*0.17)
            double_gauss_fit.SetParameter(3,amplitude*0.83)



    if run_number==357713:
        if early and Xscan:
            double_gauss_fit.SetParameter(4,s1*0.5)
            double_gauss_fit.SetParameter(1,m1*0.09)
            double_gauss_fit.SetParameter(2,s1)
            double_gauss_fit.SetParameter(5,m1*3)
            double_gauss_fit.SetParameter(0,amplitude*0.95)
            double_gauss_fit.SetParameter(3,amplitude*0.52)


        if early and not Xscan:

            double_gauss_fit.SetParameter(0,amplitude*0.9)
            double_gauss_fit.SetParameter(3,amplitude*0.1)
        if not early and not Xscan:


            double_gauss_fit.SetParameter(4,s1*0.5)
            double_gauss_fit.SetParameter(1,m1*1.2)
            double_gauss_fit.SetParameter(2,s1*1.09)
            double_gauss_fit.SetParameter(5,m1*1.8)
            double_gauss_fit.SetParameter(0,amplitude*0.95)
            double_gauss_fit.SetParameter(3,amplitude*0.52)




    if run_number==357887:
        if early and Xscan:

            double_gauss_fit.SetParameter(4,s1*1.2)
            double_gauss_fit.SetParameter(1,m1*0.98)
            double_gauss_fit.SetParameter(2,s1*0.9)
            double_gauss_fit.SetParameter(5,m1)
            double_gauss_fit.SetParameter(0,amplitude*0.52)
            double_gauss_fit.SetParameter(3,amplitude*0.49)




        if not early and not Xscan:


            double_gauss_fit.SetParameter(5,m1*10)
            double_gauss_fit.SetParameter(1,-m1*10)
            double_gauss_fit.SetParameter(4,1.1*s1)
            double_gauss_fit.SetParameter(2,0.8*s1)
            double_gauss_fit.SetParameter(0,amplitude*0.25)
            double_gauss_fit.SetParameter(3,amplitude*0.76)








    if run_number==358031:
        if Xscan and early:
            double_gauss_fit.SetParameter(4,s1*1.5)
            double_gauss_fit.SetParameter(2,s1*0.94)
            double_gauss_fit.SetParameter(5,m1*10)
            double_gauss_fit.SetParameter(0,amplitude*0.84)
            double_gauss_fit.SetParameter(3,amplitude*0.18)

        if not early and not Xscan:
#            double_gauss_fit.SetParameter(4,s1*0.83)
#            double_gauss_fit.SetParameter(2,s1*1.1)
            double_gauss_fit.SetParameter(1,m1*0.32)
#            double_gauss_fit.SetParameter(5,m1*2.1)
#            double_gauss_fit.SetParameter(0,amplitude*0.65)
            double_gauss_fit.SetParameter(3,amplitude*0.36)

        
        # else:
        #     double_gauss_fit.SetParameter(4,s1*0.9)
        #     double_gauss_fit.SetParameter(2,s1*1.2)
        #     double_gauss_fit.SetParameter(0,amplitude*0.45)
        #     double_gauss_fit.SetParameter(3,amplitude*0.56)



    if run_number==358300:
        if Xscan and early: 
            double_gauss_fit.SetParameter(4,s1*0.9)
            double_gauss_fit.SetParameter(1,m1*1.1)
            double_gauss_fit.SetParameter(2,s1*1.2)
            double_gauss_fit.SetParameter(5,m1*0.09)


            double_gauss_fit.SetParameter(0,amplitude*0.42)
            double_gauss_fit.SetParameter(3,amplitude*0.52)



    if run_number==358395:
        if early and Xscan:
            double_gauss_fit.SetParameter(2,s1*1.2)
            double_gauss_fit.SetParameter(4,s1*0.8)
            double_gauss_fit.SetParameter(0,amplitude*0.49)
            double_gauss_fit.SetParameter(3,amplitude*0.52)

        if not early and Xscan:
            double_gauss_fit.SetParameter(0,amplitude)
            double_gauss_fit.SetParameter(3,amplitude)


    if run_number==358541:
        if not early and not Xscan:
            double_gauss_fit.SetParameter(4,s1*1.05)
            double_gauss_fit.SetParameter(1,-m1*0.07)
            double_gauss_fit.SetParameter(2,s1*0.87)
            double_gauss_fit.SetParameter(5,m1)


            double_gauss_fit.SetParameter(0,amplitude*0.35)
            double_gauss_fit.SetParameter(3,amplitude*0.67)




        if not early and Xscan:

            double_gauss_fit.SetParameter(0,amplitude*0.3)
            double_gauss_fit.SetParameter(3,amplitude*0.7)


    if run_number==358615:
        if early and Xscan:
            double_gauss_fit.SetParameter(4,s1*1.3)
            double_gauss_fit.SetParameter(5,m1*1.06)
            double_gauss_fit.SetParameter(2,s1*0.9)

            double_gauss_fit.SetParameter(0,amplitude*0.54)
            double_gauss_fit.SetParameter(3,amplitude*0.46)




    if run_number==358615:
        if early and Xscan:

            double_gauss_fit.SetParameter(0,amplitude)
            double_gauss_fit.SetParameter(3,amplitude)




        if not early and not Xscan:
            double_gauss_fit.SetParameter(0,amplitude*0.3)
            double_gauss_fit.SetParameter(3,amplitude*0.75)




    if run_number==359058:
        if Xscan:

            double_gauss_fit.SetParameter(4,s1*0.8)
#            double_gauss_fit.SetParameter(1,m1*)
            double_gauss_fit.SetParameter(2,s1*1.3)
            double_gauss_fit.SetParameter(5,m1*1.1)

            double_gauss_fit.SetParameter(0,amplitude*0.7)
            double_gauss_fit.SetParameter(3,amplitude*0.36)



    if run_number==359170:
        if Xscan:

            double_gauss_fit.SetParameter(4,s1*1.6)
            double_gauss_fit.SetParameter(1,m1*0.98)
            double_gauss_fit.SetParameter(2,s1*0.98)
            double_gauss_fit.SetParameter(5,m1*1.4)

            double_gauss_fit.SetParameter(0,amplitude*0.8)
            double_gauss_fit.SetParameter(3,amplitude*0.18)


    if run_number==359171:

        if Xscan:
            double_gauss_fit.SetParameter(4,s1*1.8)
#            double_gauss_fit.SetParameter(1,m1*1.7)
#            double_gauss_fit.SetParameter(2,s1)
            double_gauss_fit.SetParameter(5,-m1*60)
            
#            double_gauss_fit.SetParameter(0,amplitude*0.98)
            double_gauss_fit.SetParameter(3,amplitude*0.031)

        if not Xscan:
            double_gauss_fit.SetParameter(0,amplitude)
            double_gauss_fit.SetParameter(3,amplitude)





    if run_number==359191:

        if not early and not Xscan:
            double_gauss_fit.SetParameter(0,amplitude*0.37)
            double_gauss_fit.SetParameter(3,amplitude*0.67)



    if run_number==359355:
        if Xscan:

#            double_gauss_fit.SetParameter(4,s1*1.1)
            double_gauss_fit.SetParameter(1,m1*0.9)
#            double_gauss_fit.SetParameter(2,s1*0.8)
            double_gauss_fit.SetParameter(5,m1*1.2)
            # double_gauss_fit.SetParameter(0,amplitude*0.35)
            # double_gauss_fit.SetParameter(3,amplitude*0.67)




    if run_number==359472:
        if not early:
            if not Xscan:

#                double_gauss_fit.SetParameter(4,s1*1.1)
                double_gauss_fit.SetParameter(1,-m1*1.7)
#                double_gauss_fit.SetParameter(2,s1*0.8)
                double_gauss_fit.SetParameter(5,m1*1.7)
        #        double_gauss_fit.SetParameter(0,amplitude*0.2)
        #        double_gauss_fit.SetParameter(3,amplitude*0.8)








    if run_number==363400:
        if Xscan:
            double_gauss_fit.SetParameter(5,m1*1.6)
            double_gauss_fit.SetParameter(4,s1*2.1)
            double_gauss_fit.SetParameter(2,s1*0.8)
            double_gauss_fit.SetParameter(0,amplitude*0.97)
            double_gauss_fit.SetParameter(3,amplitude*0.48)
        if not Xscan:
            double_gauss_fit.SetParameter(1,-m1*1.3)
            double_gauss_fit.SetParameter(1,-m1*1.3)
            double_gauss_fit.SetParameter(0,amplitude*0.38)
            double_gauss_fit.SetParameter(3,amplitude*0.96)
            double_gauss_fit.SetParameter(4,s1)
            double_gauss_fit.SetParameter(2,s1*0.7)


    if run_number==363664:
        if Xscan:
            double_gauss_fit.SetParameter(0,amplitude*0.9)
            double_gauss_fit.SetParameter(3,amplitude*0.1)




    if run_number==363710:
        if early and Xscan:
            double_gauss_fit.SetParameter(1,0.98*m1)
            double_gauss_fit.SetParameter(1,m1*1.08)

            double_gauss_fit.SetParameter(4,s1*1.3)
            double_gauss_fit.SetParameter(2,s1*0.9)
            double_gauss_fit.SetParameter(0,amplitude*0.56)
            double_gauss_fit.SetParameter(3,amplitude*0.46)
        if early and not Xscan:

            double_gauss_fit.SetParameter(0,amplitude*0.9)
            double_gauss_fit.SetParameter(3,amplitude*0.1)

        if not early and Xscan:

            double_gauss_fit.SetParameter(0,amplitude)
            double_gauss_fit.SetParameter(3,amplitude)

        if not early and not Xscan:

            double_gauss_fit.SetParameter(0,amplitude*0.5)
            double_gauss_fit.SetParameter(3,amplitude*0.67)



    if run_number==363738:
        if not Xscan:
            double_gauss_fit.SetParameter(1,-m1*0.06)
            double_gauss_fit.SetParameter(4,s1*0.9)
            double_gauss_fit.SetParameter(2,s1*1.2)
            double_gauss_fit.SetParameter(0,amplitude*0.58)
            double_gauss_fit.SetParameter(3,amplitude*0.43)

    if run_number==363830:
        if early and Xscan:
            double_gauss_fit.SetParameter(0,amplitude)
            double_gauss_fit.SetParameter(3,amplitude)
        if not early and not Xscan:
            double_gauss_fit.SetParameter(0,amplitude*0.9)
            double_gauss_fit.SetParameter(4,1.7*s1)
            double_gauss_fit.SetParameter(1,0.001*m1)
            double_gauss_fit.SetParameter(3,amplitude*0.12)


    if run_number==363910:
        double_gauss_fit.SetParameter(0,amplitude)
        double_gauss_fit.SetParameter(3,amplitude)


    if run_number==363947:
        if Xscan:
            double_gauss_fit.SetParameter(2,s1*1.2)
            double_gauss_fit.SetParameter(4,s1*0.9)
            double_gauss_fit.SetParameter(0,amplitude*0.4)
            double_gauss_fit.SetParameter(3,amplitude*0.6)



    if run_number==364030:

        if early and Xscan:
            double_gauss_fit.SetParameter(4,s1*0.9)
            double_gauss_fit.SetParameter(2,s1*1.3)
            double_gauss_fit.SetParameter(0,amplitude*0.4)
            double_gauss_fit.SetParameter(3,amplitude*0.6)



        if early and not Xscan:
            double_gauss_fit.SetParameter(4,s1*1.2)
            double_gauss_fit.SetParameter(0,amplitude*0.3)
            double_gauss_fit.SetParameter(3,amplitude*0.7)
        if not early and Xscan:
            double_gauss_fit.SetParameter(0,amplitude)
            double_gauss_fit.SetParameter(3,amplitude)
        if not early and not Xscan:
            double_gauss_fit.SetParameter(1,-m1)
            double_gauss_fit.SetParameter(0,amplitude*0.4)
            double_gauss_fit.SetParameter(3,amplitude*0.7)






    if run_number==364076:
        if not Xscan:
            double_gauss_fit.SetParameter(0,amplitude*0.42)
            double_gauss_fit.SetParameter(3,amplitude*0.6)




    if run_number==364098:
        if early and Xscan:
            double_gauss_fit.SetParameter(5,m1*1.06)
            double_gauss_fit.SetParameter(4,1.29*s1)
            double_gauss_fit.SetParameter(2,0.9*s1)

            double_gauss_fit.SetParameter(0,amplitude*0.57)
            double_gauss_fit.SetParameter(3,amplitude*0.44)



        if early and not Xscan:

            double_gauss_fit.SetParameter(0,amplitude*0.9)
            double_gauss_fit.SetParameter(3,amplitude*0.1)


        if not early:

            double_gauss_fit.SetParameter(0,amplitude)
            double_gauss_fit.SetParameter(3,amplitude*0.6)


    if run_number==364292:
        if early and Xscan:
            double_gauss_fit.SetParameter(0,amplitude*0.35)
            double_gauss_fit.SetParameter(4,1.1*s1)
            double_gauss_fit.SetParameter(3,amplitude*0.69)





    if run_number==348354:
        if early and Xscan:
            double_gauss_fit.SetParameter(2,s1*1.9)
            double_gauss_fit.SetParameter(4,s1*0.7)
            double_gauss_fit.SetParameter(3,amplitude*0.19)
            double_gauss_fit.SetParameter(0,amplitude*0.82)



    if run_number==348251:
        if early and Xscan:
            double_gauss_fit.SetParameter(5,-m1)
            double_gauss_fit.SetParameter(2,s1*0.78)
            double_gauss_fit.SetParameter(4,s1*0.8)
            double_gauss_fit.SetParameter(0,amplitude*0.37)
            double_gauss_fit.SetParameter(3,amplitude*0.7)
        if early and not Xscan:
            double_gauss_fit.SetParameter(1,-m1)
            double_gauss_fit.SetParameter(0,amplitude*0.5)
            double_gauss_fit.SetParameter(3,amplitude*0.5)


        else:
            double_gauss_fit.SetParameter(0,amplitude)
            double_gauss_fit.SetParameter(3,amplitude)



    if run_number==351969:

        if not early:
            if Xscan:
                double_gauss_fit.SetParameter(0,amplitude*0.06)

                double_gauss_fit.SetParameter(3,amplitude*0.9)
            if not Xscan:
                double_gauss_fit.SetParameter(4,s1*1.0001)
                double_gauss_fit.SetParameter(0,amplitude*0.8)
                double_gauss_fit.SetParameter(3,amplitude*0.2)

    if run_number==357821:
        if not Xscan:
            double_gauss_fit.SetParameter(0,amplitude)
            double_gauss_fit.SetParameter(3,amplitude)



        if not Xscan and not early:
#            double_gauss_fit.SetParameter(4,s1)
            double_gauss_fit.SetParameter(0,amplitude*0.18)
            double_gauss_fit.SetParameter(3,amplitude*0.83)


    if run_number==359735:
        if not early and not Xscan:
            double_gauss_fit.SetParameter(0,amplitude*0.3)
            double_gauss_fit.SetParameter(3,amplitude*0.75)


    if run_number== 355995:
        if Xscan and early:
#            double_gauss_fit.SetParameter(4,s1*2)
            double_gauss_fit.SetParameter(0,amplitude*0.5)
            double_gauss_fit.SetParameter(3,amplitude*0.48)

        if not Xscan and not early:

            double_gauss_fit.SetParameter(0,amplitude*0.4)
            double_gauss_fit.SetParameter(3,amplitude*0.6)



    if run_number==350676:

        if not early:


            double_gauss_fit.SetParameter(0,amplitude)
            double_gauss_fit.SetParameter(3,amplitude)

        else:
            
            double_gauss_fit.SetParameter(0,amplitude)
            double_gauss_fit.SetParameter(3,amplitude)


    if run_number==349011:
        if early:
            if Xscan:
                double_gauss_fit.SetParameter(4,s1*1.3)
                double_gauss_fit.SetParameter(0,amplitude*0.6)
                double_gauss_fit.SetParameter(3,amplitude*0.4)
        if not early:

            double_gauss_fit.SetParameter(0,amplitude*0.9)
            double_gauss_fit.SetParameter(3,amplitude*0.1)


    if run_number==349051:
        if early and Xscan:
            double_gauss_fit.SetParameter(0,amplitude)
            double_gauss_fit.SetParameter(3,amplitude)

        if early and not Xscan:
            double_gauss_fit.SetParameter(0,amplitude*0.2)
            double_gauss_fit.SetParameter(3,amplitude*0.8)


        if not early and not Xscan:
            double_gauss_fit.SetParameter(0,amplitude)
            double_gauss_fit.SetParameter(3,amplitude)

   


    if run_number==349268:
        if not early and not Xscan:
            double_gauss_fit.SetParameter(4,s1*1.0011)
            double_gauss_fit.SetParameter(0,amplitude*0.3)
            double_gauss_fit.SetParameter(3,amplitude*0.7)



            
    if run_number==349592:
        if not Xscan:
            double_gauss_fit.SetParameter(4,s1*0.1)
            double_gauss_fit.SetParameter(0,amplitude)
            double_gauss_fit.SetParameter(3,amplitude*0.1)



    if run_number==350160:
        if not early:
            double_gauss_fit.SetParameter(0,amplitude*0.2)
            double_gauss_fit.SetParameter(3,amplitude*0.9)
            double_gauss_fit.SetParameter(1,-m1)
            double_gauss_fit.SetParameter(4,s1*1.2)




    if run_number==358516:
        if not Xscan:
            double_gauss_fit.SetParameter(0,amplitude)
            double_gauss_fit.SetParameter(3,amplitude)



    if run_number==358656:
        if not Xscan:
            double_gauss_fit.SetParameter(0,amplitude)
            double_gauss_fit.SetParameter(3,amplitude)




    if run_number==359124:
        if Xscan:

            double_gauss_fit.SetParameter(0,amplitude*0.3)
            double_gauss_fit.SetParameter(3,amplitude*0.7)

        if not Xscan:
            double_gauss_fit.SetParameter(0,amplitude)
            double_gauss_fit.SetParameter(3,amplitude)




    if run_number==359310:
        double_gauss_fit.SetParameter(1,-m1)
        double_gauss_fit.SetParameter(4,s1*1.4)
        double_gauss_fit.SetParameter(0,amplitude*0.4)
        double_gauss_fit.SetParameter(3,amplitude*0.6)




    if run_number==359717:
        double_gauss_fit.SetParameter(0,amplitude*0.4)
        double_gauss_fit.SetParameter(3,amplitude*0.7)




    if run_number==359823:
        if not early and not Xscan:
            double_gauss_fit.SetParameter(4,s1)
            double_gauss_fit.SetParameter(2,2*s1)
            double_gauss_fit.SetParameter(0,amplitude*0.3)
            double_gauss_fit.SetParameter(3,amplitude*0.7)




    if run_number==359918:
        if Xscan:

            double_gauss_fit.SetParameter(4,s1*0.8)
            double_gauss_fit.SetParameter(0,amplitude*0.6)
            double_gauss_fit.SetParameter(3,amplitude*0.4)




        if not Xscan:
            double_gauss_fit.SetParameter(4,s1*1.3)
            double_gauss_fit.SetParameter(0,amplitude)
            double_gauss_fit.SetParameter(3,amplitude)




    if run_number==360244:
        if Xscan:

            double_gauss_fit.SetParameter(0,amplitude*0.9)
            double_gauss_fit.SetParameter(3,amplitude*0.1)
        if not Xscan:

            double_gauss_fit.SetParameter(0,amplitude)
            double_gauss_fit.SetParameter(3,amplitude*0.1)


    if run_number==360414:
        double_gauss_fit.SetParameter(0,amplitude)
        double_gauss_fit.SetParameter(3,amplitude)





    if run_number==361738:
        if not early and not Xscan:
            double_gauss_fit.SetParameter(4,s1*1.2)
            double_gauss_fit.SetParameter(0,amplitude*0.6)
            double_gauss_fit.SetParameter(3,amplitude*0.6)




    if run_number== 361862:

        if not early and not Xscan:
            double_gauss_fit.SetParameter(0,amplitude)
            double_gauss_fit.SetParameter(3,amplitude)



    if run_number==362445:
        if not early:
            if Xscan:
                double_gauss_fit.SetParameter(4,s1*1.2)
                double_gauss_fit.SetParameter(0,amplitude*0.5)
                double_gauss_fit.SetParameter(3,amplitude*0.5)
            if not Xscan:
                double_gauss_fit.SetParameter(0,amplitude)
                double_gauss_fit.SetParameter(3,amplitude*0.0)

    if run_number==362661:
        if early:
            if Xscan:
                double_gauss_fit.SetParameter(4,s1*1.2)
                double_gauss_fit.SetParameter(0,amplitude*0.5)
                double_gauss_fit.SetParameter(3,amplitude*0.5)
            if not Xscan:
                double_gauss_fit.SetParameter(0,amplitude*0.5)
                double_gauss_fit.SetParameter(3,amplitude*0.5)
        if not early:

            if Xscan:
                double_gauss_fit.SetParameter(0,amplitude*0.9)
                double_gauss_fit.SetParameter(3,amplitude*0.1)
                double_gauss_fit.SetParameter(5,-m1)
            if not Xscan:
                double_gauss_fit.SetParameter(0,amplitude*0.4)
                double_gauss_fit.SetParameter(3,amplitude*0.7)
                double_gauss_fit.SetParameter(4,s1*1.5)


    if run_number==363033:
        double_gauss_fit.SetParameter(0,amplitude)
        double_gauss_fit.SetParameter(3,amplitude)



    if run_number==363198:
        if early:
            if Xscan:
                double_gauss_fit.SetParameter(4,s1)
                double_gauss_fit.SetParameter(0,amplitude*0.5)
                double_gauss_fit.SetParameter(3,amplitude*0.6)
            if not Xscan:
                double_gauss_fit.SetParameter(0,amplitude)
                double_gauss_fit.SetParameter(3,amplitude)

        if not early:
            if Xscan:
                double_gauss_fit.SetParameter(4,s1)
                double_gauss_fit.SetParameter(0,amplitude*0.5)
                double_gauss_fit.SetParameter(3,amplitude*0.5)
            if not Xscan:
                double_gauss_fit.SetParameter(0,amplitude)
                double_gauss_fit.SetParameter(3,amplitude)


    if run_number==363979:

        if Xscan:

            double_gauss_fit.SetParameter(0,amplitude*0.3)
            double_gauss_fit.SetParameter(3,amplitude*2.0)
        if not Xscan:
            double_gauss_fit.SetParameter(0,amplitude)
            double_gauss_fit.SetParameter(3,amplitude)



    if run_number==364214:
        if not early:
            double_gauss_fit.SetParameter(4,s1*1.8)
            double_gauss_fit.SetParameter(0,amplitude*0.6)
            double_gauss_fit.SetParameter(3,amplitude*0.4)

    if run_number==357887:

        if early and Xscan:
            double_gauss_fit.SetParameter(0,amplitude*0.5)
            double_gauss_fit.SetParameter(3,amplitude*0.5)



        if early and not Xscan:
            double_gauss_fit.SetParameter(0,amplitude*0.8)
            double_gauss_fit.SetParameter(3,amplitude*0.2)



        if not early and Xscan:
            double_gauss_fit.SetParameter(0,amplitude)
            double_gauss_fit.SetParameter(3,amplitude)
        if not early and not Xscan:

            double_gauss_fit.SetParameter(0,amplitude)
            double_gauss_fit.SetParameter(3,amplitude)



    if run_number==357772:

        if Xscan:
            double_gauss_fit.SetParameter(4,s1)
            double_gauss_fit.SetParameter(0,amplitude*0.4)
            double_gauss_fit.SetParameter(3,amplitude*0.7)



        if not Xscan:

            double_gauss_fit.SetParameter(0,amplitude*0.6)
            double_gauss_fit.SetParameter(3,amplitude*0.3)



    if run_number== 357679:
        if Xscan:
            double_gauss_fit.SetParameter(0,amplitude*0.9)
            double_gauss_fit.SetParameter(3,amplitude*0.1)


        if not Xscan:
            double_gauss_fit.SetParameter(0,amplitude)
            double_gauss_fit.SetParameter(3,amplitude)


    if run_number==355331:
        if Xscan:
            double_gauss_fit.SetParameter(0,amplitude)
            double_gauss_fit.SetParameter(3,amplitude)



    if run_number==351628:
        double_gauss_fit.SetParameter(0,amplitude)
        double_gauss_fit.SetParameter(3,amplitude)



    if run_number==350361:
        if early and Xscan:
            double_gauss_fit.SetParameter(4,s1*1.2)
            double_gauss_fit.SetParameter(0,amplitude*0.6)
            double_gauss_fit.SetParameter(3,amplitude*0.4)


        if early and not Xscan:
            double_gauss_fit.SetParameter(4,s1)
            double_gauss_fit.SetParameter(0,amplitude)
            double_gauss_fit.SetParameter(3,amplitude*0.2)



    if run_number==354359:
        double_gauss_fit.SetParameter(0,amplitude)
        double_gauss_fit.SetParameter(3,amplitude)



    if run_number==349637:
        if early and Xscan:

            double_gauss_fit.SetParameter(0,amplitude)
            double_gauss_fit.SetParameter(3,amplitude)


        if not early and not Xscan:    
            double_gauss_fit.SetParameter(0,amplitude*2)
            double_gauss_fit.SetParameter(3,amplitude)
            double_gauss_fit.SetParameter(2,0.0194)
            double_gauss_fit.SetParameter(4,0.013)
            double_gauss_fit.SetParameter(5,-m1)
            
            
    if run_number== 349592:
        if Xscan:
            double_gauss_fit.SetParameter(0,1.09)
            double_gauss_fit.SetParameter(3,0.765)
            double_gauss_fit.SetParameter(2,s1)

        if not Xscan:
            double_gauss_fit.SetParameter(0,amplitude*0.9)
            double_gauss_fit.SetParameter(3,amplitude*0.0005)


    if run_number== 349309:
        if not Xscan:
            double_gauss_fit.SetParameter(0,0.6*amplitude)
            double_gauss_fit.SetParameter(3,0.4*amplitude)


    if run_number== 349268:
        if not early and not Xscan:
            double_gauss_fit.SetParameter(0,0.6*amplitude)
            double_gauss_fit.SetParameter(3,0.4*amplitude)



    if run_number==349114:
        if Xscan and early:
            double_gauss_fit.SetParameter(0,amplitude)
            double_gauss_fit.SetParameter(3,amplitude)


    if run_number==349033:
        if not early:
            double_gauss_fit.SetParameter(0,amplitude)
            double_gauss_fit.SetParameter(3,amplitude)

    if run_number==349014:
        if not Xscan:
            double_gauss_fit.SetParameter(0,amplitude)
            double_gauss_fit.SetParameter(3,amplitude)

    # legend
    legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
    legend.SetHeader("Double Gauss Fit")
    legend.AddEntry(graph,"data","P")
    legend.AddEntry(double_gauss_fit,"fit","L")

    # do the fit but save it in to variable "status"

    status= graph.Fit('double_gauss_fit','SN')#fitting of double gaussain



    # status== 4 means the fit failed
    if int(status)== 4:

        # it increases the no of iterations

        ROOT.TVirtualFitter.Fitter(graph).SetMaxIterations(10000)    

        # second time fitting
        status =graph.Fit('double_gauss_fit','SN')


        if int(status)==4:
            # if 2nd time fails again do fitting

            status = graph.Fit('double_gauss_fit','SN')
    graph.Draw("AP")

    # it gives me overall status of the fit
    overall_status = int(status)
    graph.Draw("AP")
    # chi2 from fit
    chi2_double_gauss_fit = double_gauss_fit.GetChisquare()
    graph.Draw("AP")
    # NDF
    NDF_double_gauss_fit =  double_gauss_fit.GetNDF()
    graph.Draw("AP")
    # probabilty
    Prob_double_gauss_fit = double_gauss_fit.GetProb()
    graph.Draw("AP")
    # ratio (chi2/NDF)
    chi_NDF_double_gauss_fit= chi2_double_gauss_fit/NDF_double_gauss_fit#ratio chi/NDF 
    graph.Draw("AP")
    #round the chi2 or chi2/ndf
    round_chi_double_gauss_fit = round(chi2_double_gauss_fit,2)#round off the                                         
    graph.Draw("AP")
    # round off for the legend                                

    round_chi_NDF_double_gauss_fit = round(chi_NDF_double_gauss_fit,2)#         
    graph.Draw("AP")
    # setting colors
    graph.SetLineColor( 38 )# line color
    graph.SetMarkerColor( 4 )# maker color
    graph.SetMarkerStyle( 20 )# style
    graph.SetMarkerSize( 1.1 )# size
    graph.Draw("AP")
    # colors                                                                                                                      
    canvas_name.SetFillColor(18)

    canvas_name.SetGrid()

    canvas_name.SetGridx()

    canvas_name.SetGridy()

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
        #print chi2_uncorrected,"chi2 uncorrected"


    test = ROOT.TGraphErrors(len(separation), separation,subtraction_lumi,np.zeros(separation_points), error)
    
    c = ROOT.TCanvas('fitVsdata', 'fitVsdata')
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
    legend.AddEntry(0,'#chi^{2}/ndf='+str(round_chi_NDF_double_gauss_fit), '')
    legend.SetTextSize(0.055)



    label= ROOT.TLatex()
    label.SetTextSize(0.058)
    label.DrawLatexNDC(0.75,0.90,"#bf{Double Gauss Fit}")# horizontal, vertical
    label.SetTextColor(4)

    legend.Draw()



    c.cd(2)
    test.Draw("AP")

    test.SetLineColor( 38 )
    test.SetMarkerColor( 4 )
    test.SetMarkerStyle( 20 )
    test.SetMarkerSize( 1.3 )


    test.GetXaxis().SetTitle("Separation[mm]");
    test.GetYaxis().SetTitle("lumi[data]-lumi[fit]");
    test.GetXaxis().SetTitleSize(0.064)
    test.GetYaxis().SetTitleSize(0.064)
    

    line = TLine(-0.037,0,0.037,0.0)
    line.SetLineColor(2)
    line.SetLineWidth(1)
    line.Draw()




    fileName = path + str(run_number) +"_"+str(scan)+"_"

    if Xscan== True:
        fileName += "x_scan_comp.png"
    else:
        fileName += "y_scan_comp.png"

    c.Print(fileName)
    graph.Draw("AP")

    fileName = path + str(run_number) +"_"+str(scan)+"_"
    if Xscan == True:
        fileName = fileName +"x_scan.png"
    else:
        fileName = fileName+ "y_scan.png"
        
    canvas_name.Print(fileName)





#     # covarinace matrix for the off diagnol elements
#     cov= status.GetCovarianceMatrix ()
#     corr= status.GetCorrelationMatrix ()

#     ### elements of Covariance matrix
#     e01 = cov(0,1)#1
#     e02 = cov(0,2)#2
#     e03 = cov(0,3)#3
#     e04 = cov(0,4)#4
#     e05 = cov(0,5)#5
#     e12 = cov(1,2)#6
#     e13 = cov(1,3)#7
#     e14 = cov(1,4)#8
#     e15 = cov(1,5)#9
#     e23 = cov(2,3)#10
#     e24 = cov(2,4)#11
#     e25 = cov(2,5)#12
#     e34 = cov(3,4)#13
#     e35 = cov(3,5)#14
#     e45 = cov(4,5)#15

#     # relative errors of paramters
#     # A1
#     A1= double_gauss_fit.GetParameter(0)
#     error_a1=double_gauss_fit.GetParError(0)
#     relative_a1 = (error_a1/A1)* 100

#     # A2
#     A2=double_gauss_fit.GetParameter(3)
#     error_a2=double_gauss_fit.GetParError(3)
#     relative_a2 = (error_a2/A2)* 100
    
#     #M1
#     M1= double_gauss_fit.GetParameter(1)
#     error_m1= double_gauss_fit.GetParError(1)
#     relative_m1= error_m1

#     # M2
#     M2=double_gauss_fit.GetParameter(5)
#     error_m2= double_gauss_fit.GetParError(5)
#     relative_m2= error_m2


#     # S1
#     S1= double_gauss_fit.GetParameter(2)
#     error_s1= double_gauss_fit.GetParError(2)
#     relative_s1= (error_s1/S1)*100
    
#     # S2

#     S2= double_gauss_fit.GetParameter(4)
#     error_s2= double_gauss_fit.GetParError(4)
#     relative_s2= (error_s2/S2)*100
#     print relative_s2,"s2 relative "
               







    
# #canvas
#     # sigma
#     def sigma_fit():

#         peak = double_gauss_fit.GetMaximum()#max
#  #       print peak,"peak"
#         Xmax=  double_gauss_fit.GetMaximumX()
# #        print Xmax,"xmax"
#         # parameters(derivatives)
#         # a1
#         der_a1= (1.0 / (math.sqrt(2 * math.pi) * S1) ) * np.exp( - (1.0 / 2.0) * ( ( Xmax - M1)**2 / (S1 **2) ) )# checked
#         #a2
#         der_a2= (1.0 / ( math.sqrt(2 * math.pi) * S2) ) * np.exp( - (1.0 / 2.0) * ( ( Xmax - M2)**2 / (S2 **2) ) )# checked
#         #s1
#         der_s1= - A1 * (1.0/ ( math.sqrt(2 * math.pi) * S1 **2) ) * np.exp( - (1.0 / 2.0) * ( ( Xmax - M1)**2 / (S1 **2) ) ) + A1 * (1.0/ ( math.sqrt(2 * math.pi) * S1 **4 ) )* (Xmax - M1)**2 * np.exp( - (1.0 / 2.0) * ( ( Xmax - M1)**2 / (S1 **2) ) )# checked
#         #s2
#         der_s2= - A2 * (1.0/ ( math.sqrt(2 * math.pi) * S2 **2 ) ) * np.exp( - (1.0 / 2.0) * ( ( Xmax - M2)**2 / (S2 **2) ) ) + A2 * (1.0/ ( math.sqrt(2 * math.pi) * S2 **4 ) )  * (Xmax - M2)**2 * np.exp( - (1.0 / 2.0) * ( ( Xmax - M2)**2 / (S2 **2) ) )# checked
#         #m1
#         der_m1=  A1 * (Xmax - M1) * ( 1.0/ ( math.sqrt(2 * math.pi) * S1 **3) ) * np.exp( - (1.0 / 2.0) * ( ( Xmax - M1) **2 / (S1 **2) ) )

#         #m2
#         der_m2=  A2 * (Xmax - M2) * ( 1.0/ ( math.sqrt(2 * math.pi) * S2 **3) ) * np.exp( - (1.0 / 2.0) * ( ( Xmax - M2) **2 / (S2 **2) ) )
       

#        #off- dignol elements

#         #0
#         off_diagnol_0= 2.0 * (der_a1) * (der_m1) * (e01) + 2.0 * (der_a1) * (der_s1) * (e02) + 2.0 * (der_a1) * (der_a2) * (e03) + 2.0 * (der_a1) * (der_s2) * (e04) +2.0 * (der_a1) * (der_m2) * (e05)
#         #1
#         off_diagnol_1 = 2.0 * (der_m1) * (der_s1) * e12 + 2.0 * (der_m1) * (der_a2) * e13 + 2.0 * (der_m1) * (der_s2) * e14 + 2.0 * (der_m1) * (der_m2) * e15
#         #2
#         off_diagnol_2 = 2.0 * (der_s1) * (der_a2) * e23 + 2.0 * (der_s1) * (der_s2) * e24 + 2.0 * (der_s1) * (der_m2) * e25
#         #3
#         off_diagnol_3 = 2.0 * (der_a2) * (der_s2) * e34 + 2.0 * (der_a2) * (der_m2)*e35
#         #4
#         off_diagnol_4 = 2.0 * (der_s2) * (der_m2) * e45

#         #unceratnity in peak

#         uncertanity_peak = math.sqrt( (error_a1 * der_a1)**2 + (error_a2 * der_a2)**2 + (error_s1 * der_s1)**2 + (error_s2 * der_s2)**2 + (error_m1 * der_m1)**2 + (error_m2* der_m2)**2 + off_diagnol_0 + off_diagnol_1 + off_diagnol_2 + off_diagnol_3 + off_diagnol_4 )


#         #sigma potion

#         intlimit = 0.1


#         #integral

#         integral = double_gauss_fit.Integral(-intlimit,intlimit)
#         #error on integral

#         error_integral = double_gauss_fit.IntegralError(-intlimit,intlimit)

#         #sigma
#         sigma = (1 / math.sqrt (2 * math.pi)) * integral /peak 

#         #
#         error_sigma_without=  math.sqrt( (error_integral/integral)**2 + (uncertanity_peak / peak)**2)

        
#         abs_error= sigma * error_sigma_without
#         #print abs_error,"abs_error without cov matrix" 
#         #print sigma,"sigma"




#         #derivatives of sigma

#        #  #amplitude1
#         der_sigma_a1 = (1.0 - ((A1+A2)/peak) * der_a1) / peak
#         # print A1,"A1"
#         # print A2,"A2"
#         # print peak,"peak"
#         # print der_a1,"der_a1"
# #        print der_sigma_a1,"der_sigma_a1"



# #         #amplitude 2
#         der_sigma_a2 = (1.0 - ((A1+A2)/peak) * der_a2) / peak / math.sqrt (2 * math.pi)
#         # print der_a2,"der_a2"
#         # print A1,"A1"
#         # print A2,"A2"
#         # print der_a2,"der_a2"
#         # print peak,"peak"
# #        print der_sigma_a2,"der_sigma_a2"
# #         #sigma1
       
#         der_sigma_s1 = -( (A1+A2) * der_s1) / peak**2 / math.sqrt (2 * math.pi)
#         print A1,"A1"
#         print A2,"A2"
#         print der_s1,"der-s1"
#         print peak,"peak"
#         print der_sigma_s1,"der_sigma_s1"

# #         #sigma2

# #
#         der_sigma_s2 = -( (A1+A2) * der_s2 ) / peak**2 / math.sqrt (2 * math.pi)
#         #print der_s2,"s2"
#   #      print der_sigma_s2,"der_sigma_s2"


# #         #mean 1
#         der_sigma_m1 = 2.0 *( (A1+A2) * der_m1 ) /peak**2 / math.sqrt (2 * math.pi)
#   #      print der_sigma_m1,"der_sigma_m1"
# #        print der_m1,"der_m1"
        

# #         #mean 2

#         der_sigma_m2 = 2.0 *( (A1+A2) * der_m2) / peak**2 / math.sqrt (2 * math.pi)
#  #       print der_m2,"der_m2"
#   #      print der_sigma_m2,"der_sigma_m2"
        
#         #offdiagnol terms

#         # 0
        # print "------------------------------"
        # print 2.0 * (der_sigma_a1) * (der_sigma_m1) * (e01),"1st term 2.0 * (der_sigma_a1) * (der_sigma_m1) * (e01)"
        # print 2.0 * (der_sigma_a1) * (der_sigma_s1) * (e02),"1st term 2.0 * (der_sigma_a1) * (der_sigma_s1) * (e02)"
        # print 2.0 * (der_sigma_a1) * (der_sigma_a2) * (e03),"1st term 2.0 * (der_sigma_a1) * (der_sigma_a2) * (e03)"
        # print 2.0 * (der_sigma_a1) * (der_sigma_s2) * (e04),"1st term 2.0 * (der_sigma_a1) * (der_sigma_s2) * (e04)"
        # print 2.0 * (der_sigma_a1) * (der_sigma_m2) * (e05),"1st term 2.0 * (der_sigma_a1) * (der_sigma_m2) * (e05)"

        # print (der_sigma_a1),"(der_sigma_a1)"
        # print (der_sigma_m2),"(der_sigma_m2)"
        # print e05,"e05"

        # print "-------------------------------"

#         off_sigma_0= 2.0 * (der_sigma_a1) * (der_sigma_m1) * (e01) + 2.0 * (der_sigma_a1) * (der_sigma_s1) * (e02) + 2.0 * (der_sigma_a1) * (der_sigma_a2) * (e03) + 2.0 * (der_sigma_a1) * (der_sigma_s2) * (e04) + 2.0 * (der_sigma_a1) * (der_sigma_m2) * (e05)
# #        print off_sigma_0,"off_sigma_0"




#    #     print off_sigma_0,"off_sigma_0"

# #         #1

#         off_sigma_1 = 2.0 * (der_sigma_m1) * (der_sigma_s1) * e12 + 2.0 * (der_sigma_m1) * (der_sigma_a2)*e13 + 2.0 * (der_sigma_m1) * (der_sigma_s2) * e14 + 2.0 * (der_sigma_m1) * (der_sigma_m2) * e15



#     #    print off_sigma_1,"off _sigma_1"
# #         #2
#         off_sigma_2 = 2.0 * (der_sigma_s1) * (der_sigma_a2) * e23 + 2.0 * (der_sigma_s1) * (der_sigma_s2) * e24 + 2.0 * (der_sigma_s1) * (der_sigma_m2) * e25

#      #   print off_sigma_2,"off_sigma_2"

# #         #3
#         off_sigma_3 = 2.0 * (der_sigma_a2) * (der_sigma_s2) * e34 + 2.0 * (der_sigma_a2) * (der_sigma_m2) * e35
#       #  print off_sigma_3,"off_sigma_3"
# #         #4
#         off_sigma_4 = 2 * (der_sigma_s2) * (der_sigma_m2) * e45
#        # print off_sigma_4,"off_sigma_4,"



#         sigma_uncertanity= math.sqrt( ( error_a1 * der_sigma_a1 )**2 + ( error_a2 * der_sigma_a2 )**2 + ( error_s1 * der_sigma_s1 )**2 + ( error_s2 * der_sigma_s2 )**2 + ( error_m1 * der_sigma_m1 )**2 +( error_m2 * der_sigma_m2 )**2 + off_sigma_0 + off_sigma_1 + off_sigma_2 + off_sigma_3+ off_sigma_4 )

        
#         print sigma_uncertanity,"with covarinace matrix"

#         print abs_error,"abs_error"



#         return sigma, peak, uncertanity_peak, sigma_uncertanity

#     sigma_fit = sigma_fit()

#     sigma = sigma_fit[0]
#     peak =  sigma_fit[1]

#     peak_error = sigma_fit[2]
#     sigma_uncertanity = sigma_fit[3]



    
    # path


    return  chi_NDF_double_gauss_fit, peak, chi2_double_gauss_fit#, peak_error, overall_status,sigma_uncertanity 

#my selection plots

# chi_NDF_early_x


# lumisoity for early and late
date_early = np.array([])
date_late = np.array([])
date_fom = np.array([])
# expected lumi




FOM_early= np.array([])
# early
FOM_early_low =np.array([])
FOM_early_high =np.array([])

Mu_value_early=np.array([])
#late
Mu_value_early_low=np.array([])
Mu_value_early_high=np.array([])
#late
Mu_value_late_low=np.array([])
Mu_value_late_high=np.array([])
#late
FOM_late_low =np.array([])
FOM_late_high =np.array([])

peak_errors= np.array([])
peak_errors_late= np.array([])
Mu_value_late=np.array([])

FOM_late= np.array([])
FOM=  np.array([])
intensity =np.array([])
peak_early = np.array([])
# np arrays for chi_NDF
chi_NDF_early_x = np.array([])
chi_early_x = np.array([])
chi_NDF_early_y = np.array([])
chi_NDF_late_x = np.array([])
chi_NDF_late_y = np.array([])
errors_luminosity= np.array([])
errors_luminosity_late= np.array([])

# np arrays for sigma
sigma_values_y = np.array([])
sigma_values_x_late = np.array([])
sigma_errors_x_late = np.array([])
sigma_errors_y_late = np.array([])
sigma_values_y_early = np.array([])
sigma_values_y_late = np.array([])




sigma_errors_x= np.array([])
sigma_errors_y= np.array([])

sigma_errors_y_early= np.array([])
error_on_fom = np.array([])

# runs
run_number_early = np.array([])
run_number_late = np.array([])
run_number=  np.array([])
bunches= np.array([])


####selection criteria
chi_NDF_x= np.array([])
chi_NDF_y= np.array([])
peak_error_x= np.array([])
peak_error_y= np.array([])

FOM == np.array([])

run_early_x= np.array([])
run_early_y= np.array([])
run_late_x=np.array([])
run_late_y=np.array([])


sigma1_x=np.array([])
sigma1_y=np.array([])


sigma2_x=np.array([])
sigma2_y=np.array([])


mean1_x=np.array([])
mean1_y=np.array([])


mean2_x=np.array([])
mean2_y=np.array([])
rel_FOM = np.array([])


relative_error_lumi=np.array([])
peak=np.array([])

chi_NDF_x_5= np.array([])
chi_NDF_y_5= np.array([])
error_a1=np.array([])

run= np.array([])


# fucntion for histogram

def histogram(bin,low,up,fill,Name,hist_path,hist_title):

    hist = ROOT.TH1F( "Name", "Name", bin, low, up)

#    nbin_range = hist.GetXaxis().FindBin(up) - hist.GetXaxis().FindBin(low)
#    print nbin_range,"nbin range"
#    entries=hist.Integral(0,20)
#    print entries,"entries"
#    big_bin=hist.GetMaximumBin()
#    total_entries=hist.GetEntries()
#    std=hist.GetStdDev()
#    print std 
    for i in fill:
        hist.Fill(i)

    canvas_hist= ROOT.TCanvas("Name", "Name", 800, 600)
    hist.GetXaxis().SetTitle(hist_title)
    hist.Draw()
    hist.SetFillColor(30)
    hist.SetLineColor(kRed)
    hist.SetFillStyle(3144)

    max_bin=hist.GetBinContent(hist.GetNbinsX()+1)
    
    min_bin=hist.GetBinContent(0)


    hist_leg= TLegend(0.9,0.8,0.7,0.7)#"Name")

 #   hist_leg.AddEntry(0,'Total_Entries'+str(total_entries), '')
 #   hist_leg.AddEntry(0,'Entries'+str(entries), '')
 #   hist_leg.AddEntry(0,'standard deviation'+str(std), '')
    hist_leg.AddEntry(0,'Overflow='+str(max_bin), '')

    hist_leg.AddEntry(0,'Underflow='+str(min_bin), '')
 #   hist_leg.AddEntry(0,'Max bin'+str(big_bin), '')
    hist_leg.Draw()
    

    canvas_hist.Print(hist_path)







if fit_pickle_trk:

    for filename in fit_pickle_trk:
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

            #print fom_alex,"fom Alex"
                result_x = analyse_scan(i_scan, sep_x, exp_lumi_x ,err_x ,"Emittance_scan/", run, bunches, curr_0, True, early,date )
                result_y = analyse_scan(i_scan, sep_y, exp_lumi_y ,err_y ,"Emittance_scan/", run, bunches, curr_1, False, early,date)

                n1= np.mean(output['curr'][0])/output['bunches']
                n2= np.mean(output['curr'][1])/output['bunches']
            # product of intensities
                product = n1*n2
         
                # if result_x[5]==0:# status cut (coverge status of fit)
                #     if output['bunches']>2000:# bunches should be greter then 2000

                #         if result_x[1] < 4:
                            
                #             rel_error_sigma_x = (result_x[6]/result_x[0]) *100
                #             sigma_errors_x= np.append(sigma_errors_x,rel_error_sigma_x)
                #    #     print sigma_errors_x,"sigma_errors_x"
                #             histogram(100,0.0010,10,sigma_errors_x,"sigma_error_x","myselection/der_sigma_a1_x.png","relative_error_on_width_x[%]")
                #         if result_y[1] < 4:

                #             rel_error_sigma_y  = (result_y[6]/result_y[0]) *100
                #             sigma_values_y= np.append(sigma_values_y,rel_error_sigma_y )

                #  #       print result_y[6]," y error"
                #  #       print result_y[0],"sigma values"
                #  #       print sigma_values_y,"sigma_values_y"

                #             histogram(100,0.0010,10,sigma_values_y,"sigma_error_y","myselection/der_sigma_a1_y.png","relative_error_on_width_y[%]")

                        
                   #     if result_x[1]<7:
#                        print result_x[1],"chi/NDF for x scan"



                        # sigma1_x=np.append(sigma1_x,result_x[7])
                        # histogram(100,0,12,sigma1_x,"sigma1_x","myselection/hist_sigma1x.png", "width1_x[%]")
                            
                        # sigma2_x=np.append(sigma2_x,result_x[8])
                        # histogram(100,0,40,sigma2_x,"sigma2_x","myselection/hist_sigma2x.png","width2_x[%]")
                            
                        # mean1_x=np.append(mean1_x,result_x[9])
                        # histogram(100,0.000001,0.0019,mean1_x,"mean1_x","myselection/hist_mean1x.png","mean1_x")
                            
                        # mean2_x=np.append(mean2_x,result_x[10])
                        # histogram(100,0.000001,0.04,mean2_x,"mean2_x","myselection/hist_mean2x.png","mean2_x")

                            # peak x
#                         peak_x = result_x[2]
#                         peak= np.append(peak, peak_x)

#                         peak_x_error = result_x[4]
#                         relative_error_on_peak_x = (peak_x_error/peak_x)*100
#                         peak_error_x = np.append(peak_error_x,peak_x_error)
#                             #hist for peak x
# #                        histogram(100,0.07,0.2,peak_error_x,"peak_x_error","myselection/hist_peak_x_relative_error.png","Relative_error_peak_x[%]")
#                         chi_NDF_x= np.append(chi_NDF_x, result_x[1])
# #                        histogram(100,0,20,chi_NDF_x,"chi_NDF_x","myselection/hist_chiNDF_X.png","chi/NDF_x")

#                         # width of beam x


#                         cap_sigma_x = result_x[0]# width of the beam x
#                        relative_error_width_x = (result_x[6]/result_x[0]) *100
                 #       print "-------------------------------------"

                  #      print result_x[6],"error on sigma"
                   #     print result_x[0],"sigma"

                    # #    if result_y[1]<7:

                    #     #### width of beam y
#                        relative_error_width_y = (result_y[6]/result_y[0])*100
                        


                    
                        # appending
                       #  sigma1_y=np.append(sigma1_y,result_y[7])
                       #  histogram(100,0,20,sigma1_y,"sigma1_y","myselection/hist_sigma1y.png","width1_y[%]")
                       #  # appeneding
                       #  sigma2_y=np.append(sigma2_y,result_y[8])
                       #  histogram(100,0,20,sigma2_y,"sigma2_y","myselection/hist_sigma2y.png","width2_y[%]")

                       #  # mean 1

                       #  mean1_y=np.append(mean1_y,result_y[9])
                       #  histogram(100,0.00001,0.0060,mean1_y,"mean1_y","myselection/hist_mean1y.png","mean1_y")
                                                
                       # # mean2
                       #  mean2_y=np.append(mean2_y,result_y[10])
                       #  histogram(100,0.00001,0.04,mean2_y,"mean2_y","myselection/hist_mean2y.png","mean2_y")
                                            
                        #peak_ y
 #                        peak_y = result_y[2]
#                         peak_y_error = result_y[4]
#                         relative_error_on_peak_y = (peak_y_error/peak_y)*100
#                         peak_error_y = np.append(peak_error_y,relative_error_on_peak_y)

#                         #hist for peak y
#  #                       histogram(100,0.07,0.2,peak_error_y,"peak_y_error","myselection/hist_peak_y_relative_error.png","Relative_error_peak_y[%]")


#                         # luminosity x scan
#                         expected_lumi= expected_luminosity(result_x[0],result_y[0], product)
                       
#                         abs_lumi = math.sqrt((result_x[6]/result_x[0])**2 + (result_y[6]/result_y[0])**2)
#                         lumi_error = abs_lumi * expected_lumi
#  # FOM
                        

#                         realtive_luminosity_error= (lumi_error/expected_lumi)*100
#                         relative_error_lumi= np.append(relative_error_lumi,realtive_luminosity_error)
# #                        histogram(100,0.0,10,relative_error_lumi,"lumi","myselection/hist_relative_lumi.png","relative luminosity error[%]")


#                         # FOM
#                         Figure_of_merit = peak_x/expected_lumi

#                         # errror on FOM
#                         abs_fom = math.sqrt( (lumi_error/expected_lumi)**2 + (peak_x_error /peak_x)**2)
#                         error_fom = abs_fom * Figure_of_merit

#                         # relative error on FOM
#                         relative_error_on_FOM = (error_fom/Figure_of_merit) *100
#                         # appending
#                         rel_FOM = np.append(rel_FOM , relative_error_on_FOM)

#                         # hsitogram for relative erro on FOM

#  #                       histogram(100,0.0,10,rel_FOM,"FOM","myselection/hist_realtive_FOM.png","relative FOM error[%]")

#                         FOM= np.append(FOM,Figure_of_merit)

#                         # hist for FOM
#   #                      histogram(100,0.96,1.08,FOM,"Figure_of_merit","myselection/hist_FOM.png","FOM")
#                         print len(FOM),"FOM"
#                         # chi x scan
                        
#                         chi_NDF_y= np.append(chi_NDF_y, result_y[1])
                        
#                         # hist for chi_NDF_y   
#   #                      histogram(100,0,12,chi_NDF_y,"chi_NDF_y","myselection/hist_chiNDF_Y.png","chi/NDF_y")









# plot_graph_erro("myselection/peakvsfom.png", "peakvsfom", len(FOM),peak,FOM,peak_error_x,rel_FOM,"peak","fom")

# #
#plot_graph_error("error/FOM.png", "fom", , run,fom,error_fom,"run","fom")




   #                      if rel_early<5:

   #                          bunches = np.append(bunches,output['bunches'])
   # #                     print len(bunches),"bunche----------------------------"



                        




   #                      sigma_errors_y_early = np.append(sigma_errors_y_early,result_y[5])

   #                      errors_luminosity = np.append(errors_luminosity, error_on_lumi)

     
                        
   #                      Mu_value_early= np.append(Mu_value_early,result_x[2])

   #                      chi_NDF_early_x = np.append(chi_NDF_early_x , result_x[1])

   #                      chi_NDF_early_y = np.append(chi_NDF_early_y , result_y[1])

   #                      error_on_fom = np.append(error_on_fom, error_fom)

   #                      run_number_early = np.append(run_number_early,output['run'])

   #                      FOM_early = np.append(FOM_early,Figure_of_merit_early)


   #              if output['early'] == False:


                                
                    

   #                  mu_peak_late = result_x[2]



   #                  Figure_of_merit_late = mu_peak_late/expected_lumi_late
                

                    


   #                  # filling histogram
   #                  hist_fom_late.Fill(Figure_of_merit_late)
            
   #              # sigma


   #                  abs_lumi_late = math.sqrt( (result_x[5]/ result_x[0])**2 + (result_y[5]/ result_y[0])**2)
   #                  luminosity_error_late = abs_lumi_late * expected_lumi_late



   #                  abs_fom_late = math.sqrt(( result_x[4]/result_x[2] )**2 + (luminosity_error_late / expected_lumi_late)**2)
   #                  fom_late_error =abs_fom_late * Figure_of_merit_late


   #                  relative_fom_late =(fom_late_error/ Figure_of_merit_late)*100




   #                  if relative_fom_late<5:

                        
   #                      expected_luminosity_late = np.append(expected_luminosity_late,expected_lumi_late)

   #                      errors_luminosity_late=np.append(errors_luminosity_late,luminosity_error_late)
   #                      Mu_value_late= np.append(Mu_value_late,result_x[2])
                        
   #                      peak_errors_late =np.append(peak_errors_late, result_x[4])
   #                      chi_NDF_late_x = np.append(chi_NDF_late_x , result_x[1])
   #                      chi_NDF_late_y = np.append(chi_NDF_late_y , result_y[1])
   #                      sigma_values_y_late= np.append(sigma_values_y_late,result_y[0])
   #                      sigma_errors_x_late= np.append(sigma_errors_x_late,result_x[5])
   #                      sigma_errors_y_late= np.append(sigma_errors_y_late,result_y[5])

   #                      sigma_values_x_late = np.append(sigma_values_x_late,result_x[0])


                    
   #                      late_fom_error = np.append(late_fom_error,fom_late_error)
   #                      run_number_late = np.append(run_number_late,output['run'])
   #                      FOM_late = np.append(FOM_late,Figure_of_merit_late)




   #                  hist_fom_late_error.Fill(relative_fom_late)

   #                  date_late= np.append(date_late, output['date'])            
                    
                    


#         output = dict()
        
#         output['early_run']=run_number_early
#         output['early_fom'] = FOM_early
#         output['early_fom_errors'] =error_on_fom
#         output['early_sigma_x']=sigma_values_x_early
#         output['early_sigma_y']=sigma_values_y_early
#         output['early_sigma_x_error']=sigma_errors_x_early
#         output['early_sigma_y_error']=sigma_errors_y_early 
#         output['early_lumi']= expected_luminosity_early
#         output['early_peak']=Mu_value_early
#         output['early_lumi_error']=errors_luminosity 
#         output['early_peak_error']=peak_errors
#         result.append(output)
        
# #        print result,"result"
# pickle.dump(result, open("WP_output/tightmod.pickle",'wb'))
# ROOT.TGaxis.SetMaxDigits(6)












                        
# ROOT.TGaxis.SetMaxDigits(6)
                        

# error_fom= np.concatenate([late_fom_error,error_on_fom])
# fom=np.concatenate([FOM_early,FOM_late])
# run=np.concatenate([run_number_early,run_number_late])
# sigma_x= np.concatenate([sigma_values_x_early,sigma_values_x_late])
# sigma_x_error=np.concatenate([sigma_errors_x_early,sigma_errors_x_late])



# no_sigma_x =len(sigma_x)
# plot_graph_error("error/sigma_x.png", "sigma_x", no_sigma_x, run,sigma_x,sigma_x_error,"run","sigma_x")







# sigma_y= np.concatenate([sigma_values_y_early,sigma_values_y_late])
# sigma_y_error=np.concatenate([sigma_errors_y_early,sigma_errors_y_late])


# no_sigma_y =len(sigma_y)
# plot_graph_error("error/sigma_y.png", "sigma_y", no_sigma_y, run,sigma_y,sigma_y_error,"run","sigma_y")


# chi2_x=np.concatenate([chi_NDF_early_x, chi_NDF_late_x])
# no_chi2_x =len(chi2_x)
# plot_graph("error/chi2_x.png", "chi2",no_chi2_x, run ,chi2_x,"run","chi2_x/NDF")


# chi2_y=np.concatenate([chi_NDF_early_y, chi_NDF_late_y])

# no_chi2_y =len(chi2_y)
# plot_graph("error/chi2_y.png", "chi2",no_chi2_y, run ,chi2_y,"run","chi2_y/NDF")

# # no_bunches =len(bunches)
# # plot_graph("error/bunchesvsrun.png", "bunches",no_bunches, run ,bunches,"run","bunches")



# peak_error = np.concatenate([peak_errors,peak_errors_late])
# peak=np.concatenate([Mu_value_early,Mu_value_late])
# no_peak =len(peak)
# plot_graph_error("error/peak.png", "peak_y", no_peak, run,peak,peak_error,"run","peak")


# luminosity=np.concatenate([expected_luminosity_early,expected_luminosity_late])
# luminosity_error=np.concatenate([errors_luminosity,errors_luminosity_late])
# no_luminosity =len(luminosity)
# plot_graph_error("error/luminosity.png", "luminosity", no_luminosity, run,luminosity,luminosity_error,"run","luminosity")




# no_of_run =len(run)

# no_of_run =len(run)





# no_of_sigma_x_early =len(sigma_values_x_early)
# plot_graph_error("error/sigma_values_x_early.png", "sigma_x_early", no_of_sigma_x_early, run_number_early,sigma_values_x_early,sigma_errors_x_early,"run","sigma_x")


# no_of_sigma_y_early =len(sigma_values_y_early)
# plot_graph_error("error/sigma_values_y_early.png", "sigma_y_early", no_of_sigma_y_early, run_number_early,sigma_values_y_early,sigma_errors_y_early,"run","sigma")


# no_of_lumi_x_early =len(expected_luminosity_early)
# plot_graph_error("error/lumi_valuesearly.png", "lumiearly", no_of_lumi_x_early, run_number_early,expected_luminosity_early,errors_luminosity,"run","lumi_errors")
# no_of_peak_errors_early =len(peak_errors)
# plot_graph_error("error/peak.png", "peak", no_of_peak_errors_early, run_number_early, Mu_value_early ,peak_errors,"run","mu")





                
# # 
# # #1 root canvas for Chi_NDF_early_x





# print hist_fom_early.GetStdDev(),"std of his"



# # fom late
# canvas_fom_late= ROOT.TCanvas("fom_late", "fom_late", 800, 600)
# hist_fom_late.Draw()
# hist_fom_late.SetFillColor(30)
# hist_fom_late.SetLineColor(kRed)
# hist_fom_late.SetFillStyle(3144)
# canvas_fom_late.Print("fom/hist_fom_late.png")
# #
# canvas_chi_early_y= ROOT.TCanvas("chi_NDF_early_x", "chi_NDF_early_x", 800, 600)
# hist_chi_early_y.Draw()
# 
# hist_chi_early_y.SetLineColor(kRed)
# hist_chi_early_y.SetFillStyle(3144)
# canvas_chi_early_y.Print("error/chi_early_y.png")

# # #


#hist_chi_NDF_x






# ##2 root canvas for Chi_NDF_early_y
# canvas_chi_NDF_early_y= ROOT.TCanvas("chi_NDF_early_y", "chi_NDF_early_y", 800, 600)
# hist_chi_NDF_early_y.Draw()
# hist_chi_NDF_early_y.SetFillColor(30)
# hist_chi_NDF_early_y.SetLineColor(kRed)
# hist_chi_NDF_early_y.SetFillStyle(3144)
# canvas_chi_NDF_early_y.Print("chi_NDF/chi_NDF_early_y.png")


# ## 3 root canvas for Chi_NDF_late_x
# canvas_chi_NDF_late_x= ROOT.TCanvas("chi_NDF_late_x", "chi_NDF_late_x", 800, 600)
# hist_chi_NDF_late_x.Draw()
# hist_chi_NDF_late_x.SetFillColor(30)
# hist_chi_NDF_late_x.SetLineColor(kRed)
# hist_chi_NDF_late_x.SetFillStyle(3144)
# canvas_chi_NDF_late_x.Print("chi_NDF/chi_NDF_late_x.png")


# #4 root canvas for Chi_NDF_late_y
# canvas_chi_NDF_late_y= ROOT.TCanvas("chi_NDF_late_y", "chi_NDF_late_y", 800, 600)
# hist_chi_NDF_late_y.Draw()
# hist_chi_NDF_late_y.SetFillColor(30)
# hist_chi_NDF_late_y.SetLineColor(kRed)
# hist_chi_NDF_late_y.SetFillStyle(3144)
# canvas_chi_NDF_late_y.Print("chi_NDF/chi_NDF_late_y.png")

#


# plot fom_early vs date_early




# no_of_date =len(date_fom)
# plot_graph("fom_date/fom_date.png", "fom_date",no_of_date, date_fom ,FOM,"date","fom",True)


# no_of_date_early =len(date_early)
# plot_graph("fom_date/fom_early_date.png", "fom_early_date",no_of_date_early, date_early ,FOM_early,"date_early","fom_early")


# # date _late
# no_of_date_late =len(date_late)
# plot_graph("fom_date/fom_date_late.png", "fom_late_date" , no_of_date_late, date_late ,FOM_late,"date_late","fom_late")


# #plot for early_x




#plot for early_y

# no_of_chi_NDF_early_y =len(chi_NDF_early_y)
# plot_graph("chi_NDF_run/chi_NDF_early_y.png", "chi_NDF_early_y" , no_of_chi_NDF_early_y, run_number_early , chi_NDF_early_y,"run","chi_NDF_early_x")
# #plot for late_x

# no_of_chi_NDF_late_x =len(chi_NDF_late_x)
# plot_graph("chi_NDF_run/chi_NDF_late_x.png", "chi_NDF_run_late_x", no_of_chi_NDF_late_x, run_number_late, chi_NDF_late_x,"run","chi_NDF_late_x")
# #plot for late_y

# no_of_chi_NDF_late_y =len(chi_NDF_late_y)
# plot_graph("chi_NDF_run/chi_NDF_late_y.png", "chi_NDF_run_late_y", no_of_chi_NDF_late_y, run_number_late, chi_NDF_late_y,"run","chi_NDF_late_y")


# # ##### run vs sigma x early
# # #1
# ROOT.TGaxis.SetMaxDigits(6)















# #2
# # run vs sigma_late_x 

# no_of_sigma_x_late =len(sigma_values_x_late)
# plot_graph("sigma/sigma_values_x_late.png", "sigma_x_late", no_of_sigma_x_late,run_number_late,sigma_values_x_late,"run","chi_NDF_late_x")

# # #3
# # # run vs sigma_early_y

# no_of_sigma_y_early =len(sigma_values_y_early)
# plot_graph ("sigma/sigma_values_y_early.png", "sigma_y_early",no_of_sigma_y_early,run_number_early,sigma_values_y_early,"run","chi_NDF_early_y")

# # #4
# # # run vs sigma_late_y 

# no_of_sigma_y_late =len(sigma_values_y_late)
# plot_graph ("sigma/sigma_values_y_late.png", "sigma_y_late",no_of_sigma_y_late,run_number_late,sigma_values_y_late,"run","chi_NDF_late_y")

# #fom early

# no_of_fom_early =len(FOM_early)
# plot_graph("fom/fom_early.png", "fom_early", no_of_fom_early, run_number_early,FOM_early,"run","fom_early")

# # late fom

# no_of_fom_late =len(FOM_late)
# plot_graph("fom/fom_late.png", "fom_late", no_of_fom_late, run_number_late ,FOM_late,"run","fom_late")


# no_Mu_value_early =len(Mu_value_early)
# plot_graph("Mu_FOM/Mu_FOM_early.png", "MU_FOM_early",no_Mu_value_early, Mu_value_early ,FOM_early,"MU_early","FOM_early")

# no_Mu_value_late =len(Mu_value_late)
# plot_graph("Mu_FOM/Mu_FOM_late.png", "MU_FOM_late",no_Mu_value_late, Mu_value_late ,FOM_late,"MU_late","FOM_late")

# no_Mu_value_late_low =len(Mu_value_late_low)
# plot_graph("Mu_FOM/Mu_value_late_low.png", "MU_FOM_late_low",no_Mu_value_late_low, Mu_value_late_low ,FOM_late_low,"Mu_low_late","FOM_low_mu_late", True)

# # early 
# no_Mu_value_early_low =len(Mu_value_early_low)
# plot_graph("Mu_FOM/Mu_value_early_low.png", "MU_FOM_early_low",no_Mu_value_early_low, Mu_value_early_low ,FOM_early_low,"Mu_low","FOM_low_mu", True)

# no_Mu_value_early_high =len(Mu_value_early_high)
# plot_graph("Mu_FOM/Mu_value_early_high.png", "MU_FOM_early_high",no_Mu_value_early_high, Mu_value_early_high ,FOM_early_high,"Mu_high","FOM_high_mu")

# #late



