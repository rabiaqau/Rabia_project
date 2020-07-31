

#Rabia shaheen
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
#import sympy as sym
import math
#from sympy import exp


# modules
from Functions import *
from InitialParameters import *


# these modules are for the derivative method
from peak_derivative import *
from sigma_derivative import *
from FOM_derivative import *


# this is for the numerical method(1 sigma up and 1 sigma down)

from error_propagation import *

from error_propagation_FOM import *






# importing pickle files


#fit_pickle = sorted(glob.glob("fits/integrated/fit_integrated_trk_run_362661.pickle"))

fit_pickle = sorted(glob.glob("fits_4WP/integrated/*trk*"))





    #function (analyze function)

    # parameters

    # separations (these are the 9 or 8 separation points for scanning)
    # luminosity (mu for the scans x or y )
    # error(error on data points)
    # run_number
    # bunches
    # current




def analyse_scan(scan,separations, luminosity, error, path, run_number, bunches, current, Xscan,early,date):




    # TGraph(A Graph is a graphics object made of two arrays X and Y with npoints each.)
    graph = ROOT.TGraph(len(separations),separations,luminosity)



    #TGraphErrors(A TGraphErrors is a TGraph with error bars.)
    graph = ROOT.TGraphErrors(len(separations), separations, luminosity, np.zeros(len(separations)), error ) 


    #TCanvas
    canvas = ROOT.TCanvas( 'scan_name', 'scan_name',800, 600 )



    graph_cosmetics(graph,4,4,"separations[mm]","mu")




    #Single gauss Fit (3 paramters) 
    single_gauss_fit = ROOT.TF1( "single_gauss_fit", 'gaus', -0.05, 0.05 )# single gauss (3 parameters)


    #Double gauss Fit(6 parameters)
    double_gauss_fit = ROOT.TF1("double_gauss_fit", one_dimension_double_gauss_fit, -0.05, 0.05, 6)#2d gauss (6 parameters)





   # fitting the graph with single gauss fit
    graph.Fit('single_gauss_fit','Q')


    # constant of the single gaussian function 
    # paramter 0 will represent the constant


    constant = single_gauss_fit.GetParameter(0)


    # mean of the single gaussian function
    Mean = single_gauss_fit.GetParameter(1)
#    print bunches,"bunches"

#    print bunches,"bunches"
#    print Mean,"Mean" 

    # sigma of the single gaussian function

    Sigma = single_gauss_fit.GetParameter(2)

#    print Sigma,"Sigma" 

    # getting value for the double gauss fit before we have only constant term

    Amplitude = constant * Sigma * math.sqrt(2*math.pi)

 #   print Amplitude*0.5,"Amplitude"

    # initial values for 1d double gaussian fit paramters
    
    # #amplitude 1          
    double_gauss_fit.SetParameter(0,0.5*Amplitude)
    
    # amplitude 2
    double_gauss_fit.SetParameter(3,0.5*Amplitude)
    
    # sigma 1
    double_gauss_fit.SetParameter(2,Sigma)

    # sigma 2
    double_gauss_fit.SetParameter(4,2*Sigma)
   
    # mean1
    double_gauss_fit.SetParameter(1,Mean)

    # mean2
    double_gauss_fit.SetParameter(5,Mean)


    # Please note that i have some special inital parameters value for the 1d double gaussian fit
    # for more information please see the file "InitialParameters.py"
    # here i am setting those special initial parameters





    output = SetInitialParameters(double_gauss_fit,run_number,Xscan,early,Amplitude,Sigma,Mean)
#    print output,"output from initial conditions"



   # fitting

    status = graph.Fit('double_gauss_fit','SNQ')#fitting of double gaussain



    # this step is for those scans where fit fails

    if int(status) ==4:

        # here i am increasing the number of iterations so that fit finds its convergence
        ROOT.TVirtualFitter.Fitter(graph).SetMaxIterations(10000)

        # refitting
        status = graph.Fit('double_gauss_fit','SNQ')#fitting of double gaussain

        # refitting
        if int(status)==4:
            status = graph.Fit('double_gauss_fit','SNQ')#fitting of double gaussain



    # status of the fit
    status_of_fit = int(status)


    # chi2 of 1d double gaussian fit

    chi2_double_gauss_fit = double_gauss_fit.GetChisquare()

    # NDF of the 1d double gaussian fit

    NDF_double_gauss_fit =  double_gauss_fit.GetNDF()

    # probabiltiy of the 1d double gaussian fit

    Prob_double_gauss_fit = double_gauss_fit.GetProb()

    # chi/NDF

    chi_NDF_double_gauss_fit= chi2_double_gauss_fit / NDF_double_gauss_fit
#    print chi_NDF_double_gauss_fit,"chi_NDF_double_gauss_fit"
    
    # round off chi2 of the 1d double gaussian fit

    round_chi_double_gauss_fit = round(chi2_double_gauss_fit,2)
    # round off chi2
    round_chi_NDF_double_gauss_fit = round(chi_NDF_double_gauss_fit,2)
    

    # setting the title, color and marker size
    graph_cosmetics(graph,4,4,"separations[mm]","mu")

    # drawing
    graph.Draw("AP")


    

    # luminsoity calculated from fit
    fit_luminosity = one_dimension_double_gauss_fit( [separations], double_gauss_fit.GetParameters())


    # empty arrays 
    chi2_1d_double_gauss_fit_Vector = np.array([])


    # (data-fit) luminosity difference vector

    luminosity_diffVector = np.array([])


    # luminoisty of data at every separation
    luminosity_dataVector=np.array([])


    # looping on all separations points
    # to get lumionsity of each separations points
    


    for separation in range(len(separations)):  

        # appending dataVector
        luminosity_dataVector =np.append(luminosity_dataVector,luminosity[separation])

        # difference of luminoisty per separation and fit
        luminosity_diff = (luminosity[separation]- fit_luminosity[separation])

        # appending luminosity difference
        luminosity_diffVector= np.append(luminosity_diffVector,luminosity_diff)


        # it will give me the chi2 at every point of separation
        
        # applying the formula of the chi2
        chi2_1d_double_gauss_fit = ((luminosity[separation]- fit_luminosity[separation])/error[separation])**2

        
        # appending the vector of the chi2 of every point of separation
        chi2_1d_double_gauss_fit_Vector = np.append(chi2_1d_double_gauss_fit_Vector,chi2_1d_double_gauss_fit)
    





    # tgraph for the difference in luminosity
    tgraphErrors_luminosity_diff = ROOT.TGraphErrors(len(separations),separations, luminosity_diffVector ,np.zeros(len(separations)), error)



    # calling function for 2 canvases(Functions.py)
    tgraph_canvas=graph_2_canvas_1(graph,tgraphErrors_luminosity_diff,"separation","mu","separation","lumi[data]-lumi[fit]")

    tgraph_canvas.Update()



    # path
    fileName = path + str(run_number) +"_"+str(scan)+"_"
    if Xscan== True:
        fileName += "x_scan_comp.png"
    else:
        fileName += "y_scan_comp.png"
    tgraph_canvas.Print(fileName)
    tgraph_canvas.Modified()



       
     # path
    graph.Draw("AP")
    fileName = path + str(run_number) +"_"+str(scan)+"_"
    if Xscan == True:
        fileName = fileName +"x_scan.png"
    else:
        fileName = fileName+ "y_scan.png"
        
    canvas.Print(fileName)
    





    # covarinace matrix

    cov= status.GetCovarianceMatrix ()                                                                               
  

    # correlaation matrix
    corr= status.GetCorrelationMatrix ()                                                                             





####### NUmerical method (1 sigma up and 1 sigma down)

    

    numerical_CapSigma_error = calc_with_error(width, double_gauss_fit, status,double_gauss_fit)

 #   print numerical_sigma_error,"numerical_sigma_error"

    numerical_peak_error = calc_with_error(peak_gauss, double_gauss_fit, status,double_gauss_fit)


#    print numerical_peak_error,"numerical_peak_error"






    ### elements of Covariance matrix            
  
    # 0 amp1
    # 1 mean1
    # 2 sigma1
    # 3 amp2
    # 4 sigma2
    # 5 mean2


    # i am defining a dictionary to use it latter for calculating the error on FOM
    
    covariance_element=dict()

    covariance_element['cov_amp1_mean1'] = cov(0,1)#1

    covariance_element['cov_amp1_sigma1'] = cov(0,2)#2

    covariance_element['cov_amp1_amp2'] = cov(0,3)#3
    
    covariance_element['cov_amp1_sigma2'] = cov(0,4)#4                                                                                                 
    covariance_element['cov_amp1_mean2'] = cov(0,5)#5                                                                                                 
    covariance_element['cov_mean1_sigma1'] = cov(1,2)#6                                                                                                 
    covariance_element['cov_mean1_amp2'] = cov(1,3)#7                                                                                                 
    covariance_element['cov_mean1_sigma2'] = cov(1,4)#8                                                                                                 
    covariance_element['cov_mean1_mean2'] = cov(1,5)#9                                                                                                 
    covariance_element['cov_sigma1_amp2'] = cov(2,3)#10                                                                                                
    covariance_element['cov_sigma1_sigma2'] = cov(2,4)#11                                                                                                
    covariance_element['cov_sigma1_mean2'] = cov(2,5)#12                                                                                                
    covariance_element['cov_amp2_sigma2'] = cov(3,4)#13                                                                                                
    covariance_element['cov_amp2_mean2'] = cov(3,5)#14                                                                                                
    covariance_element['cov_sigma2_mean2'] = cov(4,5)#15                     



                                                                           
    # gives me status of covarinace matrix
    staus_cov_matrix = status.CovMatrixStatus()




    # parameters information


    
    # amplitude1
    amplitude1_info = parameter_info(double_gauss_fit,0)
    amplitude1_value =amplitude1_info[0]
#    print amplitude1_value,"1"

    amplitude1_error =amplitude1_info[1]
#    print amplitude1_error,"amplitude error only diagnol terms"


    amplitude1_rel_error =amplitude1_info[2]



    # amplitude2
    amplitude2_info = parameter_info(double_gauss_fit,3)
    amplitude2_value =amplitude2_info[0]
    amplitude2_error =amplitude2_info[1]


    amplitude2_rel_error =amplitude2_info[2]

 #   print amplitude2_value,"2"

    #sigma1

    sigma1_info = parameter_info(double_gauss_fit,2)
    sigma1_value=sigma1_info[0]
    sigma1_error=sigma1_info[1]
    sigma1_rel_error=sigma1_info[2]


    # sigma2
    sigma2_info = parameter_info(double_gauss_fit,4)
    sigma2_value=sigma2_info[0]
    sigma2_error=sigma2_info[1]
    sigma2_rel_error=sigma2_info[2]



    # mean1
    mean1_info= parameter_info(double_gauss_fit,1)
    mean1_value=mean1_info[0]
    mean1_error=mean1_info[1]
    mean1_rel_error=mean1_info[2]




    #mean2
    mean2_info= parameter_info(double_gauss_fit,5)
    mean2_value=mean2_info[0]
    mean2_error=mean2_info[1]
    mean2_rel_error=mean2_info[2]


    # relative error of then paamters
    paramter_relative_error = dict()

    paramter_relative_error['amp1']=amplitude1_rel_error
    paramter_relative_error['amp2']=amplitude2_rel_error
    paramter_relative_error['sigma1']= sigma1_rel_error
    paramter_relative_error['sigma2']=sigma2_rel_error
    paramter_relative_error['mean1']=mean1_rel_error
    paramter_relative_error['mean2']=mean2_rel_error

    # dignol errors dictionary
    diagnol_error = dict()
    diagnol_error['error_amplitude1'] = amplitude1_error
    diagnol_error['error_amplitude2']=amplitude2_error
    diagnol_error['error_sigma1']=sigma1_error
    diagnol_error['error_sigma2']=sigma2_error
    diagnol_error['error_mean1']=mean1_error
    diagnol_error['error_mean2']=mean2_error



    # parameter values dictionary
    parameter_value = dict()
    parameter_value['value_amplitude1'] = amplitude1_value
    parameter_value['value_amplitude2']=amplitude2_value
    parameter_value['value_sigma1']=sigma1_value
    parameter_value['value_sigma2']=sigma2_value
    parameter_value['value_mean1']=mean1_value
    parameter_value['value_mean2']=mean2_value

    peak_position=  double_gauss_fit.GetMaximumX()




    #  cap sigma function

    def cap_sigma():

        # maximum of the fit fucntion
        peak = double_gauss_fit.GetMaximum()
        
        # position of the peak


        # limits for the integration
        
        intlimit = 0.1        
        
        # integral of the fit function( it gives me area under the curve)
        
        integral = double_gauss_fit.Integral(-intlimit,intlimit)
        
        # formula for finding the cap sigma

        Capsigma = (1 / math.sqrt (2 * math.pi)) * integral / peak 

        

        Xmax = double_gauss_fit.GetMaximumX()


        # module peak_derivtive(for more information)


        der_a1= (1.0 / (math.sqrt(2 * math.pi) * sigma1_value) ) * np.exp( - (1.0 / 2.0) * ( ( Xmax - mean1_value)**2 / (sigma1_value **2) ) )# checked        

         #a2                                                                                                                      

        der_a2= (1.0 / ( math.sqrt(2 * math.pi) * sigma2_value) ) * np.exp( - (1.0 / 2.0) * ( ( Xmax - mean2_value)**2 / (sigma2_value **2) ) )# \checked      

        # s1                                                                                                                      
        der_s1= - amplitude1_value * (1.0/ ( math.sqrt(2 * math.pi) * sigma1_value **2) ) * np.exp( - (1.0 / 2.0) * ( ( Xmax - mean1_value)**2 / (sigma1_value **2) ) ) + amplitude1_value * (1.0/ ( math.sqrt(2 * math.pi) * sigma1_value **4 ) )* (Xmax - mean1_value)**2 * np.exp( - (1.0 / 2.0) * ( ( Xmax - mean1_value)**2 / (sigma1_value **2) ) )
                                                                                                                                  
         #s2                                                                                                                      
        der_s2= - amplitude2_value * (1.0/ ( math.sqrt(2 * math.pi) * sigma2_value **2 ) ) * np.exp( - (1.0 / 2.0) * ( ( Xmax - mean2_value)**2 / (sigma2_value **2) ) ) + amplitude2_value* (1.0/ ( math.sqrt(2 * math.pi) * sigma2_value **4 ) )  * (Xmax - mean2_value)**2 * np.exp( - (1.0 / 2.0) * ( ( Xmax - mean2_value)**2 / (sigma2_value **2) ) )                                                                                                                               

         #m1                                                                                                                      

        der_m1=  amplitude1_value * (Xmax - mean1_value) * ( 1.0/ ( math.sqrt(2 * math.pi) * sigma1_value **3) ) * np.exp( - (1.0 / 2.0) * ( ( Xmax - mean1_value) **2 / (sigma1_value **2) ) ) # checked                                                                                                              
         #m2                                                                                                                      
        der_m2=  amplitude2_value * (Xmax - mean2_value) * ( 1.0/ ( math.sqrt(2 * math.pi) * sigma2_value **3) ) * np.exp( - (1.0 / 2.0) * ( ( Xmax - mean2_value) **2 / (sigma2_value **2) ) ) #checked                                

        

        off_diagnol_0= 2.0 * der_a1 * der_m1 * cov(0,1) + 2.0 * der_a1 * der_s1 * cov(0,2) + 2.0 * der_a1 * der_a2 * cov(0,3) + 2.0 * der_a1 * der_s2 * cov(0,4) + 2.0 * der_a1 * der_m2 * cov(0,5)



        off_diagnol_1 = 2.0 * der_m1 * der_s1 * cov(1,2) + 2.0 * der_m1 * der_a2 * cov(1,3) + 2.0 * der_m1 * der_s2 * cov(1,4) + 2.0* der_m1 * der_m2 * cov(1,5) 




         #2                                                                                                                                                                                                 
        off_diagnol_2 = 2.0 * der_s1 * der_a2 * cov(2,3) + 2.0 * der_s1 * der_s2 * cov(2,4) + 2.0 * der_s1 * der_m2 * cov(2,5)


         #3                                                                                                                                                                                                 
        off_diagnol_3 = 2.0 * der_a2 * der_s2 * cov(3,4) + 2.0 * der_a2 * der_m2 * cov(3,5)




         #4                                                                                                                                                                                                 
        off_diagnol_4 = 2.0 * der_s2 * der_m2 * cov(4,5)



        peak_error = math.sqrt( ( amplitude1_error * der_a1 )**2 + ( amplitude2_error * der_a2 ) **2 + ( sigma1_error * der_s1 )**2 + ( sigma2_error * der_s2)**2 + (mean1_error * der_m1)**2 + (mean2_error * der_m2)**2 + off_diagnol_0 + off_diagnol_1 + off_diagnol_2 + off_diagnol_3 + off_diagnol_4 )









#        derivative_peak_amplitude
        # peak_amplitude1=derivative_peak_amplitude.evalf(subs={amplitude:amplitude1_value,sigma:sigma1_value,mean:mean1_value, x:peak_position})


        # peak_amplitude2=derivative_peak_amplitude.evalf(subs={amplitude:amplitude2_value,sigma:sigma2_value,mean:mean2_value, x:peak_position})

        # peak_sigma1 =derivative_peak_sigma.evalf(subs={amplitude:amplitude1_value,sigma:sigma1_value,mean:mean1_value, x:peak_position})

        # peak_sigma2 = derivative_peak_sigma.evalf(subs={amplitude:amplitude1_value,sigma:sigma1_value,mean:mean1_value, x:peak_position})

        # peak_mean1=derivative_peak_mean.evalf(subs={amplitude:amplitude1_value,sigma:sigma1_value,mean:mean1_value, x:peak_position})

        # peak_mean2=derivative_peak_mean.evalf(subs={amplitude:amplitude1_value,sigma:sigma1_value,mean:mean1_value, x:peak_position})
        
#        dignol_peak= (peak_amplitude1* diagnol_error['error_amplitude1'] )**2 + (peak_amplitude2 *diagnol_error['error_amplitude2'])**2 + (peak_sigma1 *diagnol_error['error_sigma1'] )**2+(peak_sigma2 *diagnol_error['error_sigma2'] )**2 + (peak_mean1 *diagnol_error['error_mean1'] )**2 + (peak_mean2 *diagnol_error['error_mean2'] )**2




        

        ## here i am calling the derivative function to take derivative of peak


        # off diagnol terms
        
#        off_diagnol_peak_terms = offdiagnol(peak_amplitude1,peak_sigma1,peak_mean1,peak_amplitude2,peak_sigma2,peak_mean2,covariance_element['cov_amp1_mean1'], covariance_element['cov_amp1_sigma1'],covariance_element['cov_amp1_amp2'],covariance_element['cov_amp1_sigma2'],covariance_element['cov_amp1_mean2'],covariance_element['cov_mean1_sigma1'],covariance_element['cov_mean1_amp2'], covariance_element['cov_mean1_sigma2'],covariance_element['cov_mean1_mean2'],covariance_element['cov_sigma1_amp2'],covariance_element['cov_sigma1_sigma2'], covariance_element['cov_sigma1_mean2'],covariance_element['cov_amp2_sigma2'], covariance_element['cov_amp2_mean2'],covariance_element['cov_sigma2_mean2'] )



        # # here is the error on peak
 #       error_on_peak = math.sqrt( dignol_peak + off_diagnol_peak_terms)



#        print error_on_peak,"error on peak"


    
        # module sigma_derivtive.py

        #amplitude1_value --> 0
        #amplitude2_value --> 1
        #sigma1_value --> 2
        #sigma2_value --> 3
        #mean1_value --> 4
        #mean2_value --> 5

        # funtion for setting value for CapSigma
        
        Capsigma_derivative_value= set_Capsigma_derivative_value(amplitude1_value, amplitude2_value, sigma1_value, sigma2_value, mean1_value, mean2_value, peak_position)



        # diagnol capsigma terms
        diagnol_Capsigma_terms = diagnol(Capsigma_derivative_value[0],Capsigma_derivative_value[1],Capsigma_derivative_value[2],Capsigma_derivative_value[3],Capsigma_derivative_value[4],Capsigma_derivative_value[5],amplitude1_error,amplitude2_error,sigma1_error,sigma2_error,mean1_error,mean2_error)



        # off diagnol function for Capsigma
        off_diagnol_Capsigma_terms = offdiagnol(Capsigma_derivative_value[0],Capsigma_derivative_value[2], Capsigma_derivative_value[4],Capsigma_derivative_value[1], Capsigma_derivative_value[3], Capsigma_derivative_value[5],covariance_element['cov_amp1_mean1'], covariance_element['cov_amp1_sigma1'],covariance_element['cov_amp1_amp2'],covariance_element['cov_amp1_sigma2'],covariance_element['cov_amp1_mean2'],covariance_element['cov_mean1_sigma1'],covariance_element['cov_mean1_amp2'], covariance_element['cov_mean1_sigma2'],covariance_element['cov_mean1_mean2'],covariance_element['cov_sigma1_amp2'],covariance_element['cov_sigma1_sigma2'], covariance_element['cov_sigma1_mean2'],covariance_element['cov_amp2_sigma2'], covariance_element['cov_amp2_mean2'],covariance_element['cov_sigma2_mean2'] )


       # Cap sigma error

        
        Capsigma_error = math.sqrt(off_diagnol_Capsigma_terms + diagnol_Capsigma_terms )

        
        # dictionary for storing derivatives of Capsigma 


        #amplitude1_value --> 0
        #amplitude2_value --> 1
        #sigma1_value --> 2
        #sigma2_value --> 3
        #mean1_value --> 4
        #mean2_value --> 5




        Capsigma_derivative_dict = dict()
        Capsigma_derivative_dict['der_CapSigma_amplitude1_value'] = Capsigma_derivative_value[0]
        Capsigma_derivative_dict['der_CapSigma_amplitude2_value'] = Capsigma_derivative_value[1]
        Capsigma_derivative_dict['der_CapSigma_sigma1_value']     = Capsigma_derivative_value[2]
        Capsigma_derivative_dict['der_CapSigma_sigma2_value']     = Capsigma_derivative_value[3]
        Capsigma_derivative_dict['der_CapSigma_mean1_value']      = Capsigma_derivative_value[4]
        Capsigma_derivative_dict['der_CapSigma_mean2_value']      = Capsigma_derivative_value[5]

        



        return Capsigma, peak, peak_error, Capsigma_error, Capsigma_derivative_dict


    # calling the cap_sigma fucntion


    sigma_fit = cap_sigma()

    # sigma for luminosity
    sigma_for_luminosity = sigma_fit[0]

    # peak
    peak =  sigma_fit[1]

    # error on peak
    error_peak= sigma_fit[2]


    # error on cap sigma
    error_Capsigma=sigma_fit[3]

    # peak derivative dictionary 
    CapSigma_derivative_dict = sigma_fit[4]



    # analyze function will return these values

    # 0 --> sigma for luminosity
    # 1 --> chi/ndf_double gaussian fit
    # 2 -->peak
    # 3 --> chi2 double gauss fit
    # 4 --> status of fit
    # 5 --> status of covariance matrix
    # 6 --> peak error
    # 7 --> error for cap sigma
    # 8 --> covariance element
    # 9 --> dignol error
    # 10 -->parameter value
    # 11 -->position of maximum peak
    # 12 -->Capsigma_derivative_dict
    # 13 --> Numerical CapSigma_error
    # 14 --> Numerical peak_error
    # 15 --> double gauss fit
    # 16 --> status of fit

    # returning values analyse_scan
    return sigma_for_luminosity, chi_NDF_double_gauss_fit, peak, chi2_double_gauss_fit, status_of_fit, staus_cov_matrix , error_peak, error_Capsigma, covariance_element, diagnol_error, parameter_value, peak_position, CapSigma_derivative_dict, numerical_CapSigma_error, numerical_peak_error,double_gauss_fit,status




num_width_x_error = np.array([])
num_width_y_error = np.array([])

# peak x and y                                                                  
num_peak_x_error = np.array([])
num_peak_y_error = np.array([])

# FOM and nominal mu                                                            
num_FOM_error = np.array([])
num_nominal_mu_error =np.array([])



# derivative method


peak_x =np.array([])
peak_x_error=np.array([])

peak_y =np.array([])
peak_y_error=np.array([])






CapSigma_x =np.array([])
CapSigma_x_error=np.array([])



CapSigma_y =np.array([])
CapSigma_y_error=np.array([])






nominal_mu =np.array([])
nominal_mu_error =np.array([])



FOM =np.array([])
FOM_error =np.array([])





scans_date= np.array([])
run_number= np.array([])


chi2_NDF_x=np.array([])

chi2_NDF_y=np.array([])



# pickle file 

for filename in fit_pickle:
    with open(filename,'rb') as fit_pickle:
        scan = pickle.load(fit_pickle)
        #scans = scan['TightModLumi']

        for i_scan in range(len(scan)):
            output = scan[i_scan]
            #print i_scan,"scans"
            sep_x = output['sep_x']# 9 points
            exp_lumi_x = output['exp_lumi_x']# expected luminosity x
            err_x = output['err_x']# error on x scan
            curr_0 = output['curr'][0]# current of first beam 
            date = output['date']# date or year
            sep_y = output['sep_y']# sep for y scan
            exp_lumi_y = output['exp_lumi_y']# expected luminosity y
            err_y = output['err_y']# error on y scan
            curr_1 = output['curr'][1]# current of second beam
            run = output['run']# run number
            bunches= output['bunches']# bunches
            early = output['early']# early scan condition


            result_x = analyse_scan(i_scan, sep_x, exp_lumi_x ,err_x ,"Emittance_scan/", run, bunches, curr_0, True, early,date ) # calling the analyse function with x scan switch on
            
            result_y = analyse_scan(i_scan, sep_y, exp_lumi_y ,err_y ,"Emittance_scan/", run, bunches, curr_1, False, early,date)# calling the analyse function with x scan switch off

        





            # this part will use for determining the nominal mu
        

            # intensities

            n1= np.mean(output['curr'][0])/output['bunches']
            n2= np.mean(output['curr'][1])/output['bunches']
            
            # product of intensities
            product_of_intensity = n1*n2
            
            

#            print product_of_intensity,"product_of_intensity"

            # selection cuts


                

            
            if  output['bunches'] >2000:

                if result_x[4]==0 and result_y[4]==0:

                    if result_x[5] == 3 and result_y[5]==3: 

                        if  result_x[1] < 7 and  result_y[1] < 7:

                            
                            # arrays

                            chi2_NDF_x=np.append(chi2_NDF_x,result_x[1])



                            run_number =np.append(run_number,output['run'])



                            chi2_NDF_y=np.append(chi2_NDF_y,result_y[1])



                            scans_date = np.append(scans_date, output['date'])



                                
                            ###### fit results
                            fit_x = result_x[15]
                            fit_status_x = result_x[16]

                            fit_y = result_y[15]
                            fit_status_y = result_y[16]
                            
                            ##############################################
                            # sigma X  and Y scan(width of beam) and errors
                            ##############################################



                            CapSigma_x =np.append(CapSigma_x, result_x[0] )

                            CapSigma_x_error=np.append(CapSigma_x_error, result_x[7] )



                            CapSigma_y =np.append(CapSigma_y, result_y[0] )

                            CapSigma_y_error=np.append(CapSigma_y_error, result_y[7] )
                        


                    

                            #############################################
                            # peak X  and Y scan(width of beam) and errors
                            #############################################


                            peak_x =np.append(peak_x,result_x[2] )

                            peak_x_error=np.append( peak_x_error,result_x[6] )


                            peak_y =np.append(peak_y,result_y[2])

                            peak_y_error=np.append(peak_y_error,result_y[6])





                            

                            ######################################
                            # nominal mu (luminosity) and its err0r
                            #######################################

                        
                            luminosity_of_beam =expected_luminosity(result_x[0],result_y[0],product_of_intensity)

                            # error on lumi
                            absolute_luminosity = math.sqrt((result_y[7]/result_y[0])**2 + (result_x[7]/result_x[0])**2)

                            # error in luminsoity
                            error_on_luminosity = absolute_luminosity * luminosity_of_beam


                            nominal_mu =np.append(nominal_mu, luminosity_of_beam  )
                            nominal_mu_error =np.append(nominal_mu_error,error_on_luminosity   )


                          



                            #### numerical error on luminosity


                            lumi_x =  calc_with_error_lumi(width, fit_x, fit_status_x, result_y[0], product_of_intensity, fit_x)
 
                            lumi_y =  calc_with_error_lumi(width, fit_y, fit_status_y, result_x[0], product_of_intensity, fit_y)


                            num_lumi_error = math.sqrt(lumi_x + lumi_y)

#                            print lumi_error,"num method"





                            ###################################
                            # figure of merit
                            # presntation(8/june/2020)
                            # FOM =(mu(x-scan)/mu_nominal
                            # module FOM_derivative.py
                            #################################


                            ### numerical method for FOM error




                            fom_x = calc_with_error_fom(width,max_function, fit_x, fit_status_x, result_y[0], product_of_intensity, fit_x)
#                            print fom_x,"fom x"
                            fom_y = calc_with_error_fom_y(width, fit_y, fit_status_y,result_x[0], result_x[2],product_of_intensity,fit_y)

                            num_fom_error = math.sqrt(fom_x + fom_y)
#                            print num_fom_error,"num_fom_error"

                            # figure of merit



                            Figure_of_Merit= result_x[2] / luminosity_of_beam



                            # I have exlained the modified form of FOM defination
                            # FOM = (31.3249 * (A_1x +A_2x)* sigma_y) / n1n2

                            # x derivatives

                            dignol_offdiagnol_derivative_of_FOM_x = FOM_x(result_y[0],product_of_intensity,result_x[9]['error_amplitude1'], result_x[9]['error_amplitude2'],result_x[8]['cov_amp1_amp2'])


#                            print dignol_offdiagnol_derivative_of_FOM_x
                            # dignol terms for FOM w.r.t y scan parameters

                            diagnol_derivative_y = FOM_y(product_of_intensity,result_x[10]['value_amplitude1'], result_x[10]['value_amplitude2'], result_y[12]['der_CapSigma_amplitude1_value'], result_y[12]['der_CapSigma_amplitude2_value'], result_y[12]['der_CapSigma_sigma1_value'], result_y[12]['der_CapSigma_sigma2_value'], result_y[12]['der_CapSigma_mean1_value'] ,result_y[12]['der_CapSigma_mean2_value'],result_y[9]['error_amplitude1'],result_y[9]['error_amplitude2'], result_y[9]['error_sigma1'],result_y[9]['error_sigma2'],result_y[9]['error_mean1'],result_y[9]['error_mean2'])
                            



                            # FOM_Y FUNCTION returns me 6 values the last one correspond to dignol y term
                            dignol_y_terms = diagnol_derivative_y[6]
                        

                            # off dignol y terms
                        
                            off_dignol_y= offdiagnol(diagnol_derivative_y[0],diagnol_derivative_y[2],diagnol_derivative_y[4],diagnol_derivative_y[1],diagnol_derivative_y[3],diagnol_derivative_y[5],result_y[8]['cov_amp1_mean1'], result_y[8]['cov_amp1_sigma1'],result_y[8]['cov_amp1_amp2'],result_y[8]['cov_amp1_sigma2'],result_y[8]['cov_amp1_mean2'],result_y[8]['cov_mean1_sigma1'],result_y[8]['cov_mean1_amp2'], result_y[8]['cov_mean1_sigma2'],result_y[8]['cov_mean1_mean2'],result_y[8]['cov_sigma1_amp2'],result_y[8]['cov_sigma1_sigma2'], result_y[8]['cov_sigma1_mean2'],result_y[8]['cov_amp2_sigma2'], result_y[8]['cov_amp2_mean2'],result_y[8]['cov_sigma2_mean2'] ) 



                            # uncertnaity on FOM
                            uncertanity_on_FOM = math.sqrt(dignol_offdiagnol_derivative_of_FOM_x + dignol_y_terms + off_dignol_y)


                            ## appending

                            FOM=np.append(FOM , Figure_of_Merit)

                            FOM_error = np.append(FOM_error , uncertanity_on_FOM)





                            ## numeircal method



                            num_width_x_error = np.append(num_width_x_error , result_x[13] )
                            num_width_y_error = np.append( num_width_y_error , result_y[13])

                            # peak x and y                                                                                                                                  
                            num_peak_x_error = np.append(num_peak_x_error , result_x[14] )
                            num_peak_y_error = np.append(num_peak_y_error , result_y[14] )

                            # FOM and nominal mu                                                                                                                            
                            num_FOM_error = np.append(num_FOM_error , num_fom_error )


                            num_nominal_mu_error = np.append(num_nominal_mu_error , num_lumi_error )


                            my_results = []
               

                            num_errors = {'Num_CapSigmaX_error':num_width_x_error,'Num_CapSigma_y_error':num_width_y_error, 'Num_Peak_x_error':num_peak_x_error ,'Num_Peak_y_error':num_peak_y_error ,'Num_FOM_error':num_FOM_error,'Num_Nominal_mu_error':num_nominal_mu_error }



                            der_errors = {'figure_of_merit':FOM ,'Figure_Merit_error':FOM_error,'run':run_number,'Peak_x':peak_x,'Peak_x_error':peak_x_error,'Peak_y':peak_y,'Peak_y_error':peak_y_error,'CapSigmaX':CapSigma_x,'CapSigmaX_error':CapSigma_x_error,'CapSigmaY':CapSigma_y,'CapSigmaY_error':CapSigma_y_error,'Chi2_NDF_X':chi2_NDF_x,'Chi2_NDF_Y':chi2_NDF_y,'Date':scans_date,'nom_mu':nominal_mu ,'nom_mu_error':nominal_mu_error}


                            my_results=[num_errors,der_errors]
                            





                            with open('my_results/my_results_allScans_data.pkl','wb') as pickle_file:
                                pickle.dump(my_results,pickle_file)





                        










                        
                            

                            
                            
                            
                            






















    

