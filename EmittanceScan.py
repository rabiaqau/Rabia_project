
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
from Functions import *
from InitialParameters import *
from derivative import *
import sympy as sym
import math
from sympy import exp



# importing pickle files


fit_pickle = sorted(glob.glob("fits_4WP/integrated/fit_integrated_trk_run_349481.pickle"))

#fit_pickle = sorted(glob.glob("fits_4WP/integrated/*trk*"))







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


    # sigma of the single gaussian function

    Sigma = single_gauss_fit.GetParameter(2)



    # getting value for the double gauss fit before we have only constant term

    Amplitude = constant * Sigma * math.sqrt(2*math.pi)


    # initial values for 1d double gaussian fit paramters
    
    #amplitude 1          
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


    output=SetInitialParameters(double_gauss_fit,run_number,Xscan,early,Amplitude,Sigma,Mean)


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
    amplitude1_error =amplitude1_info[1]
    amplitude1_rel_error =amplitude1_info[2]

    # amplitude2
    amplitude2_info = parameter_info(double_gauss_fit,3)
    amplitude2_value =amplitude2_info[0]
    amplitude2_error =amplitude2_info[1]
    amplitude2_rel_error =amplitude2_info[2]



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



    amplitude1= double_gauss_fit.GetParameter(0)
    amplitude2= double_gauss_fit.GetParameter(3)
    

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

        Capsigma = (1 / math.sqrt (2 * math.pi)) * integral /peak 
        

        ## here i am calling the derivative function to take derivative of peak




         # value of amplitude 1
        Value_der_peak_amplitude1 = derivative_peak_amplitude.evalf(subs={amplitude:amplitude1_value,sigma:sigma1_value,mean:mean1_value,x:Xmax})

         #  print Value_der_peak_amplitude1,"Value_der_peak_amplitude1"




         # value of amplitude 2
        Value_der_peak_amplitude2 = derivative_peak_amplitude.evalf(subs={amplitude:amplitude2_value,sigma:sigma2_value,mean:mean2_value,x:Xmax})
#        print Value_der_peak_amplitude2,"Value_der_peak_amplitude2"





        
        # value of sigma1
        Value_der_peak_sigma1 = derivative_peak_sigma.evalf(subs={ amplitude:amplitude1_value, sigma:sigma1_value, mean:mean1_value, x:Xmax})

#        print Value_der_peak_sigma1,"Value_der_peak_sigma1"


        # sigma2

        Value_der_peak_sigma2 = derivative_peak_sigma.evalf(subs={ amplitude:amplitude2_value, sigma:sigma2_value, mean:mean2_value, x:Xmax})
#        print Value_der_peak_sigma2,"Value_der_peak_sigma2"



        # mean1

        Value_der_peak_mean1 = derivative_peak_mean.evalf(subs={ amplitude:amplitude1_value, sigma:sigma1_value, mean:mean1_value, x:Xmax})


#        print Value_der_peak_mean1,"Value_der_peak_mean1"


        # mean2

        Value_der_peak_mean2 = derivative_peak_mean.evalf(subs={ amplitude:amplitude2_value, sigma:sigma2_value, mean:mean2_value, x:Xmax})

#        print Value_der_peak_mean2,"Value_der_peak_mean2"





        # dignol elements of peak

        diagnol_peak_terms = ( amplitude1_error * Value_der_peak_amplitude1 )**2 + ( amplitude2_error * Value_der_peak_amplitude2 )**2 + ( sigma1_error * Value_der_peak_sigma1 )**2 + ( sigma2_error * Value_der_peak_sigma2 )**2 + ( mean1_error * Value_der_peak_mean1 )**2 + ( mean2_error * Value_der_peak_mean2 )**2
        

        # off diagnol function
        
        off_diagnol_peak_terms = offdiagnol(Value_der_peak_amplitude1, Value_der_peak_sigma1, Value_der_peak_mean1,Value_der_peak_amplitude2, Value_der_peak_sigma2, Value_der_peak_mean2,covariance_element['cov_amp1_mean1'], covariance_element['cov_amp1_sigma1'],covariance_element['cov_amp1_amp2'],covariance_element['cov_amp1_sigma2'],covariance_element['cov_amp1_mean2'],covariance_element['cov_mean1_sigma1'],covariance_element['cov_mean1_amp2'], covariance_element['cov_mean1_sigma2'],covariance_element['cov_mean1_mean2'],covariance_element['cov_sigma1_amp2'],covariance_element['cov_sigma1_sigma2'], covariance_element['cov_sigma1_mean2'],covariance_element['cov_amp2_sigma2'], covariance_element['cov_amp2_mean2'],covariance_element['cov_sigma2_mean2'] )



        # here is the error on peak
        error_on_peak = math.sqrt( diagnol_peak_terms  + off_diagnol_peak_terms )  


        #### now its time to store the values of all determinants of peak in a dictionary so that i can use it for 
        ####   determining the error on FOM (you will see) :)
        derivative_of_peak = dict()
        derivative_of_peak['der_peak_amplitude1'] = Value_der_peak_amplitude1
        derivative_of_peak['der_peak_amplitude2']=Value_der_peak_amplitude2
        derivative_of_peak['der_peak_sigma1']=Value_der_peak_sigma1
        derivative_of_peak['der_peak_sigma2']=Value_der_peak_sigma2
        derivative_of_peak['der_peak_mean1']=Value_der_peak_mean1
        derivative_of_peak['der_peak_mean2']=Value_der_peak_mean2

        



#        print error_on_peak,"error on peak"

        # now derterminants of the cap sigma for calculating the error on cap sigma X and Y
        # i will read it like that derivative of cap sigma w.r.t parameters



        # derivative of cap sigma w.r.t amplitude1
        Value_der_Capsigma_amplitude1= ( 1.0 / peak - ( (amplitude1_value + amplitude2_value) * Value_der_peak_amplitude1 )  / peak**2 ) / math.sqrt(2*math.pi)




        # derivative of cap sigma  w.r.t amplitude2

        Value_der_Capsigma_amplitude2= ( 1.0 / peak - ( (amplitude1_value + amplitude2_value) * Value_der_peak_amplitude2 )  / peak**2 ) / math.sqrt(2*math.pi)



        
        # derivative of cap sigma  w.r.t sigma1

        Value_der_Capsigma_sigma1 = -( (amplitude1_value + amplitude2_value) * Value_der_peak_sigma1 )/ math.sqrt(2*math.pi) / peak**2




        # derivative of cap sigma w.r.t sigma2

        Value_der_Capsigma_sigma2 = -( (amplitude1_value + amplitude2_value) * Value_der_peak_sigma2 ) / math.sqrt(2*math.pi) / peak**2



       # derivative of cap sigma w.r.t mean1

        Value_der_Capsigma_mean1 = ( (amplitude1_value + amplitude2_value) * Value_der_peak_mean1 ) / math.sqrt(2*math.pi) / peak**2





       # derivative of cap sigma w.r.t mean2
        Value_der_Capsigma_mean2 = ( (amplitude1_value + amplitude2_value) * Value_der_peak_mean2 ) / math.sqrt(2*math.pi) / peak**2





        # off diagnol capsigma terms
        off_diagnol_Capsigma_terms = offdiagnol(Value_der_Capsigma_amplitude1, Value_der_Capsigma_sigma1, Value_der_Capsigma_mean1,Value_der_Capsigma_amplitude2, Value_der_Capsigma_sigma2, Value_der_Capsigma_mean2,covariance_element['cov_amp1_mean1'], covariance_element['cov_amp1_sigma1'],covariance_element['cov_amp1_amp2'],covariance_element['cov_amp1_sigma2'],covariance_element['cov_amp1_mean2'],covariance_element['cov_mean1_sigma1'],covariance_element['cov_mean1_amp2'], covariance_element['cov_mean1_sigma2'],covariance_element['cov_mean1_mean2'],covariance_element['cov_sigma1_amp2'],covariance_element['cov_sigma1_sigma2'], covariance_element['cov_sigma1_mean2'],covariance_element['cov_amp2_sigma2'], covariance_element['cov_amp2_mean2'],covariance_element['cov_sigma2_mean2'] )






        # diagnol capsigma terms



        diagnol_Capsigma_terms = ( amplitude1_error * Value_der_Capsigma_amplitude1 )**2 + ( amplitude2_error * Value_der_Capsigma_amplitude2 )**2 + ( sigma1_error * Value_der_Capsigma_sigma1 )**2 + ( sigma2_error * Value_der_Capsigma_sigma2 )**2 + ( mean1_error * Value_der_Capsigma_mean1 )**2 + ( mean2_error * Value_der_Capsigma_mean2 )**2


        


       # Cap sigma error

        
        Capsigma_error = math.sqrt(off_diagnol_Capsigma_terms + diagnol_Capsigma_terms )


        # its time to store the values of cap sigma derivatives w.r.t parameters in a dictionary

        derivative_of_Capsigma = dict()
        derivative_of_Capsigma['der_Capsigma_amplitude1']=Value_der_Capsigma_amplitude1
        derivative_of_Capsigma['der_Capsigma_amplitude2']=Value_der_Capsigma_amplitude2
        derivative_of_Capsigma['der_Capsigma_sigma1']=Value_der_Capsigma_sigma1
        derivative_of_Capsigma['der_Capsigma_sigma2']=Value_der_Capsigma_sigma2
        derivative_of_Capsigma['der_Capsigma_mean1']=Value_der_Capsigma_mean1
        derivative_of_Capsigma['der_Capsigma_mean2']=Value_der_Capsigma_mean2
        





        return Capsigma, peak, error_on_peak,Capsigma_error, derivative_of_peak, derivative_of_Capsigma


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
    peak_derivative_dict = sigma_fit[4]

    # derivative of cap sigma
    derivative_of_Capsigma = sigma_fit[5]




    # 0 --> sigma for luminosity
    # 1 --> chi/ndf_double gaussian fit
    # 2 -->peak
    # 3 --> chi2 double gauss fit
    # 4 --> status of fit
    # 5 --> status of covariance matrix
    # 6 --> peak error
    # 7 --> error for cap sigma
    # 8 --> covariance element
    # 9 --> dictionary for derivative of peak
    # 10 --> dictionary for derivative of cap sigma


    # returning the values of analyse_scan
    return sigma_for_luminosity, chi_NDF_double_gauss_fit, peak, chi2_double_gauss_fit, status_of_fit, staus_cov_matrix , error_peak, error_Capsigma, covariance_element, peak_derivative_dict, derivative_of_Capsigma







#### empty arrays


chi2_NDF_x=np.array([])
chi2_NDF_y=np.array([])


# peak (relative error)
rel_error_on_peak_x = np.array([])
rel_error_on_peak_y = np.array([])

# rel error on capsigma
rel_error_capsigma_x =np.array([])
rel_error_capsigma_y=np.array([])








for filename in fit_pickle:
    with open(filename,'rb') as fit_pickle:
        scan = pickle.load(fit_pickle)
        #scans = scan['TightModLumi']

        for i_scan in range(len(scan)):
            output = scan[i_scan]
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
            



            if output['bunches'] >2000: # bunches cut


                if result_x[4]==0 and result_y[4]==0 :# converge fit
             
                    
                    if result_x[5] == 3 and result_y[5]==3: ## covariance matrix cut

                        
                        if result_x[1]<7 and result_y[1]<7:#chi/ndf cut



                           

                            #chi2/ndf X scan
                            chi2_NDF_x= np.append(chi2_NDF_x,result_x[1])

#                            print len(chi2_NDF_x),"chi/ndf"

                            histogram(100,0,20,chi2_NDF_x,"Chi_NDF_x","plots/Chi2_NDF_x.png","chi2/NDF_x")


                            #chi2/ndf Y scan
                            chi2_NDF_y=np.append(chi2_NDF_y,result_y[1])
                            histogram(100,0,20,chi2_NDF_y,"Chi_NDF_y","plots/Chi2_NDF_y.png","chi2/NDF_y")




                            # peak x relative error
                            relative_error_peak_x = (result_x[6]/result_x[2])*100


                            # appending peak x array
                            rel_error_on_peak_x = np.append(rel_error_on_peak_x,relative_error_peak_x)
                            # histogram
                            histogram(100,0,0.5,rel_error_on_peak_x,"rel_error_on_peak_x","plots/rel_error_on_peak_x.png","rel_error_on_peak_x[%]")



                            

                            # peak y
                            relative_error_peak_y = (result_y[6]/result_y[2])*100
                            # appending peak x array
                            rel_error_on_peak_y = np.append(rel_error_on_peak_y,relative_error_peak_y)
                            # histogram
                            histogram(100,0,0.5,rel_error_on_peak_y,"rel_error_on_peak_y","plots/rel_error_on_peak_y.png","rel_error_on_peak_y[%]")


                        



                                # sigma error
                            rel_sigma_x = (result_x[7]/result_x[0])*100
                            rel_error_capsigma_x = np.append(rel_error_capsigma_x,rel_sigma_x)
                                # histogram

                            histogram(100,0,3,rel_error_capsigma_x,"rel_error_on_sigma_x","plots/rel_error_on_sigma_x.png","rel_error_on_sigma_x[%]")






                                # y scan sigma
                            rel_sigma_y = (result_y[7]/result_y[0])*100
                            rel_error_capsigma_y = np.append(rel_error_capsigma_y,rel_sigma_y)

                            histogram(100,0,3,rel_error_capsigma_y,"rel_error_on_sigma_y","plots/rel_error_on_sigma_y.png","rel_error_on_sigma_y[%]")






                            
                            # now expected luminosity
                            
                            luminosity_of_beam =expected_luminosity(result_x[0],result_y[0],product)



                            # error on lumi
                            absolute_luminosity = math.sqrt((result_y[7]/result_y[0])**2 + (result_x[7]/result_x[0])**2)



                            # error in luminsoity
                            error_on_luminosity = absolute_luminosity * luminsoity_of_beam


                            
                            # relative luminosity
                            relative_luminosity= (error_on_luminosity/luminosity_of_beam)*100
                            

                            # relative lumisnoity array
                            rel_error_on_lumi = np.append(rel_error_on_lumi,rel_lumi)
 
                            #histogram

                            
                            histogram(100,0,3,rel_error_on_lumi,"rel_error_on_lumi","plots/rel_error_on_lumi.png","rel_error_on_luminosity[%]")



                           # figure of merit

                            Figure_of_Merit= result_x[2] / luminosity_of_beam

                            # appending

                            FOM=np.append(FOM,Figure_of_merit)


                             # histogram
                            histogram(100,0.9,1.05,FOM,"FOM","plots/FOM.png","Figure of merit(peak/lumi)")




                            





























    

