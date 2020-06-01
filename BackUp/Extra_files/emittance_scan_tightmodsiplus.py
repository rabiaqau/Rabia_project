


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
from functions import*

fit_pickle_trk = sorted(glob.glob("fits_4WP/integrated/fit_integrated_trk_run_348354.pickle"))
#fit_pickle = sorted(glob.glob("fits/integrated/*trk*"))


#fit_pickle_trk = sorted(glob.glob("fits_4WP/integrated/fit_integrated_trk_run_349327.pickle"))
#fit_pickle = sorted(glob.glob("fits_4WP/integrated/fit_integrated_lcd_run_348251.pickle"))
#fit_pickle_lcd = sorted(glob.glob("fits_4WP/integrated/*lcd*"))
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
    graph.Fit('single_gauss_fit','Q')



# for my fucntions

    amp1 = single_gauss_fit.GetParameter(0)

    error_amp1= double_gauss_fit.GetParError(0)
  #  print (error_amp1/amp1)*100,"consant error"
    m1 = single_gauss_fit.GetParameter(1)

    error_m1= double_gauss_fit.GetParError(1)
   # print (error_m1/m1),"single mean error"
    print bunches,"bunches"

    print m1,"m1"
    s1 = single_gauss_fit.GetParameter(2)
    error_s1= double_gauss_fit.GetParError(2)
#    print (error_s1/s1)*100 ,"sigma error"


    print s1,"s1"

    # getting value for the double gauss fit before we have only constant term

    amplitude = amp1 * s1 * math.sqrt(2*math.pi)
    print amplitude,"amplitude"


    m1_initial= m1
    s1_initial = s1
    s2_initial= 2*s1


    double_gauss_fit.SetParameter(4,s1*2)
    double_gauss_fit.SetParameter(1,m1)
    double_gauss_fit.SetParameter(2,s1)
    double_gauss_fit.SetParameter(5,m1)
    double_gauss_fit.SetParameter(0,amplitude*0.5)
    double_gauss_fit.SetParameter(3,amplitude*0.5)




                   

    status= graph.Fit('double_gauss_fit','')#fitting of double gaussain




    if int(status)== 4:

        ROOT.TVirtualFitter.Fitter(graph).SetMaxIterations(10000)    
        status2 =graph.Fit('double_gauss_fit','')
        if int(status2)==4:
            graph.Fit('double_gauss_fit','')
  
#    if first_chi2 > 100:
#        graph.Fit('double_gauss_fit','')
    
    

    # chi2 from fit before correction
    chi2_double_gauss_fit = double_gauss_fit.GetChisquare()
    print "-----------------------------"
    print chi2_double_gauss_fit,"chi2"

    print "-----------------------------"
    NDF_double_gauss_fit =  double_gauss_fit.GetNDF()
    Prob_double_gauss_fit = double_gauss_fit.GetProb()
    chi_NDF_double_gauss_fit= chi2_double_gauss_fit/NDF_double_gauss_fit#ratio chi/NDF 

    A1= double_gauss_fit.GetParameter(0)
    error_a1=double_gauss_fit.GetParError(0)
    relative_a1 = (error_a1/A1)* 100
        
    A2=double_gauss_fit.GetParameter(3)
    error_a2=double_gauss_fit.GetParError(3)
    relative_a2 = (error_a2/A2)* 100


    M1= double_gauss_fit.GetParameter(1)
    error_m1= double_gauss_fit.GetParError(1)
    relative_m1= (error_m1/M1)*100
  
    M2=double_gauss_fit.GetParameter(5)
    error_m2= double_gauss_fit.GetParError(5)
    relative_m2= (error_m2/M2)*100
    #print relative_m1,"m1 error"
#    print bunches,"bunches"

    S1= double_gauss_fit.GetParameter(2)
    error_s1= double_gauss_fit.GetParError(2)
    relative_s1= (error_s1/S1)*100
    #print relative_s1,"s1"

    S2= double_gauss_fit.GetParameter(4)
    error_s2= double_gauss_fit.GetParError(4)
    relative_s2= (error_s2/S2)*100

    print "--------------------"
    print relative_a1,"a1 error"
    print relative_a2,"a2 error"
    print relative_m1,"m1 error"
    print relative_m2,"m2 error"
    print relative_s1,"s1 error"
    print relative_s2,"s2 error"
    print "--------------------"


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
    
#canvas
    test = ROOT.TGraphErrors(len(separation), separation,subtraction_lumi,np.zeros(separation_points), error)
    
    c = ROOT.TCanvas('sub_lumi', 'sub_lumi')
    c.Divide(1,2) # 1 row, 2 columns
    c.cd(1)
    graph.GetXaxis().SetTitle("separation[mm]")
    graph.GetYaxis().SetTitle("mu")
    graph.GetXaxis().SetTitleSize(0.064)
    graph.GetYaxis().SetTitleSize(0.064)


    graph.Draw("AP")


    c.cd(2)
    test.Draw("AP")

                  
    # path
    fileName = path + str(run_number) +"_"+str(scan)+"_"
    if Xscan== True:
        fileName += "x_scan_comp.png"
    else:
        fileName += "y_scan_comp.png"
    c.Print(fileName)
    
    # sigma
    def sigma_fit():

        peak = double_gauss_fit.GetMaximum()#max

        Xmax=  double_gauss_fit.GetMaximumX()

        # parameters                                                            

        der_a1= (1.0 / (math.sqrt(2 * math.pi) * S1) ) * np.exp( - (1.0 / 2.0) * ( ( Xmax - M1)**2 / (S1 **2) ) )# checked
        
        
        der_a2= (1.0 / ( math.sqrt(2 * math.pi) * S2) ) * np.exp( - (1.0 / 2.0) * ( ( Xmax - M2)**2 / (S2 **2) ) )# checked
        
        der_s1= - A1 * (1.0/ ( math.sqrt(2 * math.pi) * S1 **2) ) * np.exp( - (1.0 / 2.0) * ( ( Xmax - M1)**2 / (S1 **2) ) ) + A1 * (1.0/ ( math.sqrt(2 * math.pi) * S1 **4 ) )* (Xmax - M1)**2 * np.exp( - (1.0 / 2.0) * ( ( Xmax - M1)**2 / (S1 **2) ) )# checked

        der_s2= - A2 * (1.0/ ( math.sqrt(2 * math.pi) * S2 **2 ) ) * np.exp( - (1.0 / 2.0) * ( ( Xmax - M2)**2 / (S2 **2) ) ) + A2 * (1.0/ ( math.sqrt(2 * math.pi) * S2 **4 ) )  * (Xmax - M2)**2 * np.exp( - (1.0 / 2.0) * ( ( Xmax - M2)**2 / (S2 **2) ) )# checked

        der_m1=  A1 * (Xmax - M1) * ( 1.0/ ( math.sqrt(2 * math.pi) * S1 **3) ) * np.exp( - (1.0 / 2.0) * ( ( Xmax - M1) **2 / (S1 **2) ) )
#        print der_m1,"der m1"

        der_m2=  A2 * (Xmax - M2) * ( 1.0/ ( math.sqrt(2 * math.pi) * S2 **3) ) * np.exp( - (1.0 / 2.0) * ( ( Xmax - M2) **2 / (S2 **2) ) )

        uncertanity_peak = math.sqrt( (error_a1 * der_a1)**2 + (error_a2 * der_a2)**2 + (error_s1 * der_s1)**2 + (error_s2 * der_s2)**2 + (error_m1 * der_m1)**2 + (error_m2* der_m2)**2)


#        
        rel= (uncertanity_peak/peak)*100
        print rel,"rel for each scan"

        # if  rel >16:
        #     print "high alert"
        #     print run_number,"run number"
        #     print "bunches"

        #     # print uncertanity_peak ,"peak uncertanity"
        #     print rel,"rel"

            # print relative_a1,"a1 error"

            # print relative_a2,"a2 error"


            # if early:
            #     print "early"

        intlimit = 0.1           

        integral = double_gauss_fit.Integral(-intlimit,intlimit)
    
        error_integral = double_gauss_fit.IntegralError(-intlimit,intlimit)

        
        sigma = (1 / math.sqrt (2 * math.pi)) * integral /peak 


        error_sigma=  math.sqrt( (error_integral/integral)**2 + (uncertanity_peak / peak)**2)
#        print peak,"peak"
#        print integral,"integral"

        
        der_int= (1 / math.sqrt (2 * math.pi))* (1.0 / peak)
        der_peak = (1 / math.sqrt (2 * math.pi))* (integral / peak**2)
        uncertanity_sigma = math.sqrt( ( der_int * error_integral )**2 + (der_peak * uncertanity_peak)**2 )
#        print uncertanity_sigma,"other way"
#        print sigma,"sigma"
        
        abs_error= sigma * error_sigma
#        print (error_sigma) *100,"abs error"


#        print abs_error/sigma,"relative error"
        relative= ( abs_error/sigma) *100

        return sigma, peak,abs_error,uncertanity_peak
    sigma_fit = sigma_fit()
#
    peak_error = sigma_fit[3]
#
    sigma_error = sigma_fit[2]
    
    sigma = sigma_fit[0]
    peak =  sigma_fit[1]
    
    
    # path
    graph.Draw("AP")
    fileName = path + str(run_number) +"_"+str(scan)+"_"
    if Xscan == True:
        fileName = fileName +"x_scan.png"
    else:
        fileName = fileName+ "y_scan.png"
        
    canvas_name.Print(fileName)

    return sigma, chi_NDF_double_gauss_fit, peak, chi2_double_gauss_fit,peak_error,sigma_error



# chi_NDF_early_x
hist_chi_early_x_first_chi = ROOT.TH1F( "chi_early_x_f", "chi_early_x_f", 100, 0, 500 )
hist_chi_early_y_first_chi = ROOT.TH1F( "chi_early_y_f", "chi_early_y_f", 100, 0, 500 )
hist_chi_late_x_first_chi=  ROOT.TH1F( "chi_late_x_f", "chi_late_x_f", 100, 0, 500 )
hist_chi_late_y_first_chi = ROOT.TH1F( "chi_late_y_f", "chi_late_y_f", 100, 0, 500 )


hist_chi_early_x_second_chi = ROOT.TH1F( "chi_early_x_r", "chi_early_x_r", 100, 0, 500 )
hist_chi_early_y_second_chi = ROOT.TH1F( "chi_early_y_r", "chi_early_y_r", 100, 0, 500 )
hist_chi_late_x_second_chi=  ROOT.TH1F( "chi_late_x_r", "chi_late_x_r", 100, 0, 500 )
hist_chi_late_y_second_chi = ROOT.TH1F( "chi_late_y_r", "chi_late_y_r", 100, 0, 500 )


hist_chi_NDF_late_x = ROOT.TH1F( "chi_NDF_late_x", "chi_NDF_late_x", 100, 0, 50 )



hist_chi_NDF_late_y = ROOT.TH1F( "chi_NDF_late_y", "chi_NDF_late_y", 100, 0, 50 )
# histogram for fom _early
hist_fom_early_error = ROOT.TH1F( "fom_error_early", "fom_error_early", 100, 0 , 100 )
hist_fom_late_error = ROOT.TH1F( "fom_error_late", "fom_error_late", 100, 0 , 100 )

hist_fom_late = ROOT.TH1F( "fom_late", "fom_late", 100, 0 , 2 )

# lumisoity for early and late
date_early = np.array([])
date_late = np.array([])
date_fom = np.array([])
# expected lumi
expected_luminosity_early = np.array([])
expected_luminosity_late = np.array([])


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

# np arrays for sigma
sigma_values_x_early = np.array([])
sigma_values_x_late = np.array([])
sigma_errors_x_late = np.array([])
sigma_errors_y_late = np.array([])
sigma_values_y_early = np.array([])
sigma_values_y_late = np.array([])
sigma_errors_x_early= np.array([])

sigma_errors_y_early= np.array([])
error_on_fom = np.array([])

# runs
run_number_early = np.array([])
run_number_late = np.array([])
run_number=  np.array([])


run_error_25= np.array([])
error_25 = np.array([])
fom_error_25 = np.array([])
late_fom_error = np.array([])
run_error_50= np.array([])
error_50 = np.array([])
fom_error_50 = np.array([])
chi_x_early=  np.array([])
chi_y_early=  np.array([])
result = []

if fit_pickle_trk:

    for filename in fit_pickle_trk:
        with open(filename,'rb') as fit_pickle:
            scan = pickle.load(fit_pickle)

            scans = scan['TightModSiPlusLumi']



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
                result_x = analyse_scan(i_scan, sep_x, exp_lumi_x ,err_x ,"Emittance_scan_tightmodsi/", run, bunches, curr_0, True, early,date )
                result_y = analyse_scan(i_scan, sep_y, exp_lumi_y ,err_y ,"Emittance_scan_tightmodsi/", run, bunches, curr_1, False, early,date)

                n1= np.mean(output['curr'][0])/output['bunches']
                n2= np.mean(output['curr'][1])/output['bunches']
            # product of intensities
                product = n1*n2    
                

                if output['early'] == True:
            


                # chi_NDF_early
                    chi_NDF_early_x = np.append(chi_NDF_early_x , result_x[1])
                    chi_NDF_early_y = np.append(chi_NDF_early_y , result_y[1])

                # run number early

    
                # appending with sigma values
                    sigma_values_x_early = np.append(sigma_values_x_early,result_x[0])
                    sigma_errors_x_early= np.append(sigma_errors_x_early,result_x[5])

                    sigma_values_y_early= np.append(sigma_values_y_early,result_y[0])
                    sigma_errors_y_early = np.append(sigma_errors_y_early,result_y[5])

                    
                    peak_errors =np.append(peak_errors, result_x[4])


                    abs_lumi = math.sqrt( (result_x[5]/ result_x[0])**2 + (result_y[5]/ result_y[0])**2)



                # expected luminosity
                    expected_lumi_early= expected_luminosity(result_x[0],result_y[0], product)


                    expected_luminosity_early= np.append(expected_luminosity_early,expected_lumi_early)


                    error_on_lumi = expected_lumi_early * abs_lumi

                    errors_luminosity = np.append(errors_luminosity, error_on_lumi)
                
                    mu_peak_early = result_x[2]
                    Mu_value_early= np.append(Mu_value_early,result_x[2])


                # fom 
                    Figure_of_merit_early = mu_peak_early/expected_lumi_early


                    abs_fom = math.sqrt(( result_x[4]/result_x[2] )**2 + (error_on_lumi / expected_lumi_early)**2)
                
                    error_fom = Figure_of_merit_early * abs_fom 
                    rel_early = (error_fom/Figure_of_merit_early) *100

                    print rel_early,"fom early"
                    hist_fom_early_error.Fill(rel_early)
                    if 8<rel_early<16:

                
                        error_on_fom = np.append(error_on_fom, error_fom)

                        run_number_early = np.append(run_number_early,output['run'])

                        FOM_early = np.append(FOM_early,Figure_of_merit_early)


                if output['early'] == False:
                # chi_NDF_ late

                    chi_NDF_late_x = np.append(chi_NDF_late_x , result_x[1])
                    chi_NDF_late_y = np.append(chi_NDF_late_y , result_y[1])

                #run number

                # expected lumi
                    expected_lumi_late= expected_luminosity(result_x[0],result_y[0], product)

                #  peak from the curve
                    mu_peak_late = result_x[2]

                    Mu_value_late= np.append(Mu_value_late,result_x[2])

                    Figure_of_merit_late = mu_peak_late/expected_lumi_late
                

                    


                    # filling histogram
                    hist_fom_late.Fill(Figure_of_merit_late)
            
                # sigma
                    sigma_values_y_late= np.append(sigma_values_y_late,result_y[0])
                    sigma_errors_x_late= np.append(sigma_errors_x_late,result_x[5])
                    sigma_errors_y_late= np.append(sigma_errors_y_late,result_y[5])

                    sigma_values_x_late = np.append(sigma_values_x_late,result_x[0])

                    abs_lumi_late = math.sqrt( (result_x[5]/ result_x[0])**2 + (result_y[5]/ result_y[0])**2)
                    luminosity_error_late = abs_lumi_late * expected_lumi_late

                    peak_errors =np.append(peak_errors, result_x[4])

                    abs_fom_late = math.sqrt(( result_x[4]/result_x[2] )**2 + (luminosity_error_late / expected_lumi_late)**2)
                    fom_late_error =abs_fom_late * Figure_of_merit_late


                    relative_fom_late =(fom_late_error/ Figure_of_merit_late)*100

                    print relative_fom_late,"relative error late"
                    if 8<relative_fom_late<16:

                        late_fom_error = np.append(late_fom_error,fom_late_error)
                        run_number_late = np.append(run_number_late,output['run'])
                        FOM_late = np.append(FOM_late,Figure_of_merit_late)




                    hist_fom_late_error.Fill(relative_fom_late)

                    date_late= np.append(date_late, output['date'])            
                    
                    


        output = dict()
        
        output['early_run']=run_number_early
        output['early_fom'] = FOM_early
        output['early_fom_errors'] =error_on_fom
        output['early_sigma_x']=sigma_values_x_early
        output['early_sigma_y']=sigma_values_y_early
        output['early_sigma_x_error']=sigma_errors_x_early
        output['early_sigma_y_error']=sigma_errors_y_early 
        output['early_lumi']= expected_luminosity_early
        output['early_peak']=Mu_value_early
        output['early_lumi_error']=errors_luminosity 
        output['early_peak_error']=peak_errors
        result.append(output)
        
#        print result,"result"
pickle.dump(result, open("WP_output/tightmod.pickle",'wb'))












ROOT.TGaxis.SetMaxDigits(6)

no_of_fom_late =len(FOM_late)
plot_graph_error("error/FOM_late_nofixing.png", "FOM_late", no_of_fom_late, run_number_late,FOM_late,late_fom_error,"run","fom_errors_late")

no_of_run_early =len(run_number_early)
plot_graph_error("error/FOM_early_error.png", "peak", no_of_run_early, run_number_early,FOM_early,error_on_fom,"run","fom")




error= np.concatenate([late_fom_error,error_on_fom])
fom=np.concatenate([FOM_early,FOM_late])
run=np.concatenate([run_number_early,run_number_late])

no_of_run_early =len(run)
plot_graph_error("error/FOM.png", "fom", no_of_run_early, run,fom,error,"run","fom")






no_of_sigma_x_early =len(sigma_values_x_early)
plot_graph_error("error/sigma_values_x_early.png", "sigma_x_early", no_of_sigma_x_early, run_number_early,sigma_values_x_early,sigma_errors_x_early,"run","sigma_x")


no_of_sigma_y_early =len(sigma_values_y_early)
plot_graph_error("error/sigma_values_y_early.png", "sigma_y_early", no_of_sigma_y_early, run_number_early,sigma_values_y_early,sigma_errors_y_early,"run","sigma")


no_of_lumi_x_early =len(expected_luminosity_early)
plot_graph_error("error/lumi_valuesearly.png", "lumiearly", no_of_lumi_x_early, run_number_early,expected_luminosity_early,errors_luminosity,"run","lumi_errors")
no_of_peak_errors_early =len(peak_errors)
plot_graph_error("error/peak.png", "peak", no_of_peak_errors_early, run_number_early, Mu_value_early ,peak_errors,"run","mu")



                
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
# canvas_chi_early_x= ROOT.TCanvas("chi_early_x", "chi_early_x", 800, 600)
# hist_chi_early_x.Draw()
# hist_chi_early_x.SetFillColor(30)
# hist_chi_early_x.SetLineColor(kRed)
# hist_chi_early_x.SetFillStyle(3144)
# canvas_chi_early_x.Print("error/chi_early_x.png") 




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

# no_Mu_value_late_high =len(Mu_value_late_high)
# plot_graph("Mu_FOM/Mu_value_late_high.png", "MU_FOM_late_high",no_Mu_value_late_high, Mu_value_late_high ,FOM_late_high,"Mu_high_late","FOM_high_mu_late")


