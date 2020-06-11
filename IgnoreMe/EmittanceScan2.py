

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
from EmittanceScan1 import *



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
                        
