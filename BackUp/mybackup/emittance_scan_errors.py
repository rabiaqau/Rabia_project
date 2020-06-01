
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



# working point selection fucntion
def workingpoint(function_name, WP,scan):

      sep_x = WP['sep_x']
      exp_lumi_x = WP['exp_lumi_x']
      err_x = WP['err_x']
      curr_0 = WP['curr'][0]
      date = WP['date']
      sep_y = WP['sep_y']
      exp_lumi_y = WP['exp_lumi_y']
      err_y = WP['err_y']
      curr_1 = WP['curr'][1]
      run = WP['run']
      bunches= WP['bunches']
      early = WP['early']

      n1= np.mean(WP['curr'][0])/WP['bunches']
      n2= np.mean(WP['curr'][1])/WP['bunches']
      product = n1*n2

      result_x = function_name(scan, WP, sep_x, exp_lumi_x ,err_x ,"Emittance_scan/", run, bunches, curr_0, True, early,date )
      result_y = function_name(scan, WP, sep_y, exp_lumi_y ,err_y ,"Emittance_scan/", run, bunches, curr_1, False, early,date)

      return result_x, result_y, product,run, early, sep_x, sep_y, exp_lumi_x, exp_lumi_y




# function to redo the fit
def redofit(fit,name):

    A1=fit.GetParameter(0)
    
    A2=fit.GetParameter(3)
        
    S1=fit.GetParameter(2)
        
    S2=fit.GetParameter(4)     
        
    M1=fit.GetParameter(1)     
    M2=fit.GetParameter(5)
    
    fit.SetParameter(0,A1)                                                                                
    fit.SetParameter(3,A2)                                                                                    
    fit.SetParameter(2,S1)                                                                                    
    fit.SetParameter(4,S2)                                                                                    
    fit.SetParameter(1,M1)                                                                                     
    fit.SetParameter(5,M2)

    status_fit2 = name.Fit('double_gauss_fit','Q')
    chi2_fit = fit.GetChisquare()
    status2= int(status_fit2)


    print status2,"status of 2nd fit"
    print chi2_fit,"chi2 after redo of fit"
    return chi2_fit




fit_pickle = sorted(glob.glob("fits_4WP/integrated/fit_integrated_trk_run_350361.pickle"))
#fit_pickle = sorted(glob.glob("fits_4WP/integrated/bad_errors/*trk*"))

#
#fit_pickle = sorted(glob.glob("fits/integrated/fit_integrated_trk_run_348251.pickle#"))
#fit_pickle = sorted(glob.glob("fits_4WP/integrated/fit_integrated_lcd_run_348251.pickle"))
#fit_pickle = sorted(glob.glob("fits_4WP/integrated/*trk*"))



def plot_graph (path_graph, name_graph, number_graph, x, y, graph_title_x,graph_title_y, islowMu = False):
# isdate= False, islowMu = False):

    canvas = ROOT.TCanvas( 'name', 'name',800, 600 )
    graph = ROOT.TGraph(number_graph , x, y)

    graph.SetMarkerColor( 2 )
    graph.SetMarkerStyle( 20 )
    graph.SetMarkerSize( 1.3 )
    graph.GetXaxis().SetTitle(graph_title_x)
    # if isdate:
    #     graph.GetXaxis().SetTimeDisplay(1)
    #     graph.GetXaxis().SetTimeFormat("%d/%m")
    graph.GetYaxis().SetTitle(graph_title_y)
    if islowMu:
        graph.GetXaxis().SetRangeUser(0,5)
    graph.GetXaxis().SetNdivisions(4,4,0)
    graph.Draw("AP")
    canvas.Modified()
    canvas.Print( path_graph )

def gaussD(x, p):                               

    return p[0] * (1.0 / ( math.sqrt(2 * math.pi) * p[2] ) ) * np.exp( - (1.0 / 2.0) * ( (x[0] - p[1])**2 / ( p[2] **2) ) ) + p[3] * (1.0 / (  math.sqrt(2 * math.pi) * p[4] ) ) * np.exp( - (1.0 / 2.0) * ( (x[0] - p[5])**2 / ( p[4] **2) ) )                         

# expected _luminosity         
freq = (1./8.89244e-05) #[s-1]
inel = 8. * 1e-24

def expected_luminosity( sigma_x, sigma_y, current ):
    return freq * current * 1e22 / ( sigma_x * sigma_y  * 2.0 * math.pi) / (freq / inel )

    # 3.main function (analyze function)
def analyse_scan(scan,WP,separation, luminosity, error, path, run_number, bunches, current, Xscan,early,date):

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
    
    # general values for the parameters

    double_gauss_fit.SetParameter(0,0.5*amp1)#amplitude          
    double_gauss_fit.SetParameter(3,0.5*amp1)#amplitude2
    double_gauss_fit.SetParameter(4,2*s1)#sigma2
    double_gauss_fit.SetParameter(1,m1)#mean1 (349592)
    double_gauss_fit.SetParameter(2,s1)#sigma1 
    double_gauss_fit.SetParameter(5,m1)#mean(349592)
    # fix paarmeter
    double_gauss_fit.FixParameter(0,0.5*amp1)#amplitude          

    if early and not Xscan:
        if run_number == 350220:
        
            double_gauss_fit.SetParameter(0,0.6*amp1)#amplitude 
            double_gauss_fit.SetParameter(3,0.4*amp1)#amplitude2


    if early and Xscan:
        if run_number== 349841:

            double_gauss_fit.SetParameter(0,0.6*amp1)#amplitude          
            double_gauss_fit.SetParameter(3,0.4*amp1)#amplitude2
        if run_number == 350013:

            double_gauss_fit.SetParameter(0,0.6*amp1)#amplitude          
            double_gauss_fit.SetParameter(3,0.4*amp1)#amplitude2

    # fit double gauss                                                        
 
    legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
    legend.SetHeader("Double Gauss Fit")
    legend.AddEntry(graph,"data","P")
    legend.AddEntry(double_gauss_fit,"fit","L")
 
    ROOT.TVirtualFitter.Fitter(graph).SetMaxIterations(10000)                   
    status= graph.Fit('double_gauss_fit','Q')#fitting of double gaussain


    if int(status) == 4:
        print "first fit has failed"
        print run_number
        
        print int(status),"status of the first fit"
    
    chi2_double_gauss_fit = double_gauss_fit.GetChisquare()
#    print chi2_double_gauss_fit,"chi of the first fit"

    if scans_TightModLumi:

        if early:

            if run_number== 350361:
                redofit(double_gauss_fit, graph)

            if run_number== 350440:
                redofit(double_gauss_fit, graph)

        if not early:

            if run_number== 359735:
                redofit(double_gauss_fit, graph)

            if run_number== 357713:
                redofit(double_gauss_fit, graph)                

            if run_number== 357620:
                redofit(double_gauss_fit, graph)

            if run_number== 351671:
                redofit(double_gauss_fit, graph)

            if run_number== 351223:
                redofit(double_gauss_fit, graph)

            if run_number== 355331:
                redofit(double_gauss_fit, graph)

        if not Xscan and not early:
            if run_number== 359124:
                redofit(double_gauss_fit, graph)

            if run_number== 350676:
                redofit(double_gauss_fit, graph)

        if early and not Xscan:

            if run_number == 358031:
                redofit(double_gauss_fit, graph)

            if run_number == 355468:
                redofit(double_gauss_fit, graph)

            if run_number == 355416:
                redofit(double_gauss_fit, graph)

            if run_number == 350682:
                redofit(double_gauss_fit, graph)

            if run_number == 362661:
                redofit(double_gauss_fit, graph)

            if run_number == 363710:
                redofit(double_gauss_fit, graph)
            
        if not Xscan and not early:

            if run_number == 363129:
                redofit(double_gauss_fit, graph)
            

        if Xscan and not early:

            if run_number == 359124:
                redofit(double_gauss_fit, graph)
            
        
            if run_number == 350310:
                redofit(double_gauss_fit, graph)
            

        if Xscan and not early:
            if run_number == 350440:
                redofit(double_gauss_fit, graph)

            if run_number == 350160:
                redofit(double_gauss_fit, graph)

        if not early:
            if run_number == 349944:
                redofit(double_gauss_fit, graph)

        if early and not Xscan:

            if run_number == 350067:
                redofit(double_gauss_fit, graph)
            
            if run_number == 349841:
                redofit(double_gauss_fit, graph)

            if run_number == 349592:
                redofit(double_gauss_fit, graph)

            if run_number == 349646:
                redofit(double_gauss_fit, graph)
            
            if run_number == 349693:
                redofit(double_gauss_fit, graph)

            if run_number == 350220:
                redofit(double_gauss_fit, graph)


    # chi2 from fit before correction
    chi2_double_gauss_fit = double_gauss_fit.GetChisquare()
 #   print chi2_double_gauss_fit,"at the end"
#    print chi2_double_gauss_fit,"chi2"
    NDF_double_gauss_fit =  double_gauss_fit.GetNDF()
 #   print NDF_double_gauss_fit,"NDF"

    Prob_double_gauss_fit = double_gauss_fit.GetProb()
    chi_NDF_double_gauss_fit= chi2_double_gauss_fit/NDF_double_gauss_fit#ratio chi/NDF 
    round_chi_double_gauss_fit = round(chi2_double_gauss_fit,2)#round off the value    
    round_chi_NDF_double_gauss_fit = round(chi_NDF_double_gauss_fit,2)# 

 
    graph.SetLineColor( 38 )
    graph.SetMarkerColor( 4 )
    graph.SetMarkerStyle( 20 )
    graph.SetMarkerSize( 1.1 )
  
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

    legend = ROOT.TLegend(0.75, 0.6, 0.9, 0.86)
   # legend.SetHeader("Double Gauss Fit")# for lagend header                    
    legend.AddEntry(graph,"data","P")# points                                      
    legend.AddEntry(double_gauss_fit," Fit","L")# 
    legend.AddEntry(0,'#chi^{2}='+str(round_chi_double_gauss_fit), '')# chi2
    legend.AddEntry(0,'#chi^{2}/ndf='+str(round_chi_NDF_double_gauss_fit), '')# chi
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
    #test.SetTitleY(0.1)
    # line                                                                                                                                     
    line = TLine(-0.037,0,0.037,0.0)
    line.SetLineColor(2)
    line.SetLineWidth(1)
    line.Draw()

    # path


    if scans_TightModLumi:
        fileName = path +"tightmodlumi" + str(run_number) +"_"+str(scan)+"_"
        if Xscan== True:
            fileName += "x_scan_comp.png"
        else:
            fileName += "y_scan_comp.png"
        c.Print(fileName)




    if scans_TightModSiPlusLumi:
        fileName = path +"tightmodsipluslumi" + str(run_number) +"_"+str(scan)+"_"
        if Xscan== True:
            fileName += "x_scan_comp.png"
        else:
            fileName += "y_scan_comp.png"
        c.Print(fileName)






    
    # sigma
    def sigma_fit():

        peak = double_gauss_fit.GetMaximum()#max

        Xmax=  double_gauss_fit.GetMaximumX()
#        print Xmax,"xmax"
        # parameters
        A1= double_gauss_fit.GetParameter(0)
 #       print A1,"A1"
        A2=double_gauss_fit.GetParameter(3)
  #      print A2,"A2"
        S1= double_gauss_fit.GetParameter(2)
  #      print S1,"s1"
        S2= double_gauss_fit.GetParameter(4)
  #      print S2,"s2"
        M1= double_gauss_fit.GetParameter(1)
  #      print M1,"m1"
        M2=double_gauss_fit.GetParameter(5)
  #      print M2,"M2"
        # errors
        error_a1=double_gauss_fit.GetParError(0)
  #      print error_a1,"error a1"
        error_a2=double_gauss_fit.GetParError(3)
  #      print error_a2,"error-a2"
        error_s1=double_gauss_fit.GetParError(2)
  #      print error_s1,"error s1"
        error_s2=double_gauss_fit.GetParError(4)
  #      print error_s2,"error s2"
        error_m1= double_gauss_fit.GetParError(1)
   #     print error_m1,"m1"
        error_m2= double_gauss_fit.GetParError(5)
    #    print error_m2,"m2"


        der_a1= (1.0 / (math.sqrt(2 * math.pi) * S1) ) * np.exp( - (1.0 / 2.0) * ( ( Xmax - M1)**2 / (S1 **2) ) )


        der_a2= (1.0 / ( math.sqrt(2 * math.pi) * S2) ) * np.exp( - (1.0 / 2.0) * ( ( Xmax - M2)**2 / (S2 **2) ) )


        der_s1= - A1 * (1.0/ ( math.sqrt(2 * math.pi) * S1 **2) ) * np.exp( - (1.0 / 2.0) * ( ( Xmax - M1)**2 / (S1 **2) ) ) + A1 * (1.0/ ( math.sqrt(2 * math.pi) * S1 **4 ) )* (Xmax - M1)**2 * np.exp( - (1.0 / 2.0) * ( ( Xmax - M1)**2 / (S1 **2) ) ) 





        der_s2= - A2 * (1.0/ ( math.sqrt(2 * math.pi) * S2 **2 ) ) * np.exp( - (1.0 / 2.0) * ( ( Xmax - M2)**2 / (S2 **2) ) ) + A2 * (1.0/ ( math.sqrt(2 * math.pi) * S2 **4 ) )  * (Xmax - M2)**2 * np.exp( - (1.0 / 2.0) * ( ( Xmax - M2)**2 / (S2 **2) ) )




        der_m1=  A1 * (Xmax - M1) * ( 1.0/ ( math.sqrt(2 * math.pi) * S1 **3) )*np.exp( - (1.0 / 2.0) * ( ( Xmax - M1) **2 / (S1 **2) ) )



        der_m2=  A2 * (Xmax - M2) * ( 1.0/ ( math.sqrt(2 * math.pi) * S2 **3) )  *np.exp( - (1.0 / 2.0) * ( ( Xmax - M2) **2 / (S2 **2) ) )

        uncertanity_peak = math.sqrt( (error_a1 * der_a1)**2 + (error_a2 * der_a2)**2 + (error_s1 * der_s1)**2 + (error_s2 * der_s2)**2 + (error_m1 * der_m1)**2 + (error_m2* der_m2)**2)

        
        intlimit = 0.1           

        integral = double_gauss_fit.Integral(-intlimit,intlimit)
    
        error_integral = double_gauss_fit.IntegralError(-intlimit,intlimit)
        
        sigma = (1 / math.sqrt (2 * math.pi)) * integral /peak 


        der_peak= (1 / math.sqrt (2 * math.pi)) * (integral / peak**2)


        der_int = (1 / math.sqrt (2 * math.pi))* (1.0 / peak)

   #     print sigma
        error_sigma=  math.sqrt( (error_integral/integral)**2 + (uncertanity_peak / peak)**2)



        uncertanity_sigma= math.sqrt( (der_peak * uncertanity_peak)**2 + (der_int * error_integral)**2)

        

        return sigma, peak, uncertanity_peak, uncertanity_sigma

    sigma_fit = sigma_fit()
    
    peak_error = sigma_fit[2]
    sigma_error = sigma_fit[3]
    
    sigma = sigma_fit[0]
    peak =  sigma_fit[1]
    
    
    # path
    graph.Draw("AP")



    if scans_TightModLumi:
        fileName = path+"tightmodlumi"+"_" + str(run_number) +"_"+str(scan)+"_"
        if Xscan == True:
            fileName = fileName +"x_scan.png"
        else:
            fileName = fileName+ "y_scan.png"
        
            canvas_name.Print(fileName)


    if scans_TightModSiPlusLumi:
        fileName = path+"tightmodsipluslumi"+"_" + str(run_number) +"_"+str(scan)+"_"
        if Xscan == True:
            fileName = fileName +"x_scan.png"
        else:
            fileName = fileName+ "y_scan.png"
        
            canvas_name.Print(fileName)



    return sigma, chi_NDF_double_gauss_fit, peak, chi2_double_gauss_fit, peak_error, sigma_error



if fit_pickle:
      for filename in fit_pickle:
            with open(filename,'rb') as fit_pickle:
                  scan = pickle.load(fit_pickle)

                  for i in range(len(scan)):
                        scans_TightModLumi = scan[i]['TightModLumi']
                        scans_TightModSiPlusLumi = scan[i]['TightModSiPlusLumi']


#        scans = scan['Tight']
#        scans = scan['TightLumi']


                  for i_scan in range(len(scans_TightModLumi)):
                        tight_mod =  workingpoint(analyse_scan, scans_TightModLumi[i_scan], i_scan)

                        if tight_mod[4]== True:

                              early_sigma_x= tight_mod[0][0]
                              print early_sigma_x,"early sigma"
                              early_sigma_x_error= tight_mod[0][5]

                              print early_sigma_x_error,"tight mod x "

                              early_sigma_y= tight_mod[1][0]
                              print early_sigma_y,"early y sigma"
                              early_sigma_y_error= tight_mod[1][5]
                              print early_sigma_y_error," y sigma "

                              early_peak =  tight_mod[0][2]
                              print early_peak,"early peak"
                              early_peak_error =  tight_mod[0][4]
                              print early_peak_error,"early peak errro"
                              early_run =  tight_mod[3]
                              early_intensity= tight_mod[2]
                              print early_intensity,"intensity"
                              #expected lumi
                              early_expected_luminosity= expected_luminosity(early_sigma_x,early_sigma_y,early_intensity)
                              print early_expected_luminosity,"expected lumi"
                              # early fom
                              early_fom=  early_peak/ early_expected_luminosity
                              print early_fom,"early fom"
                              early_lumi_quadrature = math.sqrt( (early_sigma_x_error/early_sigma_x)**2 + (early_sigma_y_error/ early_sigma_y)**2)


                              early_error_on_luminosity = early_expected_luminosity * early_lumi_quadrature

                              print early_error_on_luminosity,"error on luminoisty"
                              fom_quadrature = math.sqrt((early_error_on_luminosity/ early_expected_luminosity)**2 + (early_peak_error /early_peak )**2)

                              
                              early_error_fom = early_fom * fom_quadrature
                              print early_error_fom




                        if tight_mod[4]== False:
                              # sigma x
                              late_sigma_x= tight_mod[0][0]
                              late_sigma_x_error= tight_mod[0][5]
                              # sigma y
                              late_sigma_y= tight_mod[1][0]
                              late_sigma_y_error= tight_mod[1][5]
                              # peak
                              late_peak =  tight_mod[0][2]
                              late_peak_error =  tight_mod[0][4]
                              # intensity or run
                              late_intensity= tight_mod[2]
                              late_run = tight_mod[3]
                              # expected lumi
                              late_expected_luminosity= expected_luminosity(late_sigma_x,late_sigma_y,late_intensity)

                              late_lumi_quadrature = math.sqrt( (late_sigma_x_error/late_sigma_x)**2 + (late_sigma_y_error/late_sigma_x)**2)

                              late_error_on_luminosity = late_expected_luminosity * late_lumi_quadrature
                              print late_error_on_luminosity,"late error on luminosity"

                              fom_quadrature = math.sqrt((late_error_on_luminosity/late_expected_luminosity)**2 + (late_peak_error /late_peak )**2)

                              # late fom                                                                                        
                              late_fom=  late_peak/ late_expected_luminosity
                              
                              late_error_fom = fom_quadrature * late_fom
                              print late_error_fom,"late error fom"



####### here i will go for the tight mod lumi si plus

#                   for i_scan_TightModSiPlusLumi in range(len(scans_TightModSiPlusLumi)):
#                         tight_mod_si_plus_lumi =  workingpoint(analyse_scan, scans_TightModSiPlusLumi[i_scan_TightModSiPlusLumi], i_scan_TightModSiPlusLumi)

#                         if tight_mod_si_plus_lumi[4]== True:

#                               early_sigma_x= tight_mod_si_plus_lumi[0][0]
#                               early_sigma_x_error= tight_mod_si_plus_lumi[0][5]
                              
#                               early_sigma_y= tight_mod_si_plus_lumi[1][0]
#                               early_sigma_y_error= tight_mod_si_plus_lumi[1][5]

#                               early_peak =  tight_mod_si_plus_lumi[0][2]
#                               early_peak_error =  tight_mod_si_plus_lumi[0][4]

#                               early_run =  tight_mod_si_plus_lumi[3]
#                               early_intensity= tight_mod_si_plus_lumi[2]
#                               #expected lumi
#                               early_expected_luminosity= expected_luminosity(early_sigma_x,early_sigma_y,early_intensity)
                              
#                               # early fom
#                               early_fom=  early_peak/ early_expected_luminosity

#                               early_lumi_quadrature = math.sqrt( (early_sigma_x_error/early_sigma_x)**2 + (early_sigma_y_error/ early_sigma_x)**2)


#                               early_error_on_luminosity = early_expected_luminosity * early_lumi_quadrature

# #                              print early_error_on_luminosity,"error on luminoisty tight mod si plus"
#                               fom_quadrature = math.sqrt((early_error_on_luminosity/ early_expected_luminosity)**2 + (early_peak_error /early_peak )**2)

                              
#                               early_error_fom = early_fom * fom_quadrature
#                              print early_error_fom




#                         if tight_mod_si_plus_lumi[4]== False:
#                               # sigma x
#                               late_sigma_x= tight_mod_si_plus_lumi[0][0]
#                               late_sigma_x_error= tight_mod_si_plus_lumi[0][5]
#                               # sigma y
#                               late_sigma_y= tight_mod_si_plus_lumi[1][0]
#                               late_sigma_y_error= tight_mod_si_plus_lumi[1][5]
#                               # peak
#                               late_peak =  tight_mod_si_plus_lumi[0][2]
#                               late_peak_error =  tight_mod_si_plus_lumi[0][4]
#                               # intensity or run
#                               late_intensity= tight_mod_si_plus_lumi[2]
#                               late_run = tight_mod_si_plus_lumi[3]
#                               # expected lumi
#                               late_expected_luminosity= expected_luminosity(late_sigma_x,late_sigma_y,late_intensity)

#                               late_lumi_quadrature = math.sqrt( (late_sigma_x_error/late_sigma_x)**2 + (late_sigma_y_error/late_sigma_x)**2)

#                               late_error_on_luminosity = late_expected_luminosity * late_lumi_quadrature
# #                              print late_error_on_luminosity,"late error on luminosity tight mod si plus"

#                               fom_quadrature = math.sqrt((late_error_on_luminosity/late_expected_luminosity)**2 + (late_peak_error /late_peak )**2)

#                               # late fom                                                                                        
#                               late_fom=  late_peak/ late_expected_luminosity
                              
#                               late_error_fom = fom_quadrature * late_fom
# #                              print late_error_fom,"late error fom tight mod si plus"

            
            # for i_scan in range(len(scans_TightModLumi)):                  
            #       output = scans_TightModLumi[i_scan]

            #       sep_x = output['sep_x']
            #       exp_lumi_x = output['exp_lumi_x']
            #       err_x = output['err_x']
            #       curr_0 = output['curr'][0]
            #       date = output['date']
            
            #       sep_y = output['sep_y']
            #       exp_lumi_y = output['exp_lumi_y']
            #       err_y = output['err_y']
            #       curr_1 = output['curr'][1]

            #       run = output['run']
            #       bunches= output['bunches']
            #       early = output['early']

            # #print fom_alex,"fom Alex"
            #       result_x = analyse_scan(i_scan, sep_x, exp_lumi_x ,err_x ,"Emittance_scan/", run, bunches, curr_0, True, early,date )
            #       result_y = analyse_scan(i_scan, sep_y, exp_lumi_y ,err_y ,"Emittance_scan/", run, bunches, curr_1, False, early,date)



            #       n1= np.mean(output['curr'][0])/output['bunches']
            #       n2= np.mean(output['curr'][1])/output['bunches']
            
            # # product of intensities
            #       product = n1*n2
            #       intensity = np.append(intensity,product)
    
# #            print date,"date after loop"
#             expected_lumi_early= expected_luminosity(result_x[0],result_y[0], product)
#             expected_lumi_late= expected_luminosity(result_x[0],result_y[0], product)
#             mu_peak_early = result_x[2]
#             mu_peak_late = result_x[2]
#             mu_peak =mu_peak_early + mu_peak_late
#             expected_lumi= expected_lumi_late + expected_lumi_early
#             FOM_all = mu_peak / expected_lumi
#             FOM = np.append(FOM,FOM_all)
#  #           print FOM,"array"
#             run_number= np.append(run_number, output['run'])
#   #          print run_number,"run"
#             date_fom= np.append(date_fom, date)
# #            print date_fom,"date"
            

# #            hist_status_x.Fill(result_x[6])
# #            hist_status_y.Fill(result_y[6])






#             if output['early'] == True:
            

#                 # print result_x[4],"4x peak"
#                 # print result_y[4],"4y peak"
#                 # print result_x[5],"5x sigma"
#                 # print result_y[5],"5y sigma"

#                 # chi_NDF_early
#                 chi_NDF_early_x = np.append(chi_NDF_early_x , result_x[1])
#                 chi_NDF_early_y = np.append(chi_NDF_early_y , result_y[1])
#                 print chi_NDF_early_x,"chi x"

  #              print chi_NDF_early_y,"chi y"
#                 # run number early
#                 run_number_early = np.append(run_number_early,output['run'])
    
#                 # appending with sigma values
#                 sigma_values_x_early = np.append(sigma_values_x_early,result_x[0])
#                 sigma_errors_x_early= np.append(sigma_errors_x_early,result_x[5])

#   #              print sigma_errors_x_early,"sigma_errors_x"

#                 hist_sigma_xerror_early.Fill(result_x[5])

#                 sigma_values_y_early= np.append(sigma_values_y_early,result_y[0])
#                 sigma_errors_y_early = np.append(sigma_errors_y_early,result_y[5])


# #                print sigma_errors_x_early,"sigma_errors_x"
# #                print sigma_errors_y_early,"sigma_errors_y"
                

#                 peak_errors_early =np.append(peak_errors_early, result_x[4])
    
# #                print peak_errors_early,"early peak errors"



#                 hist_peakerror_early.Fill(result_x[4])

#                 sigma_x= result_x[0]
#                 sigma_error_x= result_x[5]

#                 sigma_y= result_y[0]
#                 sigma_error_y= result_y[5]



#                 abs_lumi = math.sqrt( (sigma_error_x/ sigma_x)**2 + (sigma_error_y/ sigma_y)**2)
# #                print sigma_error_x,"error x"
# #                print sigma_error_y,"error y"


#                 # expected luminosity
#                 expected_lumi_early= expected_luminosity(result_x[0],result_y[0], product)


#                 expected_luminosity_early= np.append(expected_luminosity_early,expected_lumi_early)
# #                print expected_luminosity_early,"expected"
#                 error_on_lumi = expected_lumi_early * abs_lumi

#                 errors_luminosity = np.append(errors_luminosity, error_on_lumi)
# #                print errors_luminosity,"error on lumi"

#                 # peak from the curve

#                 mu_peak_early = result_x[2]
#                 Mu_value_early= np.append(Mu_value_early,result_x[2])

# #                print Mu_value_early,"mu value"
#                 # fom 
#                 Figure_of_merit_early = mu_peak_early/expected_lumi_early
                

#                 abs_fom = math.sqrt(( result_x[4]/result_x[2] )**2 + (error_on_lumi / expected_lumi_early)**2)

                
            


                
#                 error_fom = Figure_of_merit_early * abs_fom 
                
#                 error_on_fom = np.append(error_on_fom, error_fom)
#  #               print error_on_fom,"error on fom early"


# #                print error_on_fom,"fom error"
# #                print error_on_fom,"errror on fom"
 
#                 # if error_fom >0.025:
#                 #     print output['run'] ,"run number early"
#                 #     print error_fom,"error on fom early"
#                 #     print Figure_of_merit_early ,"fom early"
#                 FOM_early = np.append(FOM_early,Figure_of_merit_early)
# #                print FOM_early,"fom early"
                                  
#                 # filing histogram fom

#                 hist_fom_early.Fill(Figure_of_merit_early)
                
#                 # chi_early
#                 chi_early_x = np.append(chi_early_x, result_x[3])


#             if output['early'] == False:
#                 # chi_NDF_ late

#                 # print "late"
#                 # print 
#                 # print result_y[4],"4y peak"
#                 # print result_x[5],"5x sigma"
#                 # print result_y[5],"5y sigma"

#                 chi_NDF_late_x = np.append(chi_NDF_late_x , result_x[1])
#                 chi_NDF_late_y = np.append(chi_NDF_late_y , result_y[1])

#                 #run number
#                 run_number_late = np.append(run_number_late,output['run'])
#                 # expected lumi

#                 peak_errors_late =np.append(peak_errors_late,result_x[4])
#  #               print result_x[4],"late peak"
#                 sigma_values_x_late = np.append(sigma_values_x_late,result_x[0])
#                 sigma_x_late= result_x[0]
#  #               print result_x[5],"late sigma x"
#  #               print result_y[5],"late sigma y"
                
#                 sigma_values_y_late = np.append(sigma_values_y_late,result_y[0])
#                 sigma_y_late = result_y[0]
#                 sigma_errors_x_late= np.append(sigma_errors_x_late,result_x[5])
#                 sigma_x_late_error= result_x[5]
#                 sigma_errors_y_late= np.append(sigma_errors_y_late,result_y[5])
#                 sigma_y_late_error= result_y[5]

#                 expected_lumi_late= expected_luminosity(result_x[0],result_y[0], product)


#                 #  peak from the curve
#                 mu_peak_late = result_x[2]

#                 Mu_value_late= np.append(Mu_value_late,result_x[2])

#                 sigma_errors_x_early= np.append(sigma_errors_x_early,result_x[5])


                
#                 # fom 

#                 abs_lumi_late = math.sqrt( (sigma_x_late_error/ sigma_x_late)**2 + (sigma_y_late_error/ sigma_y_late)**2)

#                 error_lumi_late =abs_lumi_late * expected_lumi_late
#   #              print error_lumi_late,"error lumi late"

#                 abs_fom_late = math.sqrt(( result_x[4]/result_x[2] )**2 + (error_lumi_late / expected_lumi_late)**2)
                

#                 Figure_of_merit_late = mu_peak_late/expected_lumi_late
#                 fom_error_late = abs_fom_late * Figure_of_merit_late

#                 error_on_fom_late = np.append(error_on_fom_late,fom_error_late)
# #                print error_on_fom_late,"fom error late"
#                 FOM_late = np.append(FOM_late,Figure_of_merit_late)
# #                print Figure_of_merit_late ,"fom late"

#                 # if fom_error_late > 0.025:
#                 #     print output['run'] ,"run number late"
#                 #     print fom_error_late,"error late"
#                 #     print Figure_of_merit_late ,"fom late"






#                 # filling histogram
#                 hist_fom_late.Fill(Figure_of_merit_late)
            
#                 # sigma
#                 sigma_values_y_late= np.append(sigma_values_y_late,result_y[0])
#                 sigma_values_x_late = np.append(sigma_values_x_late,result_x[0])
#                 date_late= np.append(date_late, output['date'])            
            
#             if output['early']==True: 

#                 hist_chi_NDF_early_x.Fill(result_x[1])
#                 hist_chi_early_x.Fill(result_x[3])
                
#                 hist_chi_NDF_early_y.Fill(result_y[1])
           

#             if output['early']==False:
                
#                 hist_chi_NDF_late_x.Fill(result_x[1])
#                 hist_chi_NDF_late_y.Fill(result_y[1])




                
# ROOT.TGaxis.SetMaxDigits(6)
# #1 root canvas for Chi_NDF_early_x

# canvas_fom_early= ROOT.TCanvas("fom_early", "fom_early", 800, 600)
# hist_fom_early.Draw()
# hist_fom_early.SetFillColor(30)
# hist_fom_early.SetLineColor(kRed)
# hist_fom_early.SetFillStyle(3144)
# canvas_fom_early.Print("fom/hist_fom_early.png")


# canvas_status_x= ROOT.TCanvas("status_x", "status_x", 800, 600)
# hist_status_x.Draw()
# hist_status_x.SetFillColor(30)
# hist_status_x.SetLineColor(kRed)
# hist_status_x.SetFillStyle(3144)
# canvas_status_x.Print("status/hist_status_x.png")

# canvas_status_y= ROOT.TCanvas("status_y", "status_y", 800, 600)
# hist_status_y.Draw()
# hist_status_y.SetFillColor(30)
# hist_status_y.SetLineColor(kRed)
# hist_status_y.SetFillStyle(3144)
# canvas_status_y.Print("status/hist_status_y.png")







# print hist_fom_early.GetStdDev(),"std of his"



# # fom late
# canvas_fom_late= ROOT.TCanvas("fom_late", "fom_late", 800, 600)
# hist_fom_late.Draw()
# hist_fom_late.SetFillColor(30)
# hist_fom_late.SetLineColor(kRed)
# hist_fom_late.SetFillStyle(3144)
# canvas_fom_late.Print("fom/hist_fom_late.png")

# ####




# canvas_peakerror_early= ROOT.TCanvas("peakerror_early", "peakerror_early", 800, 600)
# hist_peakerror_early.Draw()
# hist_peakerror_early.SetFillColor(30)
# hist_peakerror_early.SetLineColor(kRed)
# hist_peakerror_early.SetFillStyle(3144)
# canvas_peakerror_early.Print("lumi/hist_peakerror_early.png")


# canvas_sigma_xerror_early= ROOT.TCanvas("sigma_xerror_early", "sigma_xerror_early", 800, 600)
# hist_sigma_xerror_early.Draw()
# hist_sigma_xerror_early.SetFillColor(30)
# hist_sigma_xerror_early.SetLineColor(kRed)
# hist_sigma_xerror_early.SetFillStyle(3144)
# canvas_sigma_xerror_early.Print("lumi/hist_sigma_xerror_early.png")




# #
# canvas_chi_NDF_early_x= ROOT.TCanvas("chi_NDF_early_x", "chi_NDF_early_x", 800, 600)
# hist_chi_NDF_early_x.Draw()
# hist_chi_NDF_early_x.SetFillColor(30)
# hist_chi_NDF_early_x.SetLineColor(kRed)
# hist_chi_NDF_early_x.SetFillStyle(3144)
# canvas_chi_NDF_early_x.Print("chi_NDF/chi_NDF_early_x.png")

# #
# canvas_chi_early_x= ROOT.TCanvas("chi_early_x", "chi_early_x", 800, 600)
# hist_chi_early_x.Draw()
# hist_chi_early_x.SetFillColor(30)
# hist_chi_early_x.SetLineColor(kRed)
# hist_chi_early_x.SetFillStyle(3144)
# canvas_chi_early_x.Print("chi_NDF/chi_early_x.png") 




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

# ROOT.TGaxis.SetMaxDigits(6)


# # plot fom_early vs date_early




# no_of_date =len(date_fom)
# plot_graph("fom_date/fom_date.png", "fom_date",no_of_date, date_fom ,FOM,"date","fom",True)


# no_of_date_early =len(date_early)
# plot_graph("fom_date/fom_early_date.png", "fom_early_date",no_of_date_early, date_early ,FOM_early,"date_early","fom_early")


# # date _late
# no_of_date_late =len(date_late)
# plot_graph("fom_date/fom_date_late.png", "fom_late_date" , no_of_date_late, date_late ,FOM_late,"date_late","fom_late")


# #plot for early_x


# no_of_chi_NDF_early_x =len(chi_NDF_early_x)
# plot_graph("chi_NDF_run/chi_NDF_early_x.png", "chi_NDF_early_x" , no_of_chi_NDF_early_x, run_number_early , chi_NDF_early_x,"run","chi_NDF_early_x")
# #plot for early_y

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

# no_of_sigma_x_early =len(sigma_values_x_early)
# plot_graph_error("lumi/sigma_values_x_early.png", "sigma_x_early", no_of_sigma_x_early, run_number_early,sigma_values_x_early,sigma_errors_x_early,"run","sigma_x_early")


# no_of_sigma_y_early =len(sigma_values_y_early)
# plot_graph_error("lumi/sigma_values_y_early.png", "sigma_y_early", no_of_sigma_y_early, run_number_early,sigma_values_y_early,sigma_errors_y_early,"run","sigma_y_early")




# no_of_sigma_x_late =len(sigma_values_x_late)
# plot_graph_error("lumi/sigma_values_x_late.png", "sigma_x_late", no_of_sigma_x_late, run_number_early,sigma_values_x_late,sigma_errors_x_late,"run","sigma_x_late")




# no_of_lumi_x_early =len(expected_luminosity_early)
# plot_graph_error("lumi/lumi_values_x_early.png", "lumi_x_early", no_of_lumi_x_early, run_number_early,expected_luminosity_early,errors_luminosity,"run","lumi_errors")


# no_of_peak_errors_early =len(peak_errors_early)
# plot_graph_error("lumi/peakx_early.png", "peak_x_early", no_of_peak_errors_early, run_number_early,Mu_value_early,peak_errors_early,"run","peak_errors")



# no_of_fom_early =len(FOM_early)
# plot_graph_error("lumi/FOM_early.png", "FOM_early", no_of_fom_early, run_number_early,FOM_early,error_on_fom,"run","fom_errors")



# no_of_fom_late =len(FOM_late)
# plot_graph_error("lumi/FOM_late.png", "FOM_late", no_of_fom_late, run_number_late,FOM_late,error_on_fom_late,"run_late","fom_errors_late")


# no_Mu_value_early =len(Mu_value_early)
# plot_graph_erro("lumi/Mu_FOM_early_error.png", "MU_FOM_early_error",no_Mu_value_early, Mu_value_early ,FOM_early,peak_errors_early,error_on_fom,"MU_early","FOM_early")









# # #2
# # # run vs sigma_late_x 

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


