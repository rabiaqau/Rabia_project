#Rabia shaheen
#16_nov_2019

# Parts of script

# 29 Nov 
# here we reLIze that we dont know the 
# define functions which we will use for "main anlysis"
# main fucntion ("analysis scan")
# fits(single and double gauss)
# applying formula for the chi2
# then after single fit i call the beam_beam_deflection function to get value of correction
# redo the fit double gauss
# 

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

 
    # FUNCTIONS section

    # 1. Double gauss fit
def gauss2D(x, p):

    # Parameter 0: amplitude                                                                                                        
    # MAX is: p[0] / math.sqrt(2 * math.pi) / p[2] *( p[3] + (1-p[3])/p4 )  
    # Parameter 1: mean (common between two gaussian)                                                                                
    # Parameter 2: width of first gaussian                                                                                           
    # Parameter 3: Split of amplitude between two gaussians range is 0-1)                                                            
    # Parameter 4: Factor by which second gaussian is wider than first, i.e. sigma2 = sigma1 / parameter 4, p4 is at most 1          

    return p[0] / math.sqrt(2 * math.pi) / p[2] * ( p[3]  * np.exp( - ( (x[0]-p[1])**2 / (2.*p[2]**2) ) ) + ( 1.0 - p[3] ) * p[4]*np.exp( - ( (x[0]-p[1])**2/(2.*(p[2]/p[4])**2) ) ) )      
    
    # loading pickle file    
fit_pickle = sorted(glob.glob("fits/integrated/fit_integrated_trk_run_364098.pickle"))

    # 3.main function (analyze function)
def analyse_scan(scan,separation, luminosity, error, path, run_number, bunches, current, Xscan):

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
    graph.GetXaxis().SetTitle("separation")
    graph.GetYaxis().SetTitle("mu")

    #Single gauss Fit 
    single_gauss_fit = ROOT.TF1( "single_gauss_fit", 'gaus', -3, 3 )# single gauss (3 parameters)
    #Double gauss Fit
    double_gauss_fit = ROOT.TF1("double_gauss_fit", gauss2D, -3, 3, 5)#2d gauss (5 parameters)

    # fitting the graph with single gauss fit
    graph.Fit('single_gauss_fit','Q')

    # 1st method of getting single gauss  parameters
    p0 = single_gauss_fit.GetParameter(0)
    #print p0
    p1 = single_gauss_fit.GetParameter(1)
    # print p1
    p2 = single_gauss_fit.GetParameter(2)
    # print p2

    # SETTING parameters for double gauss
    double_gauss_fit.SetParameter(0,p0)#1st parameter 
    double_gauss_fit.SetParameter(1,p1)#2nd parameter
    double_gauss_fit.SetParameter(2,p2)#3rd parameter
    
    # second method of getting single gauss parameters
    #for i in range(3):            
    #    double_gauss_fit.SetParameter(i, single_gauss_fit.GetParameter(i))#get parameters of double gaussain
        
    # double gauss setting parameters
    double_gauss_fit.SetParameter(3,0.5)#4rd parameter                                                                     
    double_gauss_fit.SetParameter(4,1.3)#5th parameter 
    double_gauss_fit.SetParLimits( 3, 0.0, 1.0 )
    double_gauss_fit.SetParLimits( 4, 1.0, 5.0 )                       
    # fit double gauss                                                        
    graph.Fit('double_gauss_fit','Q')#fitting of double gaussain
    ###########
    print p0,"1st"
    print double_gauss_fit.GetParameter(0),"1"
    print double_gauss_fit.GetParameter(1),"2"
    print double_gauss_fit.GetParameter(2),"3"
    print double_gauss_fit.GetParameter(3),"4th"
    print double_gauss_fit.GetParameter(4),"5th"
    # ****** without corrected separation *****    

    # chi2 value
    # luminsoity from fit (each point)
  
    chi2_uncorrected = np.array([])
    lumi_fit = gauss2D([separation], double_gauss_fit.GetParameters())
    #print lumi_fit,"lumi from fit without"
    #print luminosity, "luminsoiuty from data"
    #print error,"error"

    
    # looping on all separation points
    # to get lumionsity of each separation points
    # lastly calculate the chi2
    subtraction_lumi = np.array([])
    for sep in range(separation_points):  
        subtraction_luminosity = (luminosity[sep]- lumi_fit[sep])
        subtraction_lumi= np.append(subtraction_lumi,subtraction_luminosity)
        
        chi2 = (luminosity[sep]- lumi_fit[sep])**2/error[sep]# chi2 formula
        chi2_uncorrected = np.append(chi2_uncorrected,chi2)
    print subtraction_lumi,"subtraction"
    #print chi2_uncorrected, "chi2_uncorrected"
    # chi2 from fit before correction
    chi2_double_gauss_fit = double_gauss_fit.GetChisquare()
    NDF_double_gauss_fit =  double_gauss_fit.GetNDF()
    print chi2_double_gauss_fit, "fit chi without correction"
    Prob_double_gauss_fit = double_gauss_fit.GetProb()
    print Prob_double_gauss_fit,"probabilty"
    graph.Draw("AP")
    test = ROOT.TGraphErrors(len(separation), separation,subtraction_lumi,np.zeros(separation_points), error)
    
    c = ROOT.TCanvas('fitVsdata', 'fitVsdata')
    c.Divide(1,2) # 1 row, 2 columns
    c.cd(1)
    graph.Draw("AP")
    c.cd(2)
    test.Draw("AP")
    test.GetXaxis().SetTitle("Separation[mm]");
    test.GetYaxis().SetTitle("subtraction_lumi");

    fileName = path + str(run_number) +"_"+str(scan)+"_"
    if Xscan== True:
        fileName += "x_scan_comp.png"
    else:
        fileName += "y_scan_comp.png"
    c.Print(fileName)
    
    #float
    if chi2_double_gauss_fit == False:
        print chi2_double_gauss_fit,"chi from fit"

    # #float
    # if NDF_double_gauss_fit == True:
    #     print NDF_double_gauss_fit,  "NDF From double gauss fit" 

   
    # #ndarray
    if chi2_uncorrected.any()==False:
        print chi2_uncorrected ,'chi values from formula'

    # #ndarray
    if separation.any()==False:
        print separation ,"separation"

    

    # grap for the 
    graph.Draw("AP")
    fileName = path + str(run_number) +"_"+str(scan)+"_"
    if Xscan == True:
        fileName = fileName +"x_scan.png"
    else:
        fileName = fileName+ "y_scan.png"
        
    canvas_name.Print(fileName)

for filename in fit_pickle:
    with open(filename,'rb') as fit_pickle:
        scans = pickle.load(fit_pickle)
        #print "Number of scans", len(scans)
        for i_scan in range(len(scans)):
            output = scans[i_scan]
            print output['chi2'], "output from alex"
            analyse_scan(i_scan,output['sep_x'],output['exp_lumi_x'],output['err_x'],"Emittance_scan/",output['run'],output['bunches'],output['curr'][0],True)
            #analyse_scan(i_scan,output['sep_y'],output['exp_lumi_y'],output['err_y'],"Emittance_scan/",output['run'],output['bunches'],output['curr'][1],False)


