#Rabia shaheen
#16_nov_2019

# Parts of script

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


    #print corrected_separation
    # pickle files   
    # fit_integrated_trk_run_364030.pickle (done) good
    # fit_integrated_trk_run_364098.pickle (done) good
    # fit_integrated_trk_run_363710.pickle (done) not very good
    # fit_integrated_trk_run_364292.pickle (done) not very good
    # fit_integrated_trk_run_363830.pickle (done) not very good
 
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
    
    # 2.gauss fucntion( for beam_beam_deflection)
 
def gauss(x, *p):    
    return p[0] * np.exp( - ( (x-p[1])**2 / (2.*p[2]**2) ) )

    # 3. Beam_beam_deflection function

    #constants
gamma = 13e6 / 2. / 938
radius = 1.5347e-15

def beam_beam_deflection( separation, current, bunches, mean, width, beta, Xscan ):
   
     # separation --> taken from pickle file
     # current    --> pickle file
     # bunches    --> pcikle file
     # mean       --> taken from the single gauss fit (p1)
     # width      --> taken from the single gauss fit (p2)
     # beta       --> contant 400
     # isX        --> tunning (confirm defination)??

     # single beam deflection angle
    single_beam_deflection_angle = beta  * radius / gamma * (current / bunches ) * 1e11 / math.tan( math.pi * ( 64.31 if Xscan else 59.32 ) )
    
    # single beam displacement
    single_beam_displacement = 2 * single_beam_deflection_angle / separation * ( 1.0 - gauss( separation, 1.0, mean,width ) )
    
    # double beam displacement
    # [0] index is for the first beam
    # [1] index for the second beam
    double_beam_displacement = np.array(single_beam_displacement[0]) + single_beam_displacement[1]
   
    # returns the separation and correction addition
    return separation + double_beam_displacement

    # loading pickle file    
fit_pickle = sorted(glob.glob("fits/integrated/fit_integrated_trk_run_364030.pickle"))

    # 3.main function (analyze function)
def analyse_scan(separation, luminosity, error, path, run_number,bunches,current,Xscan):

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
    # print p0
    p1 = single_gauss_fit.GetParameter(1)
    # print p1
    p2 = single_gauss_fit.GetParameter(2)
    # print p2

    # calling beam_beam deflection function for getting the parameters of it
    # beam_beam_deflection gives us separation plus correction

    corrected_separation = beam_beam_deflection(separation, current, bunches, p1, p2, 400, Xscan)
    # Tgraph with corrected_separation
    #graph_corrected = ROOT.TGraphErrors(separation_points, corrected_separation, luminosity, np.zeros(separation_points),error )

    # setting parameters for double gauss
    double_gauss_fit.SetParameter(0,p0)#1st parameter 
    double_gauss_fit.SetParameter(1,p1)#2nd parameter
    double_gauss_fit.SetParameter(2,p2)#3rd parameter
    
    # second method of getting single gauss parameters
    #for i in range(3):            
    #    double_gauss_fit.SetParameter(i, single_gauss_fit.GetParameter(i))#get parameters of double gaussain
        
    # double gauss setting parameters
    double_gauss_fit.SetParameter(3,1.0)#4rd parameter                                                                      
    double_gauss_fit.SetParameter(4,1.0)#5th parameter                                                                             
    # fit double gauss                                                        
    graph.Fit('double_gauss_fit','Q')#fitting of double gaussain
    
    # ****** without corrected separation *****    

    # chi2 value
    # luminsoity from fit (each point)
  
    chi2_uncorrected = np.array([])
    lumi_fit = gauss2D([separation], double_gauss_fit.GetParameters())
    
    # looping on all separation points
    # to get lumionsity of each separation points
    # lastly calculate the chi2

    for sep in range(separation_points):  
        chi2 = (luminosity[sep]- lumi_fit[sep])**2/error[sep]# chi2 formula
        chi2_uncorrected = np.append(chi2_uncorrected,chi2)

    # chi2 from fit before correction
    chi2_double_gauss_fit = double_gauss_fit.GetChisquare()
    NDF_double_gauss_fit =  double_gauss_fit.GetNDF()

    # corrected_separation = beam_beam_deflection(separation, current, bunches, p1, p2, 400, Xscan)
    # # Tgraph with corrected_separation
    # graph_corrected = ROOT.TGraphErrors(separation_points, corrected_separation, luminosity, np.zeros(separation_points),error )

    # print statements

    #float
    chi2_double_gauss_fit == True
    print chi2_double_gauss_fit,"chi from fit"

    #float
    NDF_double_gauss_fit == True
    print NDF_double_gauss_fit,  "NDF From double gauss fit" 
   
    #ndarray
    chi2_uncorrected.any()==True
    print chi2_uncorrected ,'chi values from formula'

    #ndarray
    separation.any()==True
    print separation ,"separation"

    #ndarray
    corrected_separation.any()== True
    print corrected_separation,"corrected Separation"

    # grap for the 
    graph.Draw("AP")
    canvas_name.Print(path)

for filename in fit_pickle:
    with open(filename,'rb') as fit_pickle:
        scans = pickle.load(fit_pickle)
        #print "Number of scans", len(scans)
        for i_scan in range(len(scans)):
            output = scans[i_scan]

            analyse_scan(output['sep_x'],output['exp_lumi_x'],output['err_x'],"Emittance_scan/"+str(output['run'])+"_"+str(i_scan)+"_"+"x_scan.png",output['run'],output['bunches'],output['curr'][0],True)



# single gauss fit
    # def gauss(x, p):                                                                                                           
    #     return p[0] * np.exp( - 0.5*( (x[0]-p[1])**2 / (2.*p[2]**2) ) )     
