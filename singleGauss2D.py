from Functions import *
import ROOT
#import EmittanceScan_2d
import numpy as np

def fit2D( sep_x,sep_y,lumi,error_lumi,peak,width_x,width_y,Run,scan,path,R1,R2,R3,R4):

    canvas_2d = ROOT.TCanvas( 'fit', 'fit',800, 600 )
    graph_2d = ROOT.TGraph2DErrors(len(sep_x + sep_y), sep_x, sep_y,lumi, np.zeros(len(sep_x)),np.zeros(len(sep_y)), error_lumi)
    fit_2d = ROOT.TF2("fit_2d", single_gauss_2d_fit, R1, R2,R3,R4, 5)    #2d gauss (6 parameters)                                                                                                                                             

    fit_2d.SetParameter(0 ,peak)

    # mean1(x)                                                                                                                                                                                                                                
    fit_2d.SetParameter(1 ,0.0)
    # mean2(y)

    fit_2d.SetParameter(3 , 0.0)

    # width(x)                                                                                                                                                                                                                                
    fit_2d.SetParameter(2 , width_x)
    
    # width(y)                                                                                                                                                                                                                                
    fit_2d.SetParameter(4 , width_y)

    # increaing the number of iterations 
    fit_2d.SetParNames("Amplitude","MeanX","SigmaX", "MeanY","SigmaY")                                                                                                                            


    
    ROOT.TVirtualFitter.Fitter(graph_2d).SetMaxIterations(10000)

    # fitting                                                                                                                                                                                                                                 
    graph_2d.Fit('fit_2d',"RQ")#fitting of double gaussain                                                                                                                                                                                   

    graph_2d.Draw("P")
    fit_2d.Draw("LINE")
#    canvas_2d.Draw("")
    canvas_2d.Print(path + '_' + str(Run)+'_'+str(scan)+ '.png')
#canvas_2d.Print('' + '_' + str(run)+'_'+str(i_scan)+ '.png')                                                                                                                                                                                 
    NDF=fit_2d.GetNDF()
    # print fit_2d.GetParameter(0),"peak"
    # print fit_2d.GetParameter(1),"mean1"                                                                          


    # print fit_2d.GetParameter(3),"mean2"



    # print fit_2d.GetParameter(2),"width x"
    # print fit_2d.GetParameter(4),"widthy"
    
    chi2_NDF = fit_2d.GetChisquare() / NDF
    peak_2d = fit_2d.GetMaximum()
    print peak_2d,"peak"

    return chi2_NDF, peak_2d
