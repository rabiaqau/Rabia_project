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
#from emittance_scan import *
#from emittanceScan import analyse_scan




### graph plotting function                                                                                                                    
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




## histogram plotting function                                                                                                                                                                              


def histogram(bin,low,up,fill,Name,hist_path,hist_title):

    hist = ROOT.TH1F( "Name", "Name", bin, low, up)

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

    hist_leg.AddEntry(0,'Overflow='+str(max_bin), '')

    hist_leg.AddEntry(0,'Underflow='+str(min_bin), '')

    hist_leg.Draw()


    canvas_hist.Print(hist_path)




def comp_histogram(bin,low,up,fill1,fill2,hist_path):

    canvas_hist= ROOT.TCanvas("comparision", "two comparision", 800, 600)

    hist1 = ROOT.TH1F( "h1", "h1", bin, low, up)


    for i in fill1:
        hist1.Fill(i)

    hist1.Draw("h")

    canvas_hist.Update()


    hist2 = ROOT.TH1F( "h2", "h2", bin, low, up)

    for j in fill2:
        hist2.Fill(2)

    hist2.Draw("pe,same")
    
    canvas_hist.Print(hist_path)













def initialparameter(fit, Amp1, Amp2, Mean1, Mean2 ,Sigma1, Sigma2 ):

    fit.SetParameter(0,Amp1)#amplitude                                                                                        
    fit.SetParameter(3,Amp2)#amp                                                                                              
    fit.SetParameter(4,Sigma2)#sigma2                                                                                         
    fit.SetParameter(1,Mean1)#mean1                                                                                           
    fit.SetParameter(2,Sigma1)#sigma1                                                                                         
    fit.SetParameter(5,Mean1)#mean                                                                                            

# if run_number==350013:
#     if early and not Xscan:
#         initialparameter(double_gauss_fit, 0.3*amp1, 0.9*amp1, m1, m1, s1, 1.3*s1)
#     if not early and Xscan:
#         initialparameter(double_gauss_fit, 0.711, 0.24, -m1*0.1, m1*0.5, 0.0015, 0.5)
#    # <2000                                                                                                                                     

# if run_number==359286:
#     if Xscan:
#         initialparameter(double_gauss_fit, 0.6*amp1, 0.4*amp1, m1, m1, 0.9*s1, 1.3*s1)
#     if not Xscan:
#         initialparameter(double_gauss_fit, 0.5*amp1, 0.5*amp1, m1, m1, s1, 2*s1)







# def initialparameter(run, A1, A2, M1, M2 ,S1, S2 , early, Xscan):


#     if run_number== run:
        
#         double_gauss_fit.SetParameter(0,A1)#amplitude 
#         double_gauss_fit.SetParameter(3,A2)#amp
#         double_gauss_fit.SetParameter(4,S2)#sigma2                                                                 
#         double_gauss_fit.SetParameter(1,M1)#mean1
#         double_gauss_fit.SetParameter(2,S1)#sigma1  
#         double_gauss_fit.SetParameter(5,M1)#mean   



#initialparameter(double_gauss_fit,359286,0.6*amp1, 0.4*amp1, m1, m1, 0.9*s1, 1.3*s1, True, True)
