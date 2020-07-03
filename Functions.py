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



def diagnol(amp1,amp2,sigma1,sigma2,mean1,mean2,amp1_E,amp2_E,sigma1_E,sigma2_E,mean1_E,mean2_E):


    diagnol_terms = ( amp1_E * amp1 )**2 + ( amp2_E * amp2 )**2 + ( sigma1_E * sigma1 )**2 + ( sigma2_E * sigma2 )**2 + ( mean1_E * mean1 )**2 + ( mean2_E * mean2 )**2

    return diagnol_terms




def derivative_value(derivative,amp1,sig1,mea1,amp2,sig2,mea2,position):

    
    required_value = derivative.evalf(subs={amplitude1:amp1,sigma1:sig1,mean1:mea1,amplitude2:amp2,sigma2:sig2,mean2:mea2, x:position})


    return required_value




# fit function

def one_dimension_double_gauss_fit(x, p):

    return p[0] * (1.0 / ( math.sqrt(2 * math.pi) * p[2] ) ) * np.exp( - (1.0 / 2.0) * ( (x[0] - p[1])**2 / ( p[2] **2) ) ) + p[3] * (1.0 / (  math.sqrt(2 * math.pi) * p[4] ) ) * np.exp( - (1.0 / 2.0) * ( (x[0] - p[5])**2 / ( p[4] **2) ) )







# expected _luminosity fucntion                                                                                               
frequency = (1./8.89244e-05) #[s-1]                                                                                                
inel = 8. * 1e-24

def expected_luminosity( sigma_x, sigma_y, current ):
    return frequency * current * 1e22 / ( sigma_x * sigma_y  * 2.0 * math.pi) / (frequency / inel )






### graph plotting function                                                                                                     plot_graph ("plots/FOM_date.png",Number, date, FOM, FOM_error,"fom","date",True)













def graph (path_graph, name_graph, number_graph, x, y,y_error, graph_title_x,graph_title_y, isdate = False):

    canvas = ROOT.TCanvas( 'name', 'name',800, 600 )

    graph = ROOT.TGraphErrors(number_graph , x, y, np.zeros(len(x)), y_error)

    graph.SetMarkerColor( 2 )
    graph.SetMarkerStyle( 20 )
    graph.SetMarkerSize( 1.3 )
    graph.GetXaxis().SetTitle(graph_title_x)

    if isdate:                                                                                                                
        graph.GetXaxis().SetTimeDisplay(1)                                                                                 
        graph.GetXaxis().SetTimeFormat("%d/%m")       

    graph.GetYaxis().SetTitle(graph_title_y)
    graph.GetXaxis().SetNdivisions(4,4,0)
    graph.Draw("AP")
    canvas.Modified()
    canvas.Print( path_graph )






               
def plot_graph (path_graph, number_graph, x, y, graph_title_x,graph_title_y):#, isdate = False):
# isdate= False, islowMu = False):                                                                                                             

    canvas = ROOT.TCanvas( 'name', 'name',800, 600 )
    graph = ROOT.TGraph(number_graph , x, y)

#    graph = ROOT.TGraphErrors(number_graph , x, y, np.zeros(len(x)), y_error)

    graph.SetMarkerColor( 2 )
    graph.SetMarkerStyle( 20 )
    graph.SetMarkerSize( 1.3 )
    graph.GetXaxis().SetTitle(graph_title_x)

    # if isdate:                                                                                                                               
    #     graph.GetXaxis().SetTimeDisplay(1)                                                                                                   
    #     graph.GetXaxis().SetTimeFormat("%d/%m")                                                                                 
        
    graph.GetYaxis().SetTitle(graph_title_y)

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





    
def graph_cosmetics(required_graph,color,markersize,XTitle,YTitle):


    required_graph.SetLineColor( 38 )
    required_graph.SetMarkerColor( color )
    required_graph.SetMarkerStyle( 20 )
    required_graph.SetMarkerSize( markersize )
    required_graph.GetXaxis().SetTitle(XTitle)
    required_graph.GetYaxis().SetTitle(YTitle)
    required_graph.GetXaxis().SetTitleSize(0.064)
    required_graph.GetYaxis().SetTitleSize(0.064)




def legend(entry_1,entry_2):
    leg = ROOT.TLegend(0.75, 0.6, 0.9, 0.86)
    leg.AddEntry(entry_1,"data","P")# points                                                                                   
    leg.AddEntry(entry_2," Fit","L")#                                                                                          
    leg.Draw()





def graph_2_canvas_1(graph_1,graph_2,title1_x,title1_y,title2_x,title2_y):

    canvas = ROOT.TCanvas('sub_lumi', 'sub_lumi')
    canvas.Divide(1,2) # 1 row, 2 columns                                                                                     

    canvas.cd(1)
    graph_cosmetics(graph_1,4,1.2,title1_x,title1_y)
#    legend(entry1,entry2)
    graph_1.Draw("AP")


    canvas.cd(2)
    graph_cosmetics(graph_2,4,1.2,title2_x,title2_y)
    graph_2.Draw("AP")
                                                                                                                                  
    return canvas






def parameter_info(fit,parameter):

    value=fit.GetParameter(parameter)
    parameter_error=fit.GetParError(parameter)
    relative_error = (parameter_error/value)*100
    return value,parameter_error,relative_error




def lengend(l,m,n,o,entry1,label11,label12,leg_size,latex,color):
    leg=ROOT.TLegend(l,m,n,o)
    leg.AddEntry(entry1,label11,label12)
#    leg.AddEntry(entry1,label11,label12)
    leg.SetTextSize(leg_size)
    label = ROOT.TLatex()
    label.SetTextSize(label_size)
    label.DrawLatexNDC(0.75,0.90,latex)
    label.SetTextColor(color)
    leg.Draw()


