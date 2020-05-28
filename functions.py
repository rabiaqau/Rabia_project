


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

def intensity(current1, current2, bunches):

      n1= np.mean(current1)/bunches
      n2= np.mean(current2)/bunches
      intensity = n1*n2

      return intensity


def analyse( input_data ):
    sep_x =  input_data['sep_x']
    sep_y =  input_data['sep_y']
    exp_lumi_x = input_data['exp_lumi_x']
    exp_lumi_y = input_data['exp_lumi_y']
    err_x = input_data['err_x']
    err_y = input_data['err_y']
    curr_0 = input_data['curr'][0]
    curr_1 = input_data['curr'][1]
    date =  input_data['date']
    bunches = input_data['bunches']
    early = input_data['early']
    run_number = input_data['run']

# calling function analyze                                                                                                        
# X scan true                                                                                                                     
# here i can have x values early and late                                                                                         
    result_1 = analyse_scan(sep_x, exp_lumi_x ,err_x ,"Emittance_scan/", run_number, bunches, curr_0, True, early,date )

# y scan true                                                                                                                     
# here i can have y values early and late                                                                                         

    result_2 =  analyse_scan(sep_y, exp_lumi_y ,err_y ,"Emittance_scan/", run_number, bunches, curr_1, False , early, date )
    return( result_1, result_2,early,bunches, curr_0,curr_1 )






# plotting fucntion
def plot_graph (path_graph, name_graph, number_graph, x, y, graph_title_x,graph_title_y, islowMu = False):
# isdate= False, islowMu = False):                                                                                                             

    canvas = ROOT.TCanvas( 'name', 'name',800, 600 )
    graph = ROOT.TGraph(number_graph , x, y)

    graph.SetMarkerColor( 2 )
    graph.SetMarkerStyle( 20 )
    graph.SetMarkerSize( 2 )
    graph.GetXaxis().SetTitle(graph_title_x)
    # if isdate:                                                                                                                               
    #     graph.GetXaxis().SetTimeDisplay(1)                                                                                                       #     graph.GetXaxis().SetTimeFormat("%d/%m")                                                                                              
    graph.GetYaxis().SetTitle(graph_title_y)
    if islowMu:
        graph.GetXaxis().SetRangeUser(0,5)
    graph.GetXaxis().SetNdivisions(4,4,0)
    graph.Draw("AP")
    canvas.Modified()
    canvas.Print( path_graph )



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



def plot_graph_erro (path_graph, name_graph, number_graph, x, y,x_err,y_err, graph_title_x,graph_title_y, islowMu = False):
# isdate= False, islowMu = False):                                                                                                             

    canvas = ROOT.TCanvas( 'name', 'name',800, 600 )
    print len(x),len(y),len(y_err),number_graph
    graph = ROOT.TGraphErrors(number_graph , x, y,x_err, y_err)
    
    graph.SetMarkerColor( 2 )
    graph.SetMarkerStyle( 20 )
    graph.SetMarkerSize( 0.5 )
    graph.GetXaxis().SetTitle(graph_title_x)

#    graph.GetYaxis().SetRangeUser(0.0,100)
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



def plot_graph_error (path_graph, name_graph, number_graph, x, y,y_err, graph_title_x,graph_title_y, islowMu = False):
# isdate= False, islowMu = False):                                                                                                             

    canvas = ROOT.TCanvas( 'name', 'name',800, 600 )
    print len(x),len(y),len(y_err),number_graph
    graph = ROOT.TGraphErrors(number_graph , x, y,np.zeros(len(x)), y_err)
    
    graph.SetMarkerColor( 2 )
    graph.SetMarkerStyle( 20 )
    graph.SetMarkerSize( 1 )
    graph.GetXaxis().SetTitle(graph_title_x)

#    graph.GetYaxis().SetRangeUser(0.0,100)
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










def multigraph(path, number1, number2, x1, y1,x2,y2):
    canvas = ROOT.TCanvas( 'scan', 'scan', 800, 600 )
    ROOT.TGaxis.SetMaxDigits(6)
    graph1 = ROOT.TGraph(number1 ,x1 ,y1)
    graph1.SetMarkerColor(ROOT.kRed )
    graph1.SetMarkerStyle(20)
    graph2 = ROOT.TGraph(number2 ,x2 ,y2 )
    graph2.SetMarkerColor(ROOT.kBlue)
    graph2.SetMarkerStyle(4)
    graph1.GetXaxis().SetNdivisions(4,4,0)
    graph2.GetXaxis().SetNdivisions(4,4,0)
    mg = ROOT.TMultiGraph()
    mg.Add(graph1)
    mg.Add(graph2)
    mg.Draw("ap")

    mg.GetXaxis().SetRangeUser(-1,2)
    mg.GetXaxis().SetLimits(-1,2)

    canvas.Modified()
    canvas.Print(path )






# fits
def gauss2D(x, p):                                                                                                                           
       return p[0] * ( np.exp( - (x[0]-p[1])**2/(2.*p[2]**2) ) * np.exp( - (x[1]-p[3])**2/(2.*p[4]**2)  ) )    
#
def gaussD(x, p):

    return p[0] * (1.0 / ( math.sqrt(2 * math.pi) * p[2] ) ) * np.exp( - (1.0 / 2.0) * ( (x[0] - p[1])**2 / ( p[2] **2) ) ) + p[3] * (1.0 / (math.sqrt(2 * math.pi) * p[4] ) ) * np.exp( - (1.0 / 2.0) * ( (x[0] - p[5])**2 / ( p[4] **2) ) )


freq = (1./8.89244e-05) #[s-1]                                                                        
inel = 8. * 1e-24

def expected_luminosity( sigma_x, sigma_y, current ):
    return freq * current * 1e22 / ( sigma_x * sigma_y  * 2.0 * math.pi) / (freq / inel )

