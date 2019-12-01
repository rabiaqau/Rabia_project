
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
fit_pickle = sorted(glob.glob("fits/integrated/*trk*pickle"))
for filename in fit_pickle:
#def main (filename):
    with open(filename,'rb') as fit_pickle:
        print filename
        Data = pickle.load(fit_pickle)[0]
        print Data
        sep_x = Data['sep_x']
        print sep_x
        lumi_x = Data['exp_lumi_x']
        #for index, value in enumerate(lumi_x):
            #print(index, value)
        Data_err = Data['err_x']#*100
        #print Data_err                                                                                                                                                                                  
        bunches = Data['bunches']
        no_of_sep= len(sep_x)
####percentage of error 
        # for lumi in lumi_x:
        #     #print "Luminosity", lumi
        # for error in Data_err:
        #     #print "error", error
        # def func1(error,luminosity):
        #     return (error/luminosity)*100
        # #print(func1(Data_err,lumi_x)),"percentage of error"

        # scanning data                                                                                                                                                                            
        graph = ROOT.TGraph( no_of_sep,sep_x,lumi_x)#TGraph--> class of grphical object on x and y array
        graph = ROOT.TGraphErrors( no_of_sep, sep_x, lumi_x, np.zeros(no_of_sep) , Data_err )
        canvas = ROOT.TCanvas( 'x-scan', 'x-scan',800, 600 )#TCanvas-->class for graphical display
        #single guusain function
        #gauusain fit from alex emittance_tools.py file
        def gauss(x, p): # x is vector(x-->x[0],x[1]), p is a vector-->p[0],p[1],p[2]           
            
            return p[0] * np.exp( - 0.5*( (x[0]-p[1])**2 / (2.*p[2]**2) ) )
        #####2d gusain function (5 parameters)
        def gauss2D(x, p):

            return p[0] / math.sqrt(2 * math.pi) / p[2] * ( p[3]  * np.exp( - ( (x[0]-p[1])**2 / (2.*p[2]**2) ) ) + ( 1.0 - p[3] ) * p[4] *np.exp( - ( (x[0]-p[1])**2/(2.*(p[2]/p[4])**2) ) ) )
       ## parameters setting 
        fit_3_parameters = ROOT.TF1( "fit_3_parameters", gauss, -3, 3, 3 )# 1d gauss
        fit_5_parameters = ROOT.TF1("fit_5_parameters", gauss2D, -3, 3, 5)#2d gauss
        fit_3_parameters.SetParameter(0,20)#p[0]-->mu initial value
        fit_3_parameters.SetParameter(1,0)#p[1]-->center
        fit_3_parameters.SetParameter(2,1) #p[2}-->sigma        
        #fitting 
        graph.Fit('fit_3_parameters','Q')
#supplying the parameters of single gaussain to double gaussain(idea of alex)
        for i in range(3):            
            fit_5_parameters.SetParameter(i, fit_3_parameters.GetParameter(i))#get parameters of double gaussain
    #    print (func2.SetParameter(i, func1.GetParameter(i)))
        fit_5_parameters.SetParameter(3,1.0)#3rd parameter
        fit_5_parameters.SetParameter(4,1.0)#4th parameter
        graph.Fit('fit_5_parameters','Q')#fitting of double gaussain
        
        for values in range(9):
           #print gauss2D(sep_x[values],lumi_x[values])
        #print lumi_x, "lumi"
        #print gauss2D([sep_x], fit_5_parameters.GetParameters()),"gauss"
        #print Data_err
            chi = (lumi_x[values]- gauss2D([sep_x[values]], fit_5_parameters.GetParameters()))/Data_err[values]
            print chi ,'chi values'

        # for i in range(9):
        #     print lumi_x[i]
        # for j in range(9):
        #     print gauss2D(lumi_x[j], sep_x[j])
#             print "difference", ( lumi_x[i] - gauss2D(x[0],lumi_x[i]) 
# #additnal info from the fit
#         chi_single_gauss= fit_3_parameters.GetChisquare()
#         NDF_single_gauss= fit_3_parameters.GetNDF()
#         print "chi_single_gauss:",chi_single_gauss
#         print "Effective_NDF_single_gauss(9-3):",NDF_single_gauss
#         Ratio_single_gauss = chi_single_gauss / NDF_single_gauss
#         print "Ratio_single_gauss:",Ratio_single_gauss
# #for double gaussian
        # chi_double_gauss= fit_5_parameters.GetChisquare()
        # NDF_double_gauss= fit_5_parameters.GetNDF()
        # print "chi_double_gauss:",chi_double_gauss
        # print "Effective_NDF_double_gauss(9-5):",NDF_double_gauss
        # Ratio_double_gauss = chi_double_gauss / NDF_double_gauss
        # print "Ratio_double_gauss:",Ratio_double_gauss
        # print lumi_x,"lumi"

        # def test(x, y):
        #     return (x - y)
            
        # print (test(lumi_x,gauss2D(x[0],sep_x)))                
####test
#         def gauss2D(x, p):
#             return p[0] / math.sqrt(2 * math.pi) / p[2] * ( p[3]  * np.exp( - ( (x[0]-p[1])**2 / (2.*p[2]**2) ) ) + ( 1.0 - p[3] ) * p[4] *np.exp( - ( (x[0]-p[1])**2/(2.*(p[2]/p[4])**2) ) ) )

#         func2 = ROOT.TF1("func2", gauss2D, -3, 3, 5)#2d gauss                                                                                   
#         func2.SetParameter(0,20)#p[0]-->constant                                                                                                
#         func2.SetParameter(1,0)#p[1]-->center                                                                                                   
#         func2.SetParameter(2,1)
#         func2.SetParameter(3,1.0)#3rd parameter                                                                                                 
#         func2.SetParameter(4,1.0)#4th parameter                                                                                                 
#         gr.Fit('func2','')#fitting of double gaussain
# #for double gaussian                                                                                                                            
#         chi_double_gauss= func2.GetChisquare()
#         NDF_double_gauss= func2.GetNDF()
#         print "chi_double_gauss:",chi_double_gauss
#         print "Effective_NDF_double_gauss(9-5):",NDF_double_gauss
#         Ratio_double_gauss = chi_double_gauss / NDF_double_gauss
#         print "Ratio_double_gauss:",Ratio_double_gauss
###########
        #gr.Fit('gaus','')#for built in fuction of the root
        graph.Draw("AP")
        canvas.Print( "1d_pickle.png" )

#if __name__=="__main__":
#    filename = str(sys.argv[1])
#    main(filename)


