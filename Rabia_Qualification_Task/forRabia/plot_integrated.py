import ROOT
import glob
import pickle
import numpy as np

import csv
scans_eric = {}
csvfile = open('scanTable.csv', 'rb')
for row in csvfile:
    if row[0] == "#":
        continue
    row = row.split(",")
    row = [ item.strip() for item in row ]
    if int( row[2] ) not in scans_eric.keys():
        scans_eric[int( row[2] )] = list()
    scans_eric[ int( row[2] ) ].append( [ float(row[9]), float(row[14]), float(row[16]), float(row[10]), float(row[15]), ] )


fits = sorted(glob.glob("fits/integrated/*lcd*pickle"))
results_lcd = {}
for name in fits:
    with open( name, "rb" ) as picl:
        data = pickle.load( picl )
        results_lcd[data[0]['run']] = data
        print data                
fits = sorted(glob.glob("fits/integrated/*trk*pickle"))
results_trk = {}
for name in fits:
    with open( name, "rb" ) as picl:
        data = pickle.load( picl )
        results_trk[data[0]['run']] = data

def plot( path, name, x, y, low = 0.98, high = 1.02 ):
    canvas = ROOT.TCanvas( name, name, 800, 600 ) 
    mgr = ROOT.TMultiGraph()
    for e in [0, 1]:
        gr = ROOT.TGraph( len( x[e] ),
                          x[e], y[e])
        if e:
            gr.SetMarkerStyle( ROOT.kFullCircle )
        else:
            gr.SetMarkerStyle( ROOT.kOpenCircle )
        gr.SetMarkerColor( ROOT.kBlue )
        mgr.Add( gr, "P" )
    mgr.Draw("A pmc")
    mgr.SetMinimum( low )
    mgr.SetMaximum( high )
    canvas.Modified()
    canvas.Print( path )

peak     = [ np.array([ results_trk[x][i]['peak'] for x in results_trk for i in range(len(results_trk[x])) if results_trk[x][i]['early'] == e ]) for e in [True, False] ]
fom_trk  = [ np.array([ results_trk[x][i]['fom']  for x in results_trk for i in range(len(results_trk[x])) if results_trk[x][i]['early'] == e ]) for e in [True, False] ]
plot( "neo/foms_trk.png", "foms_trk", peak, fom_trk )

peak     = [ np.array([ results_trk[x][i]['peak'] for x in results_trk for i in range(len(results_trk[x])) if results_trk[x][i]['early'] == e ]) for e in [True, False] ]
fom_lcd  = [ np.array([ results_lcd[x][i]['fom']  for x in results_trk for i in range(len(results_trk[x])) if results_trk[x][i]['early'] == e ]) for e in [True, False] ]
plot( "neo/foms_lcd.png", "foms_lcd", peak, fom_lcd )

peak       = [ np.array([ results_trk[x][i]['peak'] for x in results_trk for i in range(len(results_trk[x])) if results_trk[x][i]['early'] == e ]) for e in [True, False] ]
fom_ratio  = [ np.array([ results_trk[x][i]['fom']/results_lcd[x][i]['fom']  for x in results_trk for i in range(len(results_trk[x])) if results_trk[x][i]['early'] == e ]) for e in [True, False] ]
plot( "neo/foms_ratio.png", "foms_ratio", peak, fom_ratio )

peak           = [ np.array([ results_trk[x][i]['peak'] for x in results_trk for i in range(len(results_trk[x])) if results_trk[x][i]['early'] == e ]) for e in [True, False] ]
width_x_ratio  = [ np.array([ results_trk[x][i]['width'][0]/results_lcd[x][i]['width'][0]  for x in results_trk for i in range(len(results_trk[x])) if results_trk[x][i]['early'] == e ]) for e in [True, False] ]
plot( "neo/width_x_ratio.png", "width_x_ratio", peak, width_x_ratio )

peak           = [ np.array([ results_trk[x][i]['peak'] for x in results_trk for i in range(len(results_trk[x])) if results_trk[x][i]['early'] == e ]) for e in [True, False] ]
width_y_ratio  = [ np.array([ results_trk[x][i]['width'][1]/results_lcd[x][i]['width'][1]  for x in results_trk for i in range(len(results_trk[x])) if results_trk[x][i]['early'] == e ]) for e in [True, False] ]
plot( "neo/width_y_ratio.png", "width_y_ratio", peak, width_y_ratio )

peak         = [ np.array([ results_trk[x][i]['peak'] for x in results_trk for i in range(len(results_trk[x])) if results_trk[x][i]['early'] == e ]) for e in [True, False] ]
peak_ratio  = [ np.array([ results_trk[x][i]['peak']/results_lcd[x][i]['peak']  for x in results_trk for i in range(len(results_trk[x])) if results_trk[x][i]['early'] == e ]) for e in [True, False] ]
plot( "neo/peak_ratio.png", "peak_ratio", peak, peak_ratio )


peak         = [ np.array([ results_trk[x][i]['peak'] for x in results_trk for i in range(len(results_trk[x])) if results_trk[x][i]['early'] == e ]) for e in [True, False] ]
width_param  = [ np.array([ results_trk[x][i]['params'][0][4]  for x in results_trk for i in range(len(results_trk[x])) if results_trk[x][i]['early'] == e ]) for e in [True, False] ]
plot( "neo/width_param_x.png", "width_param_x", peak, width_param, 0.0, 2.0 )

peak         = [ np.array([ results_trk[x][i]['peak'] for x in results_trk for i in range(len(results_trk[x])) if results_trk[x][i]['early'] == e ]) for e in [True, False] ]
width_param  = [ np.array([ results_trk[x][i]['params'][1][4]  for x in results_trk for i in range(len(results_trk[x])) if results_trk[x][i]['early'] == e ]) for e in [True, False] ]
plot( "neo/width_param_y.png", "width_param_y", peak, width_param, 0.0, 2.0 )

peak         = [ np.array([ results_trk[x][i]['peak'] for x in results_trk for i in range(len(results_trk[x])) if results_trk[x][i]['early'] == e ]) for e in [True, False] ]
split_param  = [ np.array([ results_trk[x][i]['params'][0][3]  for x in results_trk for i in range(len(results_trk[x])) if results_trk[x][i]['early'] == e ]) for e in [True, False] ]
plot( "neo/split_param_x.png", "split_param_x", peak, split_param, 0.0, 1.0 )

peak         = [ np.array([ results_trk[x][i]['peak'] for x in results_trk for i in range(len(results_trk[x])) if results_trk[x][i]['early'] == e ]) for e in [True, False] ]
split_param  = [ np.array([ results_trk[x][i]['params'][1][3]  for x in results_trk for i in range(len(results_trk[x])) if results_trk[x][i]['early'] == e ]) for e in [True, False] ]
plot( "neo/split_param_y.png", "split_param_y", peak, split_param, 0.0, 1.0 )

canvas = ROOT.TCanvas("chi2", "chi2", 800, 600)
hist = ROOT.TH1F( "chi2_x", "chi2_x", 20, 0, 100 )
for x in results_lcd:
    for i in range(len(results_trk[x])):
        hist.Fill( results_trk[x][i]['chi2'][0] )
        print results_trk[x][i]['chi2']
hist.Draw()
canvas.Print("neo/chi2_x.png")

canvas = ROOT.TCanvas("chi2", "chi2", 800, 600)
hist = ROOT.TH1F( "chi2_y", "chi2_y", 20, 0, 100 )
for x in results_lcd:
    for i in range(len(results_trk[x])):
        hist.Fill( results_trk[x][i]['chi2'][1] )
        print results_trk[x][i]['chi2']
hist.Draw()
canvas.Print("neo/chi2_y.png")
