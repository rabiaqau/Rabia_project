#!/usr/bin/env python
# coding: utf-8

from collections import namedtuple
from emittance_tools import *
from functools import reduce
from matplotlib import rc
from numpy import inf
from numpy import nan
from operator import itemgetter
from scipy.interpolate import UnivariateSpline

import AtlasStyle
import AtlasUtils
import ROOT
import glob
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os.path
import pickle
import root_numpy as rnp
import scipy.integrate as integrate
import scipy.optimize
import sys
import traceback

class NoScans(Exception):
    """Run does not have any emittance scans"""
    def __init__(self):
        pass

class FitError(Exception):
    """Error during emittance scan fit"""
    def __init__(self):
        pass

# List of working points
WPs   = [ 'Tight', 'TightLumi', 'TightModLumi','TightModSiPlusLumi' ]

# Calibration constants for the tracking
calib = {  'Tight'             : 3.8573,
           'TightLumi'         : 3.6961,
           'TightModLumi'      : 1.68517,
           'TightModSiPlusLumi': 1.66573 }

# Constants useful for the beam-beam deflection potential
freq = (1./8.89244e-05) #[s-1]  
inel = 8. * 1e-24       #[mm2]
gamma = 13e6 / 2. / 938
radius = 1.5347e-15     #[mm]

# Import the reference values from Eric Torrence 
#  The are in format: X width, Y width, FoM
#
# The full CSV (Comma Separated Values) file column structure:
#    |  0 | Scan      |
#    |  1 | Fill      |
#    |  2 | Run       |
#    |  3 | Bunches   |
#    |  4 | StartTime |
#    |  5 | EndTime   |
#    |  6 | XSteps    |
#    |  7 | XPeak     |
#    |  8 | XMean     |
#    |  9 | XWidth    |
#    | 10 | XChi2     |
#    | 11 | YSteps    |
#    | 12 | YPeak     |
#    | 13 | YMean     |
#    | 14 | YWidth    |
#    | 15 | YChi2     |
#    | 16 | LRatio    |
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
    scans_eric[ int( row[2] ) ].append( [ float(row[9]), float(row[14]), float(row[16]), ] )

# Correction for beam-beam deflection, input is the separation of each
# point, the current, number of bunches, fit parameters and beta. Need
# to know whether it is an x- or y-scan to use the correct beam tune
def beam_beam_deflection( sep, curr, bunches, mean, width, beta, isX ):
    factor = beta * 10. * radius / gamma * ( curr / bunches ) * 1e11 / math.tan( math.pi * ( 64.31 if isX else 59.32 ) )
    with np.errstate(divide='ignore'):
        corr = 2 * factor / sep * ( 1.0 - gauss( sep, 1.0, mean, width ) )
    corr[ corr == inf ] = 0
    corr[ corr == nan ] = 0
    corr = numpy.nan_to_num( corr )
    return sep + corr

# Calculate expected luminosity, this is in terms of mu for a single bunch
def exp_lumi( bunches, width_x, width_y ):
    return freq * 2 / 2 * 1e22 / ( width_x * width_y  * 2.0 * math.pi) / (freq / inel )

dgauss = ROOT.TF1( "DoubleGauss", gaussD, -0.1, 0.1, 5 )
gauss2d = ROOT.TF2( "Gauss2D", gauss2D, -0.1, 0.1, -0.1, 0.1, 5 )

def scanFit( x, y, err, intprod_vec, intensity, bunches, beta, isX, isEarly ):
    # GUARDS!
    intprod_vec[ np.isnan( intprod_vec ) ] = 1.0
    intprod_vec[ np.isinf( intprod_vec ) ] = 1.0 
    intprod_vec = intprod_vec / bunches / bunches
    
    # Compute the average and standard deviation of the scan, for
    # setting initial values for the fit
    avg = numpy.average( x, weights=y/intprod_vec )
    std = math.sqrt( numpy.average( (x-avg)**2, weights=y/intprod_vec) )
    
    # Limit the parameters to avoid bad fits
    # Fix the second gaussian when doing the first fit
    # Simplifies to a simple gaussian fit
    #
    # Fit parameters are: 
    #
    # 0: Overall amplitude (constant factor in the front, this is NOT the amplitude of the dual gaussian curve!
    # 1: The mean of gaussians
    # 2: The standard deviation of the first gaussian
    # 3: The fraction of the overall amplitude taken for the first gaussian (second gets (1 - par[3]) )
    # 4: The factor by which the second gaussian is WIDER (i.e. sigma2 = par[4] * par[2])
    dgauss.SetParameter( 0, numpy.max( y / intprod_vec ) )
    dgauss.SetParLimits( 0, 0.0, 100.0 )
    dgauss.SetParameter( 1, avg )
    dgauss.SetParLimits( 1, -0.1, 0.1 )
    dgauss.SetParameter( 2, std )
    dgauss.SetParLimits( 2, 0.0, 0.5 )
    dgauss.FixParameter( 3, 1.0 )
    dgauss.SetParLimits( 3, 0.5, 1.0 )
    dgauss.FixParameter( 4, 1.0 )
    
    graph = ROOT.TGraphErrors( len(x), x, y / intprod_vec, numpy.zeros( len(err) ), err / intprod_vec )
    fit = graph.Fit( "DoubleGauss", "Q" )
    par_uncorr = [ dgauss.GetParameters()[i] for i in range(5) ]
    dgauss.SetRange( min(x), max(x) )
    
    # Using the integration out to 3 sigma for the integral but leaving the width formula for reference
    #width_uncorr = 1.0 / ( par_uncorr[3] / par_uncorr[2] + (1.0 - par_uncorr[3] ) / par_uncorr[4] / par_uncorr[2] )
    # The maximum of the curve:
    maxi = par_uncorr[0] / math.sqrt(2 * math.pi) / par_uncorr[2] *( par_uncorr[3] + (1-par_uncorr[3])*par_uncorr[4] )
    # The effective width in terms of integral out to 3 sigma and the maximum of the curve
    width_uncorr = dgauss.Integral( -5 * par_uncorr[2] - par_uncorr[1], 5 * par_uncorr[2] - par_uncorr[1] ) / maxi / math.sqrt( 2 * math.pi )
    
    # Time to repeat the fit after performing the beam-beam deflection correction
    sep_corr = beam_beam_deflection( x, intensity, bunches, par_uncorr[1], width_uncorr, beta, isX )
    avg      = numpy.average( sep_corr, weights=y/intprod_vec )
    std      = math.sqrt( numpy.average( (sep_corr-avg)**2, weights=y/intprod_vec) )
    
    graph = ROOT.TGraphErrors( len(sep_corr), sep_corr, y / intprod_vec, numpy.zeros( len(err) ), err / intprod_vec )
    
    # Release the paramters controlling the second gaussian
    dgauss.ReleaseParameter( 3 )
    dgauss.ReleaseParameter( 4 )
    dgauss.SetParLimits( 3, 0.0, 1.0 )
    dgauss.SetParLimits( 4, 0.01, 1.0 )
    dgauss.SetParameter( 4, 1.0  )
    if isEarly and isX:
        dgauss.SetParameter( 4, 1.0 / 1.4  )
    if not isEarly and not isX:
        dgauss.SetParameter( 4, 1.0 / 1.4  )
    # if not isEarly and isX:
    #     dgauss.FixParameter( 3, 1.0  )
    #     dgauss.FixParameter( 4, 1.0  )
    # if  isEarly and not isX:
    #     dgauss.FixParameter( 3, 1.0  )
    #     dgauss.FixParameter( 4, 1.0  )
        
    
    fit = graph.Fit( "DoubleGauss", "QM" )
    par_corr = [ dgauss.GetParameters()[i] for i in range(5) ]
        
    diffs = []
    for i in range( len( sep_corr ) ):
        diffs.append( ( y[i] / intprod_vec[i] - dgauss( sep_corr[i] ) ) / (err[i] / intprod_vec[i]) )
    
    par_corr = dgauss.GetParameters()

    #width_corr = 1.0 / ( par_corr[3] / par_corr[2] + (1.0 - par_corr[3] ) * par_corr[4] / par_corr[2] )
    maxi = par_corr[0] / math.sqrt(2 * math.pi) / par_corr[2] *( par_corr[3] + (1-par_corr[3]) * par_corr[4] )
    width_corr = dgauss.Integral( -5 * width_uncorr - par_corr[1], 5 * width_uncorr - par_corr[1] ) / maxi / math.sqrt( 2 * math.pi )
    #print "X-scan" if isX else "Y-scan", [ par_uncorr[i] for i in range(5) ], "->", [ par_corr[i] for i in range(5) ], dgauss.GetChisquare() / dgauss.GetNDF(), max( sep_corr ) / width_uncorr
    
    # Huge return blob...
    return ( [ par_corr[0], par_corr[1], par_corr[2], par_corr[3], par_corr[4] ], width_corr, sep_corr, 
             [ par_uncorr[0], par_uncorr[1], par_uncorr[2], par_uncorr[3], par_uncorr[4] ], width_uncorr, dgauss.GetChisquare() / dgauss.GetNDF(), diffs )

# Make a map of the train structure
def produceTrainMap( mask, limit):
    # First find the first bunch in the BCID vector that is filled
    # and use as starting point
    first =  list(mask).index(True)
    leading = numpy.zeros( len( mask ) )
    leadingSub = numpy.zeros( len( mask ) )
    leadingStruct = numpy.zeros( len( mask ) )
    leading[first] = 1 if mask[(first+1)%len( mask )] else -1
    
    for i in range( 1, len( mask ) ):
        #If a bunch is empty and the bunch "limit" BCIDs before is
        # empty set the bunch as empty This is to allow to fill the
        # full train structure despite the fact that there are up to 7
        # empty BCIDs before
        if not mask[(i + first)%len( mask )] and not mask[(i + first - limit)%len( mask )]:
            leading[ (i + first)%len( mask ) ] = 0
        else:
            # If the bunch is filled it is in one larger position in the train than the previous one
            leading[ (i + first)%len( mask ) ] = leading[(first + i - 1)%len( mask )] + 1
            # If the bunch is the first in a train, but the one after
            # it is empty then it is isolated (denoted by -1)
            if leading[ (i + first)%len( mask ) ] == 1 and not mask[ (i + first + 1)%len( mask ) ]:
                leading[ (i + first)%len( mask ) ] = -1
        
    #This version is intended to only keep track of subtrains
    for i in range( 1, len( mask ) ):
        if not mask[(i + first)%len( mask ) ]:
            leadingSub[ (i + first)%len( mask ) ] = 0
        else:
            leadingSub[ (i + first)%len( mask ) ] = leadingSub[(first + i - 1)%len( mask )] + 1
            if leadingSub[ (i + first)%len( mask ) ] == 1 and not mask[ (i + first + 1)%len( mask ) ]:
                leadingSub[ (i + first)%len( mask ) ] = -1
        
    # This version encodes how long the train that each bunch belongs to is
    for i in range( 1, len( mask ) ):
        if not mask[i]:
            continue
        else:
            leadingStruct[i] = leading[ i + list(leading[i:]).index(0) - 1]
    for i in range( 1, len( mask ) ):
        if not mask[i]:
            leading[ i ] = 0
    return (leading, leadingSub, leadingStruct )

# All information here is per LB
class raw_data_neo:
    def __init__( self, run, rootfile ):
        self.run   = run
        # Events per LB
        self.evt   = rnp.hist2array( rootfile.Get('Event_ES/h_nEvts_perLB') )
        # Events per LB in nominal physics
        self.evtNP   = rnp.hist2array( rootfile.Get('Event_NP/h_nEvts_perLB') )
        # Duration of each lumiblock
        time  = rootfile.Get('Event_ES/h_LBtime_perLB')
        self.time  = numpy.array( [ time.GetBinContent(i) for i in range( 1, time.GetNbinsX() + 1 ) ] )
        # Acquisition flag (are we meant to take data in the scan or are the beams moving?)
        self.acq   = rnp.hist2array( rootfile.Get('ScanData/h_acqFlag') )
        # Which plane is being scanned? X = 1, Y = 2
        self.plane = rnp.hist2array( rootfile.Get('ScanData/h_nominalSeparationPlane') )
        scans = self.getRanges()
        self.sep   = rnp.hist2array( rootfile.Get('ScanData/h_nominalSeparation') )
        # The code here is taking the first item in each element of scans i.e. apply the function itemgetter( 0 ), which
        # gets the first element, to all elements of scans (that is the map part)
        # Retrieves the lb-range of each scan (self.scans) and the plane being scanned (self.scanPlane)
        self.scans = map( itemgetter( 0 ), scans )
        self.scanPlane = map( itemgetter( 1 ), scans )
        self.nScans = len( self.scans )

        if self.nScans <= 1:
            raise NoScans
        self.bx    = rnp.hist2array( rootfile.Get('ScanData/h_b1DeltaXSet') ) - rnp.hist2array( rootfile.Get('ScanData/h_b2DeltaXSet') )
        self.by    = rnp.hist2array( rootfile.Get('ScanData/h_b1DeltaYSet') ) - rnp.hist2array( rootfile.Get('ScanData/h_b2DeltaYSet') ) 
        # For these histograms we need to keep track of the errors, so
        # just retrieve the histograms directly
        self.trackData = { WP : rnp.hist2array( rootfile.Get('TC_' + WP + '_ES/h_nTrack_perLB') ) for WP  in WPs }
        self.trackDataNP = { WP : rnp.hist2array( rootfile.Get('TC_' + WP + '_NP/h_nTrack_perLB') ) for WP  in WPs }
        self.trackDataDetail = { WP : rootfile.Get('TC_' + WP + '_ES/h_nTrack_perLB_perBCID') for WP  in WPs }
        self.trackDataDetailNP = { WP : rootfile.Get('TC_' + WP + '_NP/h_nTrack_perLB_perBCID') for WP  in WPs }
        # Since the histograms will be deleted once we close the ROOT
        # file make sure the histogram is no longer associated to the
        # file
        for WP in WPs:
            self.trackDataDetail[WP].SetDirectory(0)
            self.trackDataDetailNP[WP].SetDirectory(0)
        self.trackErr =  { WP : rnp.hist2array( rootfile.Get('TC_' + WP + "_ES/h_nTrack_perLB_perEvent" ) ) for WP in WPs }
        self.trackErrNP =  { WP : rnp.hist2array( rootfile.Get('TC_' + WP + "_NP/h_nTrack_perLB_perEvent" ) ) for WP in WPs }
        self.evtPerBCID = rnp.hist2array( rootfile.Get('Event_ES/h_nEvts_perLB_perBCID') )
        self.evtPerBCIDNP = rnp.hist2array( rootfile.Get('Event_NP/h_nEvts_perLB_perBCID') )
        
        self.lbs          = []
        # self.scan_lbs contains a list of lbs for each scan, e.g.
        # [ [ 23, 25, 27, 28], [ 34, 36, 38, 49 ] ]
        self.scan_lbs     = []
        self.fills        = []
        self.beta         = []
        self.date         = []
        self.scanEarly    = []
        self.scanSep      = []
        
        for i in range( len( self.scans ) ):
            #  Extract the incides of the scan
            indices = numpy.where( self.scans[i] )
            # Indices to lbs
            self.scan_lbs.append( indices[0] + 1 )
            # Determine the lb with the maximum tracks per event
            vals = [ ( idx, self.trackData[WPs[0]][idx]/self.evt[idx] ) for idx in indices[0] ]
            lb = max( vals, key=lambda x: x[1] )[0] + 1
            self.lbs.append( lb )
            # Extract the time of the LB for use later
            lbInfo = getRunLbInfo( self.run, lb )[lb]
            # A scan is early if stable beams were not declared two hours before it is performed
            # Time is in nanoseconds, hence 120 min -> 120 * 60e9 ns
            self.scanEarly.append( not getStableIov( int(lbInfo[0]) - 122 * 60e9, int(lbInfo[1]) - 120 * 60e9 ) )
            self.fills.append( getFillNumber( ( lbInfo[0], lbInfo[1] ) ) )
            # Beta star, which is in cm!
            self.beta.append( getBetaStar( ( lbInfo[0], lbInfo[1] ) ) )
            # Date in unix format (IIRC seconds since Jan 1st 1970 normally, in the database it is in nanoseconds again...) 
            self.date.append( lbInfo[0] / 1e9 )
            self.scanSep.append( [ self.bx[point] if self.scanPlane[i] == 1 else self.by[point] for point in indices ][0] )

        self.scanPairs = []
        for i in range( self.nScans - 1 ):
            if any( [ i in x for x in self.scanPairs ] ):
                continue
            # Pair the scans, selection criterion is: one x- and y-scan not separated by more than 30 lumiblocks
            if self.scanPlane[i] + self.scanPlane[i+1] == 3 and self.lbs[i+1] - self.lbs[i] < 30:
                self.scanPairs.append( ( i, i+1 ) )
                continue
        # Select the bunches which actually have some track data
        self.bunchMask = [ self.evtPerBCID[lb] > 0 for lb in self.lbs ]
        first =  list(self.bunchMask[0]).index(True)
        
        # Produce the train structure map
        self.leading = numpy.zeros( len( self.bunchMask[0] ) )
        self.leadingSub = numpy.zeros( len( self.bunchMask[0] ) )
        self.leadingStruct = numpy.zeros( len( self.bunchMask[0] ) )
        self.leading[first] = 1 if self.bunchMask[0][(first+1)%len( self.bunchMask[0] )] else -1
        
        # Please see above: def produceTrainMap( mask, limit):
        # For an explanation
        # maximum difference between trains that is treated as a "subtrain" separation
        limit = 7
        ( self.leading, self.leadingSub, self.leadingStruct ) = produceTrainMap(self.bunchMask[0], limit )
                
        # Correct the intensity of the first bunch 
        # Intensity from FBCT (nominal measurement)
        self.intensity = [ [ getPerBunchIntensity( self.run, lb ) for lb in sorted( self.scan_lbs[i] ) ] for i in range( self.nScans ) ]
        # Intensity from BPTX (for reference of first vs later bunches in the train)
        self.intensityBPTX = [ [ getPerBunchIntensity( self.run, lb, 0 ) for lb in sorted( self.scan_lbs[i] ) ] for i in range( self.nScans ) ]
        for i in range( self.nScans ):
            for j in range( len( self.scan_lbs[i] ) ):
                #Correct the leading bunches ( self.leadingSub == 1 ) to have the same FBCT/BPTX ratio as the later bunches in the train 
                self.intensity[i][j][0][ self.leadingSub == 1 ] = self.intensity[i][j][0][ self.leadingSub == 1 ] *  sum( self.intensityBPTX[i][j][0][ self.leadingSub == 1 ] ) / sum( self.intensity[i][j][0][ self.leadingSub == 1 ] ) * sum( self.intensity[i][j][0][ self.leadingSub  > 1 ] ) / sum( self.intensityBPTX[i][j][0][ self.leadingSub  > 1 ] )
                #Same for beam 2
                self.intensity[i][j][1][ self.leadingSub == 1 ] = self.intensity[i][j][1][ self.leadingSub == 1 ] *  sum( self.intensityBPTX[i][j][1][ self.leadingSub == 1 ] ) / sum( self.intensity[i][j][1][ self.leadingSub == 1 ] ) * sum( self.intensity[i][j][1][ self.leadingSub  > 1 ] ) / sum( self.intensityBPTX[i][j][1][ self.leadingSub  > 1 ] )
        # Determine the correction factor to normalise the FBCT to the DCCT24 readings
        self.corrFact  = [ [ ( getBcidIntIntensity( self.run, self.scan_lbs[i][index] )[0] / numpy.sum(self.intensity[i][index][0]), 
                                  getBcidIntIntensity( self.run, self.scan_lbs[i][index] )[1] / numpy.sum(self.intensity[i][index][1]) ) 
                                for index in range( len( self.scan_lbs[i] ) ) ] for i in range( self.nScans ) ]
        
    def getTrkLumi( self, tmplb, mask, WP=WPs[2] ):
        trk = numpy.array( [ self.trackDataDetail[WP].GetBinContent( self.trackDataDetail[WP].GetBin( int(tmplb), int(i+1) ) ) for i in numpy.where( mask )[0] ] )
        err = numpy.array( [ self.trackDataDetail[WP].GetBinError(   self.trackDataDetail[WP].GetBin( int(tmplb), int(i+1) ) ) for i in numpy.where( mask )[0] ] )
        evt = numpy.array( [ self.evtPerBCID[tmplb-1][i] for i in numpy.where( mask )[0] ] )
        err = err[evt>0]
        trk = trk[evt>0]
        evt = evt[evt>0]
        err_calib = ( numpy.square(err) / evt - numpy.square( trk / evt ) ) / evt
        return ( sum(trk) / sum(evt) / calib[WP], numpy.sqrt( numpy.average( err_calib ) ) / math.sqrt( sum( mask ) ) )

    def getTrkLumiStable( self, tmplb, WP=WPs[2] ):
        mask = self.evtPerBCIDNP[tmplb-1]>0
        
        trk = numpy.array( [ self.trackDataDetailNP[WP].GetBinContent( self.trackDataDetailNP[WP].GetBin( int(tmplb), int(i+1) ) ) for i in numpy.where( mask )[0] ] )
        err = numpy.array( [ self.trackDataDetailNP[WP].GetBinError( self.trackDataDetailNP[WP].GetBin( int(tmplb), int(i+1) ) ) for i in numpy.where( mask )[0] ] )
        evt = numpy.array( [ self.evtPerBCIDNP[tmplb-1][i] for i in numpy.where( mask )[0] ] )
        err_calib = numpy.sqrt( ( numpy.square(err) / evt - numpy.square( trk / evt ) ) / evt )
        return ( sum(trk) / sum(evt) / calib[WP], numpy.sqrt( numpy.average( numpy.square( err_calib ) ) ) / math.sqrt( sum( mask ) ) )

    # Perform the scan analysis
    #
    # Possible to provide a selection function to select which bunches to include
    # useTrk: Select whether LUCID or tracking data should be used
    def analyse( self, excl_fun = lambda x: x, useTrk = False, WP=WPs[2] ):
        new_bunchMask = [ excl_fun( self.bunchMask[i] ) for i in range( len( self.bunchMask ) ) ]
        # Equivalent to:
        # new_bunches = []
        # for i in range( self.nScans ):
        #    new_bunches.append( sum( new_bunchMask[i] ) if sum( new_bunchMask[i] ) > 0 else 1 )
        new_bunches   = [ sum( new_bunchMask[i] ) if sum( new_bunchMask[i] ) > 0 else 1 
                          for i in  range( self.nScans ) ]
        if not useTrk:
            self.lumi_data = [ [ getPerBunchLumi( self.run, lb, new_bunchMask[ i ] ) 
                                 for lb in self.scan_lbs[ i ] ] 
                               for i in range( self.nScans ) ]
        else:
            self.lumi_data = [ [ self.getTrkLumi( lb, new_bunchMask[ i ], WP )
                             for lb in self.scan_lbs[ i ] ] for i in range( self.nScans ) ]
        
        self.int_b1 = [ [ numpy.sum(self.intensity[i][index][0][new_bunchMask[i]]) * self.corrFact[i][index][0] 
                          for index in range( len( self.scan_lbs[i] ) ) ] 
                        for i in range( self.nScans ) ]
        
        self.int_b2 = [ [ numpy.sum(self.intensity[i][index][1][new_bunchMask[i]]) * self.corrFact[i][index][1] 
                          for index in range( len( self.scan_lbs[i] ) ) ] 
                        for i in range( self.nScans ) ]

        results = []
        if len(self.scanPairs):
            scan_num = 0
            # Go through all the scan pairs
            for pair in self.scanPairs:
                output = dict()
                try:
                    # Prepare the x and y data for each scan
                    
                    # Figure out which scan is the x and which is the y one
                    # x_scan = pair[0] if self.scanPlane[pair[0]] == 1 else pair[1]
                    # y_scan = pair[0] if self.scanPlane[pair[0]] == 2 else pair[1]
                    if self.scanPlane[pair[0]] == 1:
                        x_scan = pair[0]
                        y_scan = pair[1]
                    else:
                        x_scan = pair[1]
                        y_scan = pair[0]
                    
                    
                    # Extract the lumi values for the x-scan
                    # each element of self.lumi_data[x_scan] is a pair of lumi, error values for each lumiblock
                    data_x = np.array( [ lumi[0] for lumi in self.lumi_data[x_scan] ] )
                    err_x  = np.array( [ lumi[1] for lumi in self.lumi_data[x_scan] ] )
                    sep_x  = self.scanSep[x_scan]
                    print sep_x
                    print data_x

                    data_y = np.array( [ lumi[0] for lumi in self.lumi_data[y_scan] ] )
                    err_y  = np.array( [ lumi[1] for lumi in self.lumi_data[y_scan] ] )
                    sep_y  = self.scanSep[y_scan]
                      
                    print data_y 
                    print sep_y
                    # Check for zero intensities and if that is the case for any LB abort the whole thing
                    # by raising an exception!
                    if np.any( np.array( self.int_b1[x_scan] ) * np.array( self.int_b2[x_scan] ) == 0 ):
                        print np.array( self.int_b1[x_scan] ) * np.array( self.int_b2[x_scan] )
                        raise "Zero intensity in x-scan!"

                    if np.any( np.array( self.int_b1[y_scan] ) * np.array( self.int_b2[y_scan] ) == 0 ):
                        print np.array( self.int_b1[y_scan] ) * np.array( self.int_b2[y_scan] )
                        raise "Zero intensity in y-scan!"
                    
                    ( curve_x, width_x, sep_corr_x, curve_uncorr_x, width_uncorr_x, chi2_x, diffs_x ) = scanFit( sep_x, data_x, err_x,
                                                                                                                 np.array( self.int_b1[x_scan]) * self.int_b2[x_scan],
                                                                                                                 ( self.int_b1[x_scan][4] + self.int_b2[x_scan][4] ) / 2.0, new_bunches[x_scan], self.beta[pair[0]], True, self.scanEarly[pair[0]] )
                    
                    ( curve_y, width_y, sep_corr_y, curve_uncorr_y, width_uncorr_y, chi2_y, diffs_y ) = scanFit( sep_y, data_y, err_y,
                                                                                                                 np.array( self.int_b1[y_scan]) * self.int_b2[y_scan], 
                                                                                                                 ( self.int_b1[y_scan][4] + self.int_b2[y_scan][4] ) / 2.0, new_bunches[y_scan], self.beta[pair[0]], False, self.scanEarly[pair[0]] )
                    # 2D fits for peak determination
                    peak_corr_spec = ( curve_x[0] + curve_y[0] ) / 2.0
                    exclusion_points = 2
                    data2d = numpy.concatenate( ( (data_x / ( np.array(self.int_b1[x_scan]) * self.int_b2[x_scan] / new_bunches[x_scan]**2 ) )[exclusion_points:-exclusion_points],
                                                  (data_y / ( np.array(self.int_b1[y_scan]) * self.int_b2[y_scan] / new_bunches[y_scan]**2 ) )[exclusion_points:-exclusion_points] ) )
                    err2d  = numpy.concatenate( ( (err_x  / ( np.array(self.int_b1[x_scan]) * self.int_b2[x_scan] / ( new_bunches[x_scan] )**2 ))[exclusion_points:-exclusion_points], 
                                                  (err_y  / ( np.array(self.int_b1[y_scan]) * self.int_b2[y_scan] / ( new_bunches[y_scan] )**2 ))[exclusion_points:-exclusion_points] ) )
                    plot = ROOT.TGraph2DErrors( len( [ sep for sep in sep_corr_x ][exclusion_points:-exclusion_points] + [ 0.0 for sep in sep_corr_y ][exclusion_points:-exclusion_points] ),
                                                np.array([ sep for sep in sep_corr_x ][exclusion_points:-exclusion_points] + [ 0.0 for sep in sep_corr_y ][exclusion_points:-exclusion_points]),
                                                np.array([ 0.0 for sep in sep_corr_x ][exclusion_points:-exclusion_points] + [ sep for sep in sep_corr_y ][exclusion_points:-exclusion_points]),
                                                data2d,
                                                np.zeros( len( [ sep for sep in sep_corr_x ][exclusion_points:-exclusion_points] + [ 0.0 for sep in sep_corr_y ][exclusion_points:-exclusion_points] ) ),
                                                np.zeros( len( [ sep for sep in sep_corr_x ][exclusion_points:-exclusion_points] + [ 0.0 for sep in sep_corr_y ][exclusion_points:-exclusion_points] ) ),                                                          
                                                err2d )
                    gauss2d.SetParameter( 0, peak_corr_spec )
                    gauss2d.SetParLimits( 0, 0., 100. )
                    gauss2d.SetParameter( 1, curve_x[1] )
                    gauss2d.SetParameter( 2, width_x )
                    gauss2d.SetParLimits( 2, 0.0, 0.1 )
                    gauss2d.SetParameter( 3, curve_y[1] )
                    gauss2d.SetParameter( 4, width_y )
                    gauss2d.SetParLimits( 4, 0.0, 0.1 )

                    plot.Fit( gauss2d, "Q" )
                    peak_corr_spec = gauss2d.GetParameters()[0]
                    #width_x = gauss2d.GetParameters()[2]
                    #width_y = gauss2d.GetParameters()[4]
                    
                except Exception as e:
                    print "Error in scan ", scan_num, e
                    traceback.print_exc()
                    raise FitError 

                # The output is put into a dictionary and then appended to the vector of results (one dictionary per scan)
                output['fom']        = peak_corr_spec / exp_lumi( new_bunches[pair[1]], width_x,        width_y )
                output['params']        = ( curve_x, curve_y )
                output['params_uncorr']        = ( curve_uncorr_x, curve_uncorr_y )
                output['fom_uncorr'] = peak_corr_spec / exp_lumi( new_bunches[pair[1]], width_uncorr_x, width_uncorr_y )
                output['mean'] = ( curve_x[1], curve_y[1] )
                output['width'] = ( width_x, width_y )
                output['chi2'] = ( chi2_x, chi2_y, gauss2d.GetChisquare() / gauss2d.GetNDF() )
                output['width_uncorr'] = ( width_uncorr_x, width_uncorr_y )
                output['peak'] = peak_corr_spec * self.int_b1[x_scan][4] * self.int_b2[x_scan][4] / ( new_bunches[y_scan] )**2 
                output['peak_spec'] = peak_corr_spec
                output['peak_exp'] = exp_lumi( new_bunches[pair[1]], width_x, width_y )
                output['curr'] = ( ( self.int_b1[x_scan], self.int_b2[x_scan] ), ( self.int_b1[y_scan], self.int_b2[y_scan] ) )
                output['bunches'] = new_bunches[y_scan]
                output['early'] = self.scanEarly[pair[0]]
                output['date']  = self.date[0]
                output['diffs']  = [ diffs_x, diffs_y ]
                output['run']    = self.run
               #for plot of bell shape 
               #x-scan
                output['sep_x']    = sep_x
                output['err_x']    = err_x
                output['exp_lumi_x']    = data_x
                output['sep_y']    = sep_y
                output['err_y']    = err_y
                output['exp_lumi_y']    = data_y
                results.append( output )
        return results
        print type(results) 
        
        
    # Perform the analysis over the folded bunches
    def analysePerBunchFolded( self, excl_fun = lambda x: x, corr = [0.0, 1.0], useTrk=True ):
        bunches = set( self.leading[ numpy.logical_and( self.leading > 0, self.leadingStruct > 110 ) ] )
        output = {}
        if len( self.leading[ self.leading < 0 ] ) > 0:
            print "++++ Found isolated bunches"
            new_mask = lambda x: np.logical_and( excl_fun(x), self.leading == -1 )
            output[-1] = [ x for x in self.analyse( excl_fun = new_mask, useTrk = useTrk )]
        else:
            print "++++ No isolated bunches"
        
        for i in bunches:
            new_mask = lambda x: np.logical_and( np.logical_and( excl_fun(x), self.leading == i ), self.leadingStruct > 110.0 )
            output[i] = [ x for x in self.analyse( excl_fun = new_mask, useTrk = useTrk )]
        

        # Some plots from the bunch folded that can be useful for debugging the folded output Until the end of the
        # function, the important stuff is already done
        canvas = ROOT.TCanvas( "fom_folded_"+str(self.run), "fom_folded_"+str(self.run), 800, 600 )
        mgr = ROOT.TMultiGraph()
        for fit in range( len( output[1] ) ):
            gr = ROOT.TGraph( len( output ), 
                              np.array( [ i for i in sorted( output.keys() ) ] ), 
                              np.array( [ output[i][fit]['fom'] for i in sorted( output.keys() ) ] ) )
            mgr.Add(gr,"P")
            gr.SetMarkerColor( ROOT.kBlue + fit*3 )
        mgr.Draw("A")
        canvas.Print("neo/foms_folded_"+str(self.run)+".pdf")

        canvas = ROOT.TCanvas( "width_x_folded_"+str(self.run), "width_x_folded_"+str(self.run), 800, 600 )
        gr = ROOT.TGraph( len( output ), 
                          np.array( [ i for i in sorted( output.keys() ) ] ), 
                          np.array( [ output[i][0]['width'][0] for i in sorted( output.keys() ) ] ) )
        gr.Draw("AP")
        canvas.Print("neo/width_x_folded_"+str(self.run)+".pdf")

        canvas = ROOT.TCanvas( "width_y_folded_"+str(self.run), "width_y_folded_"+str(self.run), 800, 600 )
        gr = ROOT.TGraph( len( output ), 
                          np.array( [ i for i in sorted( output.keys() ) ] ), 
                          np.array( [ output[i][0]['width'][1] for i in sorted( output.keys() ) ] ) )
        gr.Draw("AP")
        canvas.Print("neo/width_y_folded_"+str(self.run)+".pdf")

        canvas = ROOT.TCanvas( "mean_y_folded_"+str(self.run), "mean_y_folded_"+str(self.run), 800, 600 )
        gr = ROOT.TGraph( len( output ), 
                          np.array( [ i for i in sorted( output.keys() ) ] ), 
                          np.array( [ output[i][0]['params'][1][1] for i in sorted( output.keys() ) ] ) )
        gr.Draw("AP")
        canvas.Print("neo/mean_y_folded_"+str(self.run)+".pdf")


        canvas = ROOT.TCanvas( "mean_x_folded_"+str(self.run), "mean_x_folded_"+str(self.run), 800, 600 )
        gr = ROOT.TGraph( len( output ), 
                          np.array( [ i for i in sorted( output.keys() ) ] ), 
                          np.array( [ output[i][0]['params'][0][1] for i in sorted( output.keys() ) ] ) )
        gr.Draw("AP")
        canvas.Print("neo/mean_x_folded_"+str(self.run)+".pdf")


        canvas = ROOT.TCanvas( "bunches_folded_"+str(self.run), "bunches_folded_"+str(self.run), 800, 600 )
        gr = ROOT.TGraph( len( output ), 
                          np.array( [ i for i in sorted( output.keys() ) ] ), 
                          np.array( [ output[i][0]['bunches'] for i in sorted( output.keys() ) ] ) )
        gr.Draw("AP")
        canvas.Print("neo/bunches_folded_"+str(self.run)+".pdf")

        canvas = ROOT.TCanvas( "foms_peak_folded_"+str(self.run), "foms_peak_folded_"+str(self.run), 800, 600 )
        mgr = ROOT.TMultiGraph()
        for fit in range( len( output[1] ) ):
            gr = ROOT.TGraph( len( output ), 
                              np.array( [ output[i][fit]['peak'] for i in sorted( output.keys() ) ] ), 
                              np.array( [ output[i][fit]['fom'] for i in sorted( output.keys() ) ] ) )
            mgr.Add( gr, "P" )
            gr.SetMarkerColor( ROOT.kBlue + fit*3 )
        mgr.Draw("A")
        mgr.SetMinimum( 0.95 )
        mgr.SetMaximum( 1.05 )
        canvas.Modified()
        canvas.Print("neo/foms_peak_folded_"+str(self.run)+".pdf")

        canvas = ROOT.TCanvas( "peak_exp_folded_"+str(self.run), "peak_exp_folded_"+str(self.run), 800, 600 )
        gr = ROOT.TGraph( len( output ), 
                          np.array( [ i for i in sorted( output.keys() ) ] ), 
                          np.array( [ output[i][0]['peak_exp'] for i in sorted( output.keys() ) ] ) )
        gr.Draw("AP")
        canvas.Print("neo/peak_exp_folded_"+str(self.run)+".pdf")

        canvas = ROOT.TCanvas( "peak_folded_"+str(self.run), "peak_folded_"+str(self.run), 800, 600 )
        gr = ROOT.TGraph( len( output ), 
                          np.array( [ i for i in sorted( output.keys() ) ] ), 
                          np.array( [ output[i][0]['peak_spec'] for i in sorted( output.keys() ) ] ) )
        gr.Draw("AP")
        canvas.Print("neo/peak_folded_"+str(self.run)+".pdf")

        return output
    
    # Full per bunch analysis
    def analysePerBunchNew( self, excl_fun = lambda x: x, corr = [0.0, 1.0], useTrk=True ):
        output = {}
        
        for i in range(60, 110):
            if self.leading[i] == 0:
                continue
            bunch_set = numpy.zeros( len( self.leading ) )
            bunch_set[i] = True
            new_mask = lambda x: np.logical_and( excl_fun(x), bunch_set )
            output[i] = [ x for x in self.analyse( excl_fun = new_mask, useTrk = useTrk )]
        
        return output

    # Create a list of scanning periods based on the plane and acquiring information
    def getRanges( self ):
        groups = list()
        planes = list()
        for i,j in groupby( enumerate( self.plane ), lambda (k,l): l if l > 0.0 else self.plane[k-1] if ( k > 0 and k < len(self.plane)-1 and self.plane[k-1] == self.plane[k+1] ) else 0.0 ):
            if i > 0.0:
                lbs = map( itemgetter(0), list( j ) )
                groups.append( lbs)
                planes.append( int(i) )
        scans = list()
        for group in zip( groups, planes ):
            mask = np.zeros_like( self.acq, bool)
            np.put( mask, group[0], 1 )
            mask = np.logical_and( mask, self.acq )
            if ( np.sum( mask ) > 7 ):
                scans.append( ( mask, group[1] ) )
        return scans
    
    def __repr__(self):
        return "Run %i: %i scans" % ( self.run, self.nScans )
    def __str__(self):
        return "Run %i: %i scans" % ( self.run, self.nScans )


params = [ 0.001509, 
           1.009506 ]

# Produce the fit data for the file passed as the first argument
#
# e.g. python2 emittance_analysis_neo.py /eos/atlas/atlascerngroupdisk/data-prep/lumi_track_counting/HistFiles/TCLUMI2018_FALL18/data18_13TeV.00360402.calibration_PixelBeam.merge.HIST_IDTRKLUMI.c1278_m2023_c1303_m2042/data18_13TeV.00360402.calibration_PixelBeam.merge.HIST_IDTRKLUMI.c1278_m2023_c1303_m2042._0001.1
def main( filename ):
    run = int( filename.split('/')[-1].split(".")[1] ) 
    print " == Analysing run", run,"=="
    data = raw_data_neo( run, ROOT.TFile(filename) )
    # out_folded = data.analysePerBunchFolded( lambda x : np.array( [ x[i] if i > 180 else 0.0 for i in range(len(x)) ] ), params, useTrk=True )
    # with open("fits/folded/fit_folded_run_"+ str( run ) + ".pickle", 'wb') as fit_pickle:
    #     pickle.dump(out_folded, fit_pickle, protocol=pickle.HIGHEST_PROTOCOL)
    out = data.analyse( excl_fun = lambda x : x, useTrk=False )
    with open("fits/integrated/fit_integrated_lcd_run_"+ str( run ) + ".pickle", 'wb') as fit_pickle:
        pickle.dump(out, fit_pickle, protocol=pickle.HIGHEST_PROTOCOL)
    out = data.analyse( excl_fun = lambda x : x, useTrk=True )
    with open("fits/integrated/fit_integrated_trk_run_"+ str( run ) + ".pickle", 'wb') as fit_pickle:
        pickle.dump(out, fit_pickle, protocol=pickle.HIGHEST_PROTOCOL)

if __name__ == "__main__":
    filename = str(sys.argv[1])
    #print filename
    main( filename )

