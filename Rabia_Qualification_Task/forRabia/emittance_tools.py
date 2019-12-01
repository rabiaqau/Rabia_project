from itertools import groupby
from operator import itemgetter
from scipy.optimize import curve_fit
from struct import unpack
import math
import numpy
import numpy as np
import os
import sys

# Create the database access tools
from PyCool import cool
dbSvc=cool.DatabaseSvcFactory.databaseService()
dbconnTrig    = dbSvc.openDatabase("COOLONL_TRIGGER/CONDBR2",True)
dbconnTdaq    = dbSvc.openDatabase("COOLONL_TDAQ/CONDBR2"   ,True)
dbconnTdaqMon = dbSvc.openDatabase("COOLONL_TDAQ/MONP200"   ,True)
dbconnTrigOff = dbSvc.openDatabase("COOLOFL_DCS/CONDBR2"    ,True)

# Once object for each folder in the database
lbfolderTrig=dbconnTrig.getFolder('/TRIGGER/LUMI/LBLB')
dbfolderLumi=dbconnTrig.getFolder('/TRIGGER/LUMI/OnlPrefLumi')
dbfolderTdaq=dbconnTdaq.getFolder('/TDAQ/OLC/LHC/LBDATA')
dbfolderTdaqBunch=dbconnTdaq.getFolder('/TDAQ/OLC/LHC/BUNCHDATA')
dbfolderTdaqBunchLumi=dbconnTdaqMon.getFolder('/TDAQ/OLC/BUNCHLUMIS')
dbfolderTdaqBunchLumiCalib=dbconnTdaq.getFolder('/TDAQ/OLC/CALIBRATIONS')
dbAccounting = dbconnTrigOff.getFolder('/LHC/DCS/FILLSTATE')

import cPickle as pickle

def writeCond():
     outfile = open(cond_pickle,'wb')
     pickle.dump(cond_dict,outfile)
     outfile.close()

cond_pickle = "condPickle.pkl"
cond_dict = dict()
if os.path.isfile(cond_pickle):
     infile = open(cond_pickle,'rb')
     cond_dict = pickle.load(infile)
     infile.close()
else:
     cond_dict["lumi"] = dict()
     cond_dict["int"]  = dict()
     writeCond()


# Import all ghost charge information, in terms of fill number
import csv
ghost_charge = {}
with open('ghost_charge.txt', 'rb') as csvfile:
     spamreader = csv.reader(csvfile, delimiter=' ', quotechar='#')
     for row in spamreader:
         ghost_charge[ int(row[0]) ] = float(row[2])


def getFillNumber( run_iov ):

     citr = dbAccounting.browseObjects( int(run_iov[0]), int(run_iov[1]), cool.ChannelSelection.all() )
     fill = 0
     for cobj in citr:
          fill = cobj.payload()['FillNumber']

     return int( fill )

def getBetaStar( run_iov ):

     citr = dbAccounting.browseObjects( int(run_iov[0]), int(run_iov[1]), cool.ChannelSelection.all() )
     fill = 0
     for cobj in citr:
          fill = cobj.payload()['BetaStar']

     return int( fill )

def getFillDate( run_iov ):

     citr = dbAccounting.browseObjects( int(run_iov[0]), int(run_iov[1]), cool.ChannelSelection.all() )
     fill = 0
     for cobj in citr:
          fill = cobj.payload()['IOV']

     return int( fill )

def getRunLbInfo( run, lb ):
    iovstart = ( run << 32 ) + lb 
    iovend =   ( run << 32 ) + lb 
    citr=lbfolderTrig.browseObjects(iovstart,iovend,cool.ChannelSelection.all())
    
    lbs =  { int( cobj.since() & 0xFFFFFFFF): ( cobj.payload()['StartTime'], cobj.payload()['EndTime']  ) for cobj in citr }

    return lbs


def getLumi( run, lb ):
    iovstart = ( run << 32 ) + lb
    iovend =   ( run << 32 ) + lb

    citr=dbfolderLumi.browseObjects(iovstart,iovend,cool.ChannelSelection(0))
    
    data = dict()
    for cobj in citr:
        data[ int( cobj.since() & 0xFFFFFFF ) ] = cobj.payload()['LBAvInstLumi']
    

    return data[lb]

# Channel 147 is LUCID C12, hardcoded since for this year it is the main
# luminometer to use!
def getLumiCalib( run, lb, channel=147 ):
    lbInfo = getRunLbInfo( run, lb )
    ( iovstart, iovend ) = ( lbInfo[lb][0], lbInfo[lb][1] )

    citr=dbfolderTdaqBunchLumiCalib.browseObjects(iovstart,iovend,cool.ChannelSelection(channel))
    
    data = dict()
    for cobj in citr:
        blob = cobj.payload()['Parameters']
        data[ lb ] = (cobj.payload()['Function'], [ unpack( 'f', blob.read(4))[0] for i in range( cobj.payload()['NumOfParameters'] ) ] )

    return  data[lb]

# Get the per bunch luminosity
#
# Channel 147 is LUCID C12, hardcoded since for this year it is the main
# luminometer to use!    
def getPerBunchLumi( run, lb, mask ):
     channel = 147

     if ( run, lb ) not in cond_dict["lumi"].keys():
          lbInfo = getRunLbInfo( run, lb )
          ( iovstart, iovend ) = ( lbInfo[lb][0], lbInfo[lb][1] )
          
          citr = dbfolderTdaqBunchLumi.browseObjects( iovstart, iovend, cool.ChannelSelection(channel) )
          
          for cobj in citr:
               if int( cobj.payload()['RunLB'] & 0xFFFFFFF ) != lb:
                    continue
               blob = cobj.payload()['BunchRawInstLum']
               if len(blob)!=1+3564*4:
                    print "BAD info in COOL!"
                    return 0
               head1 =  unpack( 'B', blob.read(1) )[0]
               cond_dict["lumi"][( run,lb )] = np.array( [ unpack( 'f', blob.read(4))[0] for i in range(0, 3564) ] )
     
     data = cond_dict["lumi"][( run,lb )][mask]
     
     if np.any( data > 1.0 ):
          print "====== Out of range data in run", run, " and lb", lb
     calibdata = - np.log( 1 - data ) *                                               80000. / 6890.4
     caliberr  =   1.0 / ( 1 - data ) * numpy.sqrt( data * ( 1 - data ) / 100000. ) * 80000. / 6890.4
     
     calibdata = np.average( calibdata )
     caliberr  = np.sqrt( np.average( np.square( caliberr ) ) )
     
     caliberr  = 1.009506 * caliberr  + 2 * calibdata * caliberr * 0.001509
     calibdata = calibdata * 1.009506 - calibdata * calibdata * 0.001509
     
     return ( calibdata, caliberr )

# Get per bunch luminosity without integrating over the bunches
def getPerBunchLumiAll( run, lb, mask ):
     channel = 147

     if ( run, lb ) not in cond_dict["lumi"].keys():
          lbInfo = getRunLbInfo( run, lb )
          ( iovstart, iovend ) = ( lbInfo[lb][0], lbInfo[lb][1] )
          
          citr = dbfolderTdaqBunchLumi.browseObjects( iovstart, iovend, cool.ChannelSelection(channel) )
          
          for cobj in citr:
               if int( cobj.payload()['RunLB'] & 0xFFFFFFF ) != lb:
                    continue
               blob = cobj.payload()['BunchRawInstLum']
               if len(blob)!=1+3564*4:
                    print "BAD info in COOL!"
                    return 0
               head1 =  unpack( 'B', blob.read(1) )[0]
               cond_dict["lumi"][( run,lb )] = np.array( [ unpack( 'f', blob.read(4))[0] for i in range(0, 3564) ] )
     
     data = cond_dict["lumi"][( run,lb )][mask]
     if np.any( data > 1.0 ):
          print "====== Out of range data in run", run, " and lb", lb
     calibdata = - np.log( 1 - data ) *                                               80000. / 6890.4
     caliberr  =   1.0 / ( 1 - data ) * numpy.sqrt( data * ( 1 - data ) / 100000. ) * 80000. / 6890.4
     
     return ( calibdata, caliberr )

def getStable( run, lb ):
    lbInfo = getRunLbInfo( run, lb )
    run_iov = ( lbInfo[lb][0], lbInfo[lb][1] )


    citr = dbAccounting.browseObjects( int(run_iov[0]), int(run_iov[1]), cool.ChannelSelection.all() )
    stable = False
    for cobj in citr:
         stable = cobj.payload()['StableBeams']

    return stable

def getStableIov( iov_start, iov_end ):
    citr = dbAccounting.browseObjects( int(iov_start), int(iov_end), cool.ChannelSelection.all() )
    stable = False
    for cobj in citr:
         stable = cobj.payload()['StableBeams']
         return stable
    
    return stable

# Get the ghost fraction for the particular fill
# This is unbunched charge that is not contributing to luminosity,
# but contributes to the current measurements
def getGhostFrac( run, lb ):
    lbInfo = getRunLbInfo( run, lb )
    run_iov = ( lbInfo[lb][0], lbInfo[lb][1] )

    fill = getFillNumber( run_iov ), 
    return ghost_charge[ int( fill[0] ) ] if int( fill[0] ) in ghost_charge.keys() else 0.0 


def getBcidIntIntensity( run, lb ):
    lbInfo = getRunLbInfo( run, lb )
    run_iov = ( lbInfo[lb][0], lbInfo[lb][1] )

    citr = dbfolderTdaq.browseObjects( int(run_iov[0]), int(run_iov[1]), cool.ChannelSelection(3) )
    data = { int( cobj.payload()['RunLB'] & 0xFFFFFFFF) : ( cobj.payload()['Beam1IntensityAll'] / 1e11, cobj.payload()['Beam2IntensityAll'] / 1e11 ) for cobj in citr }

    return data[lb]

def getPerBunchIntensity( run, lb, channel = 1 ):
    if ( run, lb ) in cond_dict["int"].keys():
         return cond_dict["int"][(run,lb)]
    
    lbInfo = getRunLbInfo( run, lb )
    run_iov = ( lbInfo[lb][0], lbInfo[lb][1] )

    citr = dbfolderTdaqBunch.browseObjects( int(run_iov[0]), int(run_iov[1]), cool.ChannelSelection( channel ) )
    
    data = dict()
    for cobj in citr:
        run = int( cobj.payload()['RunLB'] & 0xFFFFFFFF)
        blob1, blob2 = cobj.payload()['B1BunchIntensities'], cobj.payload()['B2BunchIntensities']
        head1, head2 =  unpack( 'B', blob1.read(1) )[0], unpack( 'B', blob2.read(1))[0]
        data[run] = ( np.array( [ unpack( 'f', blob1.read(4))[0] for i in range(0, 3564) ] ),
                      np.array( [ unpack( 'f', blob2.read(4))[0] for i in range(0, 3564) ] ) )
    cond_dict["int"][(run,lb)] = data[lb]
    
    return data[lb]

# Parameter 0: amplitude
# Parameter 1: mean (common between two gaussian)
# Parameter 2: width of first gaussian
# Parameter 3: Split of amplitude between two gaussians range is 0-1)
# Parameter 4: Factor by which second gaussian is wider than first, i.e. sigma2 = sigma1 * parameter 4 p4 is at least 1
# MAX is: p[0] / math.sqrt(2 * math.pi) / p[2] *( p[3] + (1-p[3])/p4 )
def gaussD(x, p):
     return p[0] / math.sqrt(2 * math.pi) / p[2] * ( p[3]  * numpy.exp( - ( (x[0]-p[1])**2 / (2.*p[2]**2) ) ) + ( 1.0 - p[3] ) / p[4] *numpy.exp( - ( (x[0]-p[1])**2/(2.*(p[2]*p[4])**2) ) ) )

def gauss2D(x, p):
     return p[0] * ( numpy.exp( - (x[0]-p[1])**2/(2.*p[2]**2) ) * numpy.exp( - (x[1]-p[3])**2/(2.*p[4]**2)  ) )

def gauss(x, *p):
     return p[0] * numpy.exp( - ( (x-p[1])**2 / (2.*p[2]**2) ) )

def dualGauss(x, *p):
     X, Y = x 
     A, mu_x, sigma_x, mu_y, sigma_y = p
     return A*numpy.exp( - ( (X-mu_x)**2/(2.*sigma_x**2) )  ) * numpy.exp( - ( (Y-mu_y)**2/(2.*sigma_y**2) )  )
