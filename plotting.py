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
from functions import*

def workingpoint(WP):
    WP[0]['early_peak']


pickle_in_tightmod = open("WP_output/tightmod.pickle","rb")
scan_tightmod = pickle.load(pickle_in_tightmod)
print scan_tightmod[0]['early_peak'],"scan"
