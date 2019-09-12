# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 17:03:08 2019

@author: giha

Generate pulse height/integral spectra and calibrate light output
"""

# Libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.optimize import leastsq
from scipy.optimize import curve_fit



#sys.path.extend([pywaves_directory])

# Load .p datafile to be processed
from loader import DataLoad

pFile = 'D:\\X6data.p'

csdata = DataLoad(pFile)

### Editables ###


#################


#################
csEdge = 478 # keVee
naEdge = 341 # keVee

#################

# Input list of pulse heights/integrals
def EdgeCal(spec, src = 'cs'):
    if src == 'cs':
        edge = csEdge
    
    binedges = np.arange(0,100,1)
        # Histogram time bins
    
    binedges = np.arange(0,100,1)
    
    specHist, temp = np.histogram(spec, bins = binedges)
        
    # Plot histograms
    centers = (binedges[:-1]+binedges[1:])/2
    width = binedges[1]-binedges[0]
    
    plt.close()
    plt.bar(centers, specHist, align='center', alpha = 0.75, width = width, label = 'glass')
    
    plt.xlabel(r'Pulse Integral $(V\cdot ns)$')
    plt.ylabel(r'Counts')
    
    # plt.axis([-bound,bound,0,0.200])

    
    
    
    
    return centers, specHist

def MovingAvg(spec):
    N = 3
    cumsum, moving_aves = [0], []
    for i, x in enumerate(spec, 1):
        cumsum.append(cumsum[i-1] + x)
        if i>=N:
            moving_ave = (cumsum[i] - cumsum[i-N])/N
            #can do stuff with moving_ave here
            moving_aves.append(moving_ave)
    return moving_aves



centers, specHist = EdgeCal(csdata[0][:,0])

deriv = np.diff(MovingAvg(specHist))

zero_crossing = np.where(np.diff(np.sign(deriv)))[0][-1]
centers_crossing_loc = zero_crossing + 2
localmax = specHist[centers_crossing_loc]





plt.plot(centers[1:-1], MovingAvg(specHist) )
plt.plot(centers[1:-2], deriv )




