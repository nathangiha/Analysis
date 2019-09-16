# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 17:03:08 2019

@author: giha

Generate pulse height/integral spectra and calibrate light output
"""

# Libraries
import numpy as np
import matplotlib
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

# Input list of pulse heights/integrals
def EdgeCal(spec, src = 'cs', histLabel='', xCal=0):
    if src == 'cs':
        edge = csEdge
    if src == 'na':
        edge = naEdge
    
    
    maxrange = 100
    steps = 75
    binedges = np.arange(0,maxrange,maxrange/steps)
    
    specHist, temp = np.histogram(spec, bins = binedges)
        
    # Make histograms
    centers = (binedges[:-1]+binedges[1:])/2
    width = binedges[1]-binedges[0]
    
    

    # Find Compton edge
    deriv = np.diff(MovingAvg(specHist))
    
    zero_crossing = len(centers)-1
    i = 1
    while (specHist[zero_crossing]< max(specHist)/10):
        zero_crossing = np.where(np.diff(np.sign(deriv)))[0][-i]
        i = i +1
        
    # 
    centers_crossing_loc = zero_crossing + 2
    localmax = specHist[centers_crossing_loc]
    
    f = 0.8 # Percentage of max for edgefinding
    
    ind_exceed = (np.where(specHist > localmax*f) )[-1]
    ind_exceed = ind_exceed[-1]
    
    x0=centers[ind_exceed]
    x1=centers[ind_exceed+1]
    y0=specHist[ind_exceed]
    y1=specHist[ind_exceed+1]
    
    edgeLoc = x0 + (f*localmax-y0)*(x1-x0)/(y1-y0)
    
    # Make PI to LO conversion
    conv_factor = edge / edgeLoc


    # Plot everything
    font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 12}

    matplotlib.rc('font', **font)

    #plt.close()
    #plt.figure()
    if xCal ==0:
        plt.bar(centers, specHist, align='center', alpha = 0.5, width = width,
         label = histLabel)
        
        plt.xlabel(r'Pulse Integral $(V\cdot ns)$')
        plt.ylabel(r'Counts')
        p = plt.plot(centers[1:-1], MovingAvg(specHist) )
        plt.axvline(x=edgeLoc, label= histLabel + ' edge @ '+str(round(edgeLoc,2))
        + ' $V\cdot ns$', c = p[-1].get_color())
        plt.legend()
    
    else:
        centers = centers * conv_factor
        plt.bar(centers, specHist, align='center', alpha = 0.5,
         width = width*conv_factor, label = histLabel, edgecolor='k')
        plt.legend()
        
    

    return conv_factor
'''
plt.close('all')
g0 = EdgeCal(csdata[0][:,0], histLabel='G#0')
g1 = EdgeCal(csdata[1][:,0], histLabel='G#1')
plt.savefig('cs_cal.png')
'''
plt.close('all')
S1 = EdgeCal(X2[1][:,0], histLabel='S1', xCal=1)   
G1 = EdgeCal(X6[0][:,0], histLabel='G#1', xCal=1)         
G2 = EdgeCal(X5[0][:,0], histLabel='G#2', xCal=1)
plt.savefig('pixel2_linedup.png')




#EdgeCal(csdata)


# plt.axis([-bound,bound,0,0.200])
# centers, specHist = EdgeCal(csdata[0][:,0])



#plt.plot(centers[1:-1], MovingAvg(specHist) )
#plt.plot(centers[1:-2], deriv )




