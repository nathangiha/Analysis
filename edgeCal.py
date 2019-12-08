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

pFile = 'D:\\X13data.p'

#csdata = DataLoad(pFile)

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
def EdgeCal(spec, measTime = 1800, src = 'cs', histLabel='', xCal=0, integral = True):
    if src == 'cs':
        edge = csEdge
    if src == 'na':
        edge = naEdge
    
    
    maxrange = max(spec)
    steps = 500
    width = maxrange/steps
    binedges = np.arange(0,maxrange,width)    
    specHist, temp = np.histogram(spec, bins = binedges)
    
    # Adjust bounds based on histogram content
    sumCounts = sum(specHist)
    n = 0   
    while sum(specHist[0:n]) < (sumCounts * 0.999):
        n=n+1
    
    steps = 100
    maxNew = n*width    
    binedges = np.arange(0, maxNew, maxNew/steps)
    specHist, temp = np.histogram(spec, bins = binedges)
    
    #specHist = specHist / measTime
    
        
    # Make histograms
    centers = (binedges[:-1]+binedges[1:])/2
    width = binedges[1]-binedges[0]
    
    

    # Find Compton edge
    deriv = np.diff(MovingAvg(specHist))
    
    zero_crossing = len(centers)-1
    i = 1
    while (specHist[zero_crossing]< max(specHist)/10):
        zero_crossing = np.where(np.diff(np.signbit(deriv)))[0][-i]
        print(str(zero_crossing) + ' ' + str(specHist[zero_crossing]))
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
    
    '''
    # Readjust axes based on edge location
    binedgesNew = np.arange(0,edgeLoc*1.5, edgeLoc*1.5/steps)
    specHistNew, temp = np.histogram(spec, bins = binedgesNew)
    centersNew = (binedgesNew[:-1]+binedgesNew[1:])/2
    widthNew = binedgesNew[1]-binedgesNew[0]
    '''
    
    
    
    
    # Make PI to LO conversion
    conv_factor = edge / edgeLoc


    # Plot everything
    font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 12}

    matplotlib.rc('font', **font)
    
    
    # Normalize spectrum to edge max
    #specHist = specHist / localmax

    #plt.close()
    #plt.figure()
    if xCal ==0:
        plt.bar(centers, specHist, align='center', alpha = 0.5, width = width,
         label = histLabel)
        p = plt.plot(centers[1:-1], MovingAvg(specHist) )

        if integral == True:
            plt.xlabel(r'Pulse Integral $(V\cdot ns)$')
            plt.axvline(x=edgeLoc, label= histLabel + ' edge @ '+str(round(edgeLoc,2))
            + ' $V\cdot ns$', c = p[-1].get_color())
    
        else:
            plt.xlabel(r'Pulse Height $(mV)$')
            plt.axvline(x=edgeLoc, label= histLabel + ' edge @ '+str(round(edgeLoc,2))
            + ' $mV$', c = p[-1].get_color())

        plt.ylabel(r'Counts (normalized)')
        plt.legend()
    
    else:
        centers = centers * conv_factor
        plt.bar(centers, specHist, align='center', alpha = 0.5,
         width = width*conv_factor, label = histLabel, edgecolor='k')
        plt.xlabel(r'Light output (keVee)')
        plt.ylabel(r'Counts (normalized)')

        plt.legend()
        
    #plt.xlim(0,maxNew*1.2)
    
    return conv_factor, list(specHist), deriv
'''
plt.close('all')
g0 = EdgeCal(csdata[0][:,0], histLabel='G#0')
g1 = EdgeCal(csdata[1][:,0], histLabel='G#1')
plt.savefig('cs_cal.png')
'''
'''
plt.close('all')
S1 = EdgeCal(X2[1][:,0], histLabel='S1', xCal=0)   
G1 = EdgeCal(X6[0][:,0], histLabel='G#1', xCal=0)         
G2 = EdgeCal(X6[1][:,0], histLabel='G#2', xCal=0)
plt.savefig('pixel2_linedup.png')
'''
'''
plt.close('all')
S1 = EdgeCal(X2[1][:,0], histLabel='S1', xCal=0)   
G1 = EdgeCal(X6[0][:,0], histLabel='G#1', xCal=0)         
G2 = EdgeCal(X5[0][:,0], histLabel='G#2', xCal=0)
plt.savefig('pixel2.png')
'''
'''
plt.close('all')
S1 = EdgeCal(X2[1][:,0], histLabel='S1', xCal=0)   
G1 = EdgeCal(X3[0][:,0], histLabel='G#1.1', xCal=0)         
G1 = EdgeCal(X6[0][:,0], histLabel='G#1.2', xCal=0)         
'''
'''
plt.close('all')
G1 = EdgeCal(X6[0][:,0], histLabel='G#1 pre', xCal=0)         
G2 = EdgeCal(X6[1][:,0], histLabel='G#2 pre', xCal=0)
G1 = EdgeCal(X8[0][:,0], histLabel='G#2 post', xCal=0)         
G2 = EdgeCal(X8[1][:,0], histLabel='G#2 post', xCal=0)
'''
#EdgeCal(csdata)


# plt.axis([-bound,bound,0,0.200])
# centers, specHist = EdgeCal(csdata[0][:,0])



#plt.plot(centers[1:-1], MovingAvg(specHist) )
#plt.plot(centers[1:-2], deriv )




