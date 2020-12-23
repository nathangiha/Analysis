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
from scipy.signal import find_peaks



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
def _gaussian(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def lin_interp(x, y, i, half):
    return x[i] + (x[i+1] - x[i]) * ((half - y[i]) / (y[i+1] - y[i]))

def FWHM(x, hist):
    half = np.max(hist)/2.0
    signs = np.sign(np.add(hist, -half))
    zero_crossings = (signs[0:-2] != signs[1:-1])
    zero_crossings_i = np.where(zero_crossings)[0]
    left = lin_interp(x, hist, zero_crossings_i[0], half)
    right = lin_interp(x, hist, zero_crossings_i[1], half)
    return right - left


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


def listEdgeCal(listCal):
    numCals = len(listCal)
    calFactors = np.zeros(numCals)
    
    plt.figure()
    for i in range(numCals):
        calFactors[i] = EdgeCal(listCal[i][0,:], histLabel = 'Cal. '+str(i), xCal=0, integral = True )[0]
        
    return calFactors
    


def fitBackscatterPeak(histogram, continuum, binEdges, unevenFactor = 1, userPeakErg = 0, histLabel = '', 
                       plot = True, percentMaxThreshold = 0.1):

    peakErg = 477.3734081988662
        
    if userPeakErg > 0:
        peakErg = userPeakErg
        
    plotLinewidth = 1
    
    centers = (binEdges[:-1] + binEdges[1:]) / 2

    '''
    normLowerBound = 200
    normUpperBound = 400
  
    binsInBound = np.logical_and(centers > normLowerBound, centers < normUpperBound)
    '''
    
    histSum = 1 # np.sum(histogram[binsInBound])
    contSum = 1 # np.sum(continuum[binsInBound])
    
    histNormed = histogram / histSum
    contNormed = continuum / unevenFactor
    
    
    histErrNormed = np.sqrt(histogram) / histSum
    contErrNormed = np.sqrt(continuum) / contSum
    subtErr = np.sqrt(np.square(histErrNormed) + 
                      np.square(contErrNormed))

    backscatterHistSubt = histNormed - contNormed
    
    
    
    hist = backscatterHistSubt
        
    
    minPeakHeight =  np.max(hist) * percentMaxThreshold
    peaks = find_peaks(hist, height = minPeakHeight)[0]
    
    peakLoc = peaks[-1]
    peakMax = hist[peakLoc]
        # Find FW at Nth max on right, and find point
    N = 100
    fwNthMax = peakMax / N
    peakRightBound = np.argwhere(hist > fwNthMax)[-1][0]
    
    rightSidePeak = hist[peakLoc:peakRightBound]
    leftSidePeak = np.flip(rightSidePeak)[:-1]
    
    artifPeakInd = centers[2 * peakLoc - peakRightBound + 1 : peakRightBound ]
    artifPeak = np.hstack((leftSidePeak, rightSidePeak))
    
    
    popt0, pcov0 = curve_fit(_gaussian,
                           artifPeakInd,
                           artifPeak,
                           p0 = [np.max(artifPeak), centers[peakLoc], 10],
                           maxfev = 10000)
    
    boundSigma = 4
    
    peakRightBound = np.argmax(centers > int(popt0[1] + boundSigma * np.abs(popt0[2])))
    peakLeftBound  = np.argmax(centers > int(popt0[1] - boundSigma * np.abs(popt0[2]))) - 1

    peak = hist[peakLeftBound:peakRightBound]
    peakIndices = centers[peakLeftBound:peakRightBound]

    popt, pcov = curve_fit(_gaussian,
                           peakIndices,
                           peak,
                           p0 = popt0,
                           maxfev = 10000)
    
    
    
    
    fig, ax = plt.subplots(constrained_layout = True,
                           figsize = (6,4))
    plotOrig, = plt.step(centers,
                         histNormed,
                         linewidth = plotLinewidth,
                         c = 'k',
                         where = 'mid',
                         label = 'Gated')
    
    plt.errorbar(centers,
                 histNormed,
                 yerr = histErrNormed,
                 c = 'k',
                 elinewidth = plotLinewidth,
                 capsize = 2,
                 fmt = 'none')


    plotCont, = plt.step(centers,
                         contNormed,
                         linewidth = plotLinewidth,
                         c = 'r',
                         where = 'mid',
                         label = 'Chance coinc.')
    
    plt.errorbar(centers,
                 contNormed,
                 yerr = contErrNormed,
                 c = 'r',
                 elinewidth = plotLinewidth,
                 capsize = 2,
                 fmt = 'none')


    plotSpec, = plt.step(centers,
                         hist,
                         linewidth = plotLinewidth,
                         c = 'b',
                         where = 'mid',
                         label = 'Chance subt.')
    
    plt.errorbar(centers,
                 hist,
                 yerr = subtErr,
                 c = 'b',
                 elinewidth = plotLinewidth,
                 capsize = 2,
                 fmt = 'none')



    fittedPeakLoc = popt[1]
    calibFactor = peakErg / fittedPeakLoc

    try:
        fwhm = FWHM(centers, backscatterHistSubt) * calibFactor
    except:
        fwhm = popt[2] * 2.355
        print('Could not measure FWHM')

    ergRes = 100 * fwhm / peakErg
    
        
        # plt.plot(artifPeakInd, artifPeak, label = 'Mirrored peak')

    xSmooth = np.linspace(peakIndices[0], peakIndices[-1], 1000)
    plotFit, = plt.plot(xSmooth,
                       _gaussian(xSmooth, *popt),
                       'g',
                       linewidth = plotLinewidth * 3,
                       label = 'Gaussian fit,\n' + str(round(ergRes, 2)) + '% energy res.')
    
    # plt.axvline(x = fittedPeakLoc)
    plt.axhline(0,
                c = 'k',
                linewidth = plotLinewidth * 0.5)

    plt.xlabel('Energy (keVee)')
    plt.ylabel('Counts (Norm.) / keVee')
    '''
    plt.title(r'FWHM = ' + str(round(fwhm)) + r' $keVee$, ' + r'Res. = ' + 
              str(round(100 * fwhm / peakErg, 2)) + '%\n' + 
              r'cal = ' + str(round(calibFactor, 4)) + r' $keV (V\cdot ns)^{-1}$')
    '''
    plt.xlim(0, 700)
    
    yLimStep = 100
    # plt.ylim(0, yLimStep * (np.floor(peakMax / yLimStep) + 1))
    
    
    handles,labels = ax.get_legend_handles_labels()

    handles = [handles[2], handles[0], handles[1], handles[3]]
    labels = [labels[2], labels[0], labels[1], handles[3]]
    ax.legend(handles,
              labels,
              loc = 2,
              fontsize = 'medium',
              frameon = False)
    
    
    return calibFactor, ergRes, fwhm, fittedPeakLoc, popt



def peakCal(spec, measTime = 1, src = 'cs', userPeakErg = 0, histLabel = '', xCal = 0, integral = True,
            plot = True, percentMaxThreshold = 0.1, userBinEdges = np.array([])):
    if src == 'cs':
        peakErg = 661.7
    elif src == 'na':
        peakErg = 511
    else:
        peakErg = 100
        
    if userPeakErg > 0:
        peakErg = userPeakErg
        
    plotLinewidth = 1
    
    # Plot initial histogram
    maxRange = np.max(spec)
    steps = 500
    width = maxRange / steps
    binEdges = np.arange(0, maxRange, width)
    specHist, temp = np.histogram(spec, bins = binEdges)
    
    # Adjust bounds based on histogram content
    sumCounts = sum(specHist)
    n = 0   
    while sum(specHist[0:n]) < (sumCounts * 0.999):
        n = n + 1
    
    steps = 200
    maxNew = n * width + 5

    if (userBinEdges.size != 0):
        binEdges = userBinEdges
    else:
        binEdges = np.arange(1, maxNew, maxNew / steps)
        
    centers = (binEdges[:-1] + binEdges[1:]) / 2
    width = binEdges[1] - binEdges[0]
    
    specHist, temp = np.histogram(spec, bins = binEdges)
    
    # Normalize spectrum by measurement time, if specified
    specHist = specHist / measTime
    
    specHistErr = np.sqrt(specHist)
    
    if plot:
        fig, ax = plt.subplots(constrained_layout = True,
                               figsize = (6,4))
        plotSpec, = plt.step(centers,
                            specHist,
                            linewidth = plotLinewidth,
                            c = 'k',
                            where = 'mid',
                            label = 'Backscatter\nspec.')
        
        plt.errorbar(centers,
                     specHist,
                     yerr = specHistErr,
                     c = 'k',
                     elinewidth = plotLinewidth,
                     capsize = 2,
                     fmt = 'none')

    
    # Pick off peak location, taking the right-most peak that exceeds one-tenth
    # of the largest bin
    
    minPeakHeight =  np.max(specHist) * percentMaxThreshold
    peaks = find_peaks(specHist, height = minPeakHeight)[0]
    
    peakLoc = peaks[-1]
    
    # plt.axvline(x = centers[peakLoc])
    
    peakMax = specHist[peakLoc]
    
    # Find FW at Nth max on right, and find point
    N = 100
    fwNthMax = peakMax / N
    peakRightBound = np.argwhere(specHist > fwNthMax)[-1][0]
    
    rightSidePeak = specHist[peakLoc:peakRightBound]
    leftSidePeak = np.flip(rightSidePeak)[:-1]
    
    artifPeakInd = centers[2 * peakLoc - peakRightBound + 1 : peakRightBound ]
    artifPeak = np.hstack((leftSidePeak, rightSidePeak))
    
    
    popt0, pcov0 = curve_fit(_gaussian,
                           artifPeakInd,
                           artifPeak,
                           p0 = [np.max(artifPeak), centers[peakLoc], 10],
                           maxfev = 10000)
    
    boundSigma = 4
    
    peakRightBound = np.argmax(centers > int(popt0[1] + boundSigma * np.abs(popt0[2])))
    peakLeftBound  = np.argmax(centers > int(popt0[1] - boundSigma * np.abs(popt0[2]))) - 1
    
    # print(popt0)
    # print(peakRightBound, peakLeftBound)
    # print(centers[peakRightBound], centers[peakLeftBound])
    '''
    peakToNth = peakRightBound - peakLoc
    peakLeftBound = peakLoc - peakToNth
    '''    
    # Linear interpolate continuum and remove
    peak = specHist[peakLeftBound:peakRightBound]
    continuum = np.linspace(peak[0], peak[-1], len(peak))
    

    peakSubt = peak - continuum
    
    
    specHistFixed = specHist
    specHistFixedErr = specHistErr

    specHistFixed[peakLeftBound:peakRightBound] = peakSubt
    specHistFixedErr[peakLeftBound:peakRightBound] \
        = np.sqrt(peak + continuum)

    peakIndices = centers[peakLeftBound:peakRightBound]

    

    popt1, pcov1 = curve_fit(_gaussian,
                           peakIndices,
                           peakSubt,
                           p0 = popt0,
                           maxfev = 10000)



    fittedPeakLoc = popt1[1]
    calibFactor = peakErg / fittedPeakLoc

    try:
        fwhm = FWHM(peakIndices, peakSubt) * calibFactor
    except:
        fwhm = popt1[2] * 2.355
        print('Could not measure FWHM')

    ergRes = 100 * fwhm / peakErg
    
    if plot:
        # plt.plot(artifPeakInd, artifPeak, label = 'Mirrored peak')
        plt.plot(peakIndices,
                 continuum,
                 'k--',
                 linewidth = 0.5 * plotLinewidth)

        plotContSubt, = plt.step(peakIndices,
                                 peakSubt,
                                 label = 'Continuum\nsubt. peak',
                                 linewidth = plotLinewidth,
                                 c = 'b',
                                 where = 'mid')
        plt.errorbar(peakIndices,
                     peakSubt,
                     yerr = specHistFixedErr[peakLeftBound:peakRightBound],
                     c = 'b',
                     elinewidth = plotLinewidth,
                     capsize = 2,
                     fmt = 'none')


        xSmooth = np.linspace(peakIndices[0], peakIndices[-1], 1000)
        plotFit, = plt.plot(xSmooth,
                           _gaussian(xSmooth, *popt1),
                           'r',
                           linewidth = plotLinewidth * 2,
                           label = 'Gaussian fit,\n' + str(round(ergRes, 2)) + '% energy res.')
        
        # plt.axvline(x = fittedPeakLoc)
    
        plt.xlabel('Energy (keVee)')
        plt.ylabel('Counts / keVee')
        '''
        plt.title(r'FWHM = ' + str(round(fwhm)) + r' $keVee$, ' + r'Res. = ' + 
                  str(round(100 * fwhm / peakErg, 2)) + '%\n' + 
                  r'cal = ' + str(round(calibFactor, 4)) + r' $keV (V\cdot ns)^{-1}$')
        '''
        plt.xlim(0, 700)
        
        yLimStep = 100
        plt.ylim(0, yLimStep * (np.floor(peakMax / yLimStep) + 1))
        
        
        handles,labels = ax.get_legend_handles_labels()

        handles = [handles[2], handles[1], handles[0]]
        labels = [labels[2], labels[1], labels[0]]
        ax.legend(handles,
                  labels,
                  loc = 2,
                  fontsize = 'medium',
                  frameon = False)
        '''
        plt.legend(handles = [plotSpec, plotContSubt, plotFit],
                   loc = 2,
                   fontsize = 'medium',
                   frameon = False)
        '''
    
    return calibFactor, ergRes, fwhm, fittedPeakLoc
    
    
    
    
# Input list of pulse heights/integrals
def EdgeCal(spec, measTime = 1800, src = 'cs', histLabel='', xCal=0, integral = True, edgeFrac=0.62):
    print('Running EdgeCal...')
    
    
    if src == 'cs':
        edge = csEdge
    if src == 'na':
        edge = naEdge
    
    
    maxrange = np.max(spec) # mikwa says mikwassup
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
    maxNew = n*width+5    
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
        # print(str(zero_crossing) + ' ' + str(specHist[zero_crossing]))
        i = i +1
    # 
    centers_crossing_loc = zero_crossing + 2
    localmax = specHist[centers_crossing_loc]
    
    f = edgeFrac # Percentage of max for edgefinding
    
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
    specHist = specHist / localmax

    #plt.close()
    #plt.figure()
    if xCal ==0:
        # plt.bar(centers, specHist, align='center', alpha = 0.5, width = width,
        # label = histLabel)
        plt.bar(centers, specHist, align='center', alpha = 0.5, width = width)
        p = plt.plot(centers[1:-1], MovingAvg(specHist) )

        if integral == True:
            plt.xlabel(r'Pulse Integral $(V\cdot ns)$')
            plt.axvline(x=edgeLoc, label= histLabel + ' edge @ '+str(round(edgeLoc,2))
            + ' $V\cdot ns$', c = p[-1].get_color())
    
        else:
            plt.xlabel(r'Pulse Height $(mV)$')
            plt.axvline(x=edgeLoc, label= histLabel + ' edge @ '+str(round(edgeLoc,2))
            + ' $mV$', c = p[-1].get_color())

        plt.xlim(left=0)
        plt.ylabel(r'Counts (normalized)')
        plt.tight_layout()
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




