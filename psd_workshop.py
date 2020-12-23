# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 15:49:58 2020

@author: giha
"""
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import LogNorm
import numpy as np
from scipy.optimize import leastsq
from scipy.optimize import curve_fit
from scipy.signal import find_peaks, find_peaks_cwt, argrelmax
from scipy.stats import zscore

def _gaussian(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def _2gaussian(x, a1, x01, sigma1,  a2, x02, sigma2):
    return a1*np.exp(-(x-x01)**2/(2*sigma1**2)) +\
           a2*np.exp(-(x-x02)**2/(2*sigma2**2)) 

def DiscLineN(x, a, b,c,):
    return a*np.exp( -b*x) + c
'''
def DiscLineG(x, a, b,c,):
    return -a*np.exp( -b*x) + c
'''
def DiscLineG(x, a, b, c, d):
    return a*x**3 + b*x**2 + c*x + d

plt.close('all')
pi = tagPSD[0]
ratio = tagPSD[1]

pi = psdDataLO[3][0]
ratio = psdDataLO[3][1]
cal = 1
ergbin = 50
minerg = 0
maxerg = 1000
plot = True
    
# plt.close('all')
# Convert pulse integral (or height, I guess) to electron equivalent erg
erg = pi * cal

# Sort erg and ratio vectors based on increaseing erg
sort = np.argsort(erg)
erg = erg[sort]
ratio = ratio[sort]

# Create vector of bin edge ends
ergbin = np.arange(minerg,maxerg+ergbin,ergbin)

# Number of energy bins
ergBinNum=len(ergbin)-1

# Initialize FOM vector
ergvec = np.zeros(ergBinNum)
fomvec = np.zeros(ergBinNum)
locmaxvec = np.zeros((ergBinNum,2))
locmaxamp = np.zeros((ergBinNum,2))
parsmat = np.zeros((ergBinNum,2,3))

# Perform plotting for each energy bin
for i in range(ergBinNum):
    '''
    print( str( ergbin[i] ) )
    '''
    # Define energy bounds
    lowerErg = ergbin[i]
    upperErg = ergbin[i+1]
    
    # Find start and end indices of bins
    lowerInd = np.argmax(erg>lowerErg)
    upperInd = np.argmax(erg>upperErg)
    
    # Create histogram
    binRatio = ratio[lowerInd:upperInd]
    binNum = 1000
    upperRatio = 1
    lowerRatio = 0
    
    
    xSmooth = np.arange(lowerRatio,upperRatio, 1e-5)
    binEdges = np.linspace(lowerRatio,upperRatio, binNum)
    
    ratioHist = np.histogram(binRatio, bins=binEdges)[0]
    
    # Make bar graph stuff
    centers = (binEdges[:-1]+binEdges[1:])/2
    width = binEdges[1]-binEdges[0]
    
    # Fit a double Gaussian to the slice
    
    # Find local maxima
    smoothedRatioHist = np.array(MovingAvg(ratioHist))
    histMax = np.max(smoothedRatioHist)
    # locmax = argrelmax(smoothedRatioHist)
    locmax = find_peaks(smoothedRatioHist, distance = 5, height = histMax / 10)[0]
    
    
    xg = locmax[0]+1
    try:
        xn = locmax[1]+1
    except:
        xn = xg + 10
    
            
    x01 = centers[xg]
    x02 = centers[xn]
    
    a1 = ratioHist[xg]
    a2 = ratioHist[xn]
    
    locmaxvec[i][0] = x01
    locmaxvec[i][1] = x02
    locmaxamp[i][0] = a1
    locmaxamp[i][1] = a2

    sigma1 = 0.0025
    sigma2 = 0.0025
    
    # Fit double gaussian
    popt, pcov = curve_fit(_2gaussian, centers, ratioHist,   \
                            p0 =[a1, x01, sigma1,  a2, x02, sigma2], maxfev=10000)
    
    gauss_2 = _2gaussian(xSmooth, *popt)
    # Separate into 2 gaussians
    pars_g = popt[0:3]
    pars_n = popt[3:6]
    
    if pars_g[0] < pars_n[0]:
        continue
    
    # Add parameters to output
    parsmat[i,0,:] = pars_g
    parsmat[i,1,:] = pars_n
    
    gauss_g = _gaussian(xSmooth, *pars_g)
    gauss_n = _gaussian(xSmooth, *pars_n)
    
    fom = (pars_n[1] - pars_g[1]) / (2.355*(abs(pars_g[2])+abs(pars_n[2])))     
    fomvec[i] = fom
    ergvec[i] = (lowerErg+upperErg)/2
    
    if plot:
        # Plot everything
        plt.figure()

        plt.bar(centers, ratioHist, align='center', alpha = 0.75, width = width,
                label = (str(lowerErg)+'-'+str(upperErg)+' keVee slice') +\
                '\n FOM = ' + str(round(fom,3) ) )
        '''
        plt.bar(centers[2:], smoothedratioHist, align='center', alpha = 0.75, width = width,
                label = (str(lowerErg)+'-'+str(upperErg)+' keVee slice') +\
                '\n FOM = ' + str(round(fom,3) ) )
        '''
        plt.plot(xSmooth, gauss_2, 'k')
        plt.plot(xSmooth, gauss_g, 'r')
        plt.plot(xSmooth, gauss_n, 'b')
        plt.xlabel(r'Tail/Total')
        plt.ylabel(r'Counts')
        plt.legend()
        plt.tight_layout()
        plt.xlim(x01 - 4 * pars_g[2], x02 + 4 * pars_n[2])
        plt.savefig('PSD_'+str(lowerErg)+'_'+str(upperErg), dpi=500)

    
