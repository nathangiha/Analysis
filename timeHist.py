# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 12:12:43 2019

@author: giha

Create time histogram from dataset and evaluate time resolution
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

pFile = 'D:\\X7data.p'

#data = DataLoad(pFile)
def TimeHist():
    ttt = data[1][:,5] - data[0][:,5]
    cfd = data[1][:,4] - data[0][:,4]
    tdif = ttt + cfd
    '''
    for i in range(len(cfd)):
        threshold = 10
        if data[0][i,0] < threshold or data[1][i,0] < threshold:
            tdif[i] = 0
    
    tdif = tdif[np.nonzero(tdif)]
    '''
    bound = 2.5
    # plt.hist(tdif, bins=np.arange(-bound,bound,0.002))
    
    # Plot histograms

    # Histogram time bins
    binedges = np.arange(-bound, bound, step = 0.002)
    
    timeHist, temp = np.histogram(tdif, bins = binedges)
    
    # Normalize to time
    measTime = 50400 # seconds
    timeHist = timeHist / measTime
    # Determine time resolution
    halfmax = max(timeHist)/2
    d = np.sign(halfmax - timeHist[0:-1]) - np.sign(halfmax - timeHist[1:])
    left_idx = np.where(d > 0)[0][0]
    right_idx = np.where(d < 0)[-1][-1]
    
    timeRes = 1000*(binedges[right_idx] - binedges[left_idx])

    
    # Plot histograms
    centers = (binedges[:-1]+binedges[1:])/2
    width = binedges[1]-binedges[0]
    
   # plt.close()
    plt.figure()
    plt.bar(centers, timeHist, align='center', alpha = 0.75, width = width,
     label = 'Data,\nFWHM ='+ str(round(timeRes))+ ' ps')
    
    plt.xlabel(r'$\Delta t$')
    plt.ylabel(r'Count Rate $(s^{-1})$')
    
    plt.axis([-bound,bound,0,max(timeHist)])
    
    #plt.show()
    
    # Fit Gaussian
    mean = sum(centers*timeHist)/len(centers)
    sigma = sum(timeHist*(centers-mean)**2)/len(centers)
    def gaussian(x, a, x0, sigma):
        return a*np.exp(-(x-x0)**2/(2*sigma**2))
    popt, pcov = curve_fit(gaussian, centers, timeHist, p0 = [10000, mean, sigma])
    x = binedges
    y = gaussian(x, *popt)


    plt.plot( x, y, 'r-',
     label='Gaussian fit,\nFWHM ='+ str(round(popt[2]*2350))+ ' ps',
     linewidth = 1)

    '''
    plt.plot( x, y, 'r-',
     label='Gaussian fit, FWHM ='+ str(round(popt[2]*2350))+ 'ps',
     linewidth = 1)
    '''
    
    plt.legend(loc = 2)

    
    
    return popt, left_idx, right_idx, d
plt.close()
j = TimeHist()












