# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 12:12:43 2019

@author: giha

Create time histogram from dataset and evaluate time resolution
"""

# Libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit



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

#sys.path.extend([pywaves_directory])

# Load .p datafile to be processed

pFile = 'D:\\X7data.p'



#data = DataLoad(pFile)
def doTimeHist(tdif):
    print('Running TimeHist...')
    '''
    ttt = data[1][:,5] - data[0][:,5]
    cfd = data[1][:,4] - data[0][:,4]
    tdif = ttt + cfd
    '''
    '''
    for i in range(len(cfd)):
        threshold = 10
        if data[0][i,0] < threshold or data[1][i,0] < threshold:
            tdif[i] = 0
    
    tdif = tdif[np.nonzero(tdif)]
    '''
    bound = 2
    # plt.hist(tdif, bins=np.arange(-bound,bound,0.002))
    
    # Plot histograms

    # Histogram time bins
    binedges = np.arange(-bound, bound, step = 0.03)
    
    # Translate to zero
    tdif = tdif - np.average(tdif)
    
    
    
    timeHist, temp = np.histogram(tdif, bins = binedges)
    
    # binedges = binedges - binedges[np.argmax(timeHist)]
    
    centers = (binedges[:-1] + binedges[1:]) / 2
    width = binedges[1] - binedges[0]
    
    fwhm = FWHM(centers, timeHist)
    

    
    # Normalize to time
    #measTime = 1800 # seconds
    #timeHist = timeHist / measTime

    
    # Plot histograms
    
   # plt.close()
   

    plt.figure(figsize = (6, 4))

    plt.hist(tdif,
             bins = binedges,
             histtype = 'step',
             color = 'k',
             label = 'Data,\nFWHM = '+ str(round(fwhm * 1000, 1)) + ' ps',
             zorder = 0)
    
    # plt.bar(centers, timeHist, align='center', alpha = 0.75, width = width,
        #     label = 'Data,\nFWHM = '+ str(round(fwhm * 1000, 2)) + ' ps')
    
    plt.xlabel(r'$\Delta t$ (ns)')
    plt.ylabel(r'Counts/$\Delta t$')
    
    
    yLimStep = 500
    plt.xlim(-bound, bound)
    plt.xticks(np.arange(-bound, bound + 1e-3, step = 1))    


    plt.ylim(0, yLimStep * (np.floor(max(timeHist) / yLimStep) + 1))
    
    #plt.show()
    
    # Fit Gaussian
    a = np.max(timeHist)
    x0 = centers[np.argmax(timeHist)]
    sigma = 0.15 #sum(timeHist*(centers-x0)**2)/len(centers)
    
    def _gaussian(x, a, x0, sigma):
        return a*np.exp(-(x - x0)**2 / (2 * sigma**2))
    
    popt, pcov = curve_fit(_gaussian, centers, timeHist, p0 = [a, x0, sigma])
    xSmooth = np.linspace(centers[0],
                          centers[-1],
                          10000)
    
    gaussFit = _gaussian(xSmooth, *popt)

    fwhmFit = FWHM(xSmooth, gaussFit)

    plt.plot(xSmooth, gaussFit, 'b--',
             label='Gaussian fit,\nFWHM = '+ str(round(fwhmFit * 1000, 1)) + ' ps',
             linewidth = 1,
             zorder = 1)

    '''
    plt.plot( x, y, 'r-',
     label='Gaussian fit, FWHM ='+ str(round(popt[2]*2350))+ 'ps',
     linewidth = 1)
    '''
    
    plt.legend(loc = 1,
               fontsize = 'small',
               frameon = False)
    
    plt.ticklabel_format(axis = 'y',
                         style = 'scientific',
                         scilimits = (3, 4),
                         useOffset = False,
                         useMathText = True)
    
    plt.show()

    # plt.tight_layout()


    
    
    return None
#j = TimeHist()












