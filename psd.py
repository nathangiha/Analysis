# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 14:31:18 2019

@author: giha

Do PSD
"""
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import LogNorm
import numpy as np
import heapq
from scipy.optimize import leastsq
from scipy.optimize import curve_fit
from scipy.signal import argrelextrema

#data = X14

#ph = data[0][:,0]
#conv_factor = 0.18835

data = B1cf

intg = data[0,:]
ratio = data[3,:]
conv_factor = 6.823

# Convert to LO, electron energy equivalent
intg = intg * conv_factor

# Remove pulses above energy limit
ergLimit = 1500
clipped = np.where(intg < ergLimit)[0]


#tail = data[0][clipped,2]
#total = data[0][clipped,3]


intg = intg[clipped]
ratio = ratio[clipped]
del(clipped)
#del(tail, total)
test = 0

def PSD_hist(binnum=300):
    plt.close()
#    psd = plt.hist2d(ph,ratio, bins=(binnum,50*binnum ), cmap=plt.cm.jet,
 #                    norm=LogNorm())
    psd = plt.hist2d(intg,ratio, bins=(binnum,binnum ), cmap=plt.cm.jet
                    , norm=LogNorm())
    plt.colorbar()
    plt.xlim(0, ergLimit)
    plt.ylim(0.8,1) 
    plt.xlabel(r'Pulse Integral (keVee)')
    plt.ylabel(r'Tail/Total')

def FOM_plot(pi, ratio):
    fomvecs = PSD_ergslice(pi,1,ratio, plot=0)
    erg = fomvecs[0]
    fom = fomvecs[1]
    
    ergstep = erg[1]-erg[0]
    
    plt.figure()
    plt.plot(erg,fom, 'b-')
    plt.scatter(erg,fom, c='b')
    plt.xlabel('Energy $(keVee)$')
    #plt.xlabel(r'Energy ($' + str(np.round(ergstep)) +' keVee$ bins)')
    plt.ylabel(r'$\frac{\delta}{FWHM_\gamma + FWHM_n} $')
    plt.legend()
    plt.tight_layout()
    i=3
    plt.annotate(str(round(erg[i]-ergstep/2,0)) + '-' +\
                 str(round(erg[i]+ergstep/2,0))+ ' keVee, ' \
                 + str(round(fom[i],2)), (erg[i], fom[i]))
    #for i in range(len(fom)):
    #    plt.annotate(str(round(erg[i],0)) + ', ' + str(round(fom[i],2)), (erg[i], fom[i]))

    

def PSD_ergslice( pi, cal, ratio, ergbin=50, maxerg = 1000, plot=1 ):
    
    plt.close('all')
    # Convert pulse integral (or height, I guess) to electron equivalent erg
    erg = pi * cal
    
    # Sort erg and ratio vectors based on increaseing erg
    sort = np.argsort(erg)
    erg = erg[sort]
    ratio = ratio[sort]
    
    # Create vector of bin edge ends
    ergbin = np.arange(0,1000+ergbin,ergbin)
    
    # Number of energy bins
    ergBinNum=len(ergbin)-1
    
    # Initialize FOM vector
    ergvec = np.zeros(ergBinNum)
    fomvec = np.zeros(ergBinNum)
    
    # Perform plotting for each energy bin
    for i in range(ergBinNum):
        
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
        binEdges = np.arange(lowerRatio,upperRatio, (upperRatio-lowerRatio)/binNum)
        
        ratioHist = np.histogram(binRatio, bins=binEdges)[0]
        
        # Make bar graph stuff
        centers = (binEdges[:-1]+binEdges[1:])/2
        width = binEdges[1]-binEdges[0]
        
        # Fit a double Gaussian to the slice
        
        # Find local maxima
        locmax = argrelextrema(ratioHist, np.greater)       
        
        xg = locmax[0][0]
        xn = locmax[0][1]
               
        x01 = centers[xg]
        x02 = centers[xn]
        a1 = ratioHist[xg]
        a2 = ratioHist[xn]
        
        print(str(i) + ' ' + str(x01) + ' ' +str(x02))
        
        
        '''
 
        # Guesses
        a1 = max(ratioHist)
        a2 = a1 /11
        x01 = centers[np.argmax(ratioHist)]
        x02 = 0.8891
        '''
        sigma1 = 0.01
        sigma2 = 0.01
        
        # Fit double gaussian
        popt, pcov = curve_fit(_2gaussian, centers, ratioHist,   \
                                p0 =[a1, x01, sigma1,  a2, x02, sigma2])
        
        gauss_2 = _2gaussian(xSmooth, *popt)
        # Separate into 2 gaussians
        pars_g = popt[0:3]
        pars_n = popt[3:6]
        
        gauss_g = _gaussian(xSmooth, *pars_g)
        gauss_n = _gaussian(xSmooth, *pars_n)
        
        fom = (pars_n[1] - pars_g[1]) / (2.355*(abs(pars_g[2])+abs(pars_n[2])))     
        fomvec[i] = fom
        ergvec[i] = (lowerErg+upperErg)/2
        
        if plot != 0:
            # Plot everything
            plt.figure()

            plt.bar(centers, ratioHist, align='center', alpha = 0.75, width = width,
                    label = (str(lowerErg)+'-'+str(upperErg)+' keVee slice') +\
                    '\n FOM = ' + str(round(fom,3) ) )
            
            plt.plot(xSmooth, gauss_2, 'k')
            plt.plot(xSmooth, gauss_g, 'r')
            plt.plot(xSmooth, gauss_n, 'b')
            plt.xlabel(r'Tail/Total')
            plt.ylabel(r'Counts')
            plt.legend()
            plt.tight_layout()
            plt.xlim(0,1)
            plt.savefig('PSD_'+str(lowerErg)+'_'+str(upperErg), dpi=500)

    return ergvec, fomvec
        
def _gaussian(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))
    
def _2gaussian(x, a1, x01, sigma1,  a2, x02, sigma2):
    return a1*np.exp(-(x-x01)**2/(2*sigma1**2)) +\
           a2*np.exp(-(x-x02)**2/(2*sigma2**2)) 
    
    
    
    
    
    
    
    