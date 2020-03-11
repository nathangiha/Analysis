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
from scipy.optimize import leastsq
from scipy.optimize import curve_fit
from scipy.signal import argrelmax

#data = X14

#ph = data[0][:,0]
#conv_factor = 0.18835

#data = B1cf



def CalNClip(intg, ratio, conv_factor, ergLimit = 1500):
    # Convert to LO
    intconv = intg* conv_factor


    # Remove pulses above energy limit
    clipped = np.where(intconv < ergLimit)[0]


    intconv = intconv[clipped]
    ratio = ratio[clipped]
    
       
    del(clipped)
    return [intconv, ratio]

# Moving average
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

def PSD_hist(pi,ratio,binnum=300,ergLimit=1500, discLine=0, binsize=50):
    plt.figure()
    

    
    
    intp = pi[(ratio >=0) & (ratio <=1) ]
    ratiop = ratio[(ratio >=0) & (ratio <=1) ]
    
#    psd = plt.hist2d(ph,ratio, bins=(binnum,50*binnum ), cmap=plt.cm.jet,
 #                    norm=LogNorm())
    psd = plt.hist2d(intp,ratiop, bins=(binnum,binnum*10 ),\
        range = [[0,ergLimit],[0,1]], cmap=plt.cm.jet, norm=LogNorm())
    
    maxInd = np.where(psd[0] == np.amax(psd[0]))
    distLoc = psd[2][maxInd[1][0]]
    plt.clf()
    plt.xlim(0, ergLimit)
    plt.ylim(distLoc-0.05,distLoc+0.05)


    # If we are plotting a discrimination line (where the value of discLine is
    # the number of standard deviations away from the gamma Gaussian)
    if discLine >0:
        slices = PSD_ergslice(pi,ratio, ergbin = binsize, maxerg = ergLimit, plot=0)
        erg = slices[0]
        parsmat = slices[2]
        
        discPoints = parsmat[:,0,1]+ discLine * parsmat[:,0,2]
        #plt.autoscale(False)
        #plt.scatter(erg, discPoints)
        #plt.scatter(erg, discPoints, zorder = 1)
        popt, pcov = curve_fit(DiscLine, erg, discPoints,   \
                                p0 =[1, 0.1,1])
        
        eSmooth = np.arange(0,ergLimit,1)
        discFit = DiscLine(eSmooth, *popt)
        
        plt.plot(eSmooth, discFit, 'k', label = str(discLine)+ r'$\sigma$')

    
    
    
    psd = plt.hist2d(intp,ratiop, bins=(binnum,binnum ),\
        range = [[0,ergLimit],[distLoc-0.05,distLoc+0.05]], cmap=plt.cm.jet,\
        zorder = 0, norm=LogNorm())
    plt.colorbar()
    plt.xlabel(r'Pulse Integral (keVee)')
    plt.ylabel(r'Tail/Total')
    plt.legend()
    plt.tight_layout()
 

def FOM_plot(pi, ratio, binsize=50, maxerg= 1000):
    fomvecs = PSD_ergslice(pi,ratio, ergbin = binsize, maxerg = maxerg, plot=0)
    erg = fomvecs[0]
    fom = fomvecs[1]
    
    # Neutralize outliers
    fom = fom[ (fom >= 5) or (fom <= 0) ]    
    erg = erg[ (fom > 5) or (fom < 0) ]  
    ergstep = erg[1]-erg[0]
    
    plt.figure()
    plt.plot(erg,fom, 'b-')
    plt.scatter(erg,fom, c='b')
    plt.xlabel('Energy $(keVee)$')
    #plt.xlabel(r'Energy ($' + str(np.round(ergstep)) +' keVee$ bins)')
    plt.ylabel(r'$\frac{\delta}{FWHM_\gamma + FWHM_n} $')
    
    plt.ylim(bottom=0)

    
    plt.legend()
    plt.tight_layout()
    i=3
    plt.annotate(str(round(erg[i]-ergstep/2,0)) + '-' +\
                 str(round(erg[i]+ergstep/2,0))+ ' keVee, ' \
                 + str(round(fom[i],2)), (erg[i], fom[i]))
    #for i in range(len(fom)):
    #    plt.annotate(str(round(erg[i],0)) + ', ' + str(round(fom[i],2)), (erg[i], fom[i]))

    

def PSD_ergslice( pi, ratio, ergbin=50, maxerg = 1000, cal=1, plot=1 ):
    
    plt.close('all')
    # Convert pulse integral (or height, I guess) to electron equivalent erg
    erg = pi * cal
    
    # Sort erg and ratio vectors based on increaseing erg
    sort = np.argsort(erg)
    erg = erg[sort]
    ratio = ratio[sort]
    
    # Create vector of bin edge ends
    ergbin = np.arange(0,maxerg+ergbin,ergbin)
    
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
        smoothedratioHist = np.array(MovingAvg(ratioHist))
        locmax = argrelmax(smoothedratioHist)
        '''
        xg = np.argmax(smoothedratioHist)+3
        
        locmaxvec[i] = locmax[0][0]
        
        try:
            xn = locmax[0][np.argmin(  np.abs(locmax[0] - xg ) ) +1]+3
        except:
            xn = xg+13
        '''
        
        
        
        xg = locmax[0][0]+1
        j = 0
        
        # While loop acts as high-pass filter in histogram, filtering out
        # low-amplitude noise before choosing maxima
        
        
        while (ratioHist[xg]< 15):
            xg = locmax[0][j]+1
            j = j+1
        try:
            xn = locmax[0][j+1]+1
        except:
            xn = xg+20
            
            
        # Band-aid solution to the argrelmax function sucking sometimes
        if xn > xg+50:
            xn = xg+20
    
        '''    
            
        xg = locmax[0][0]+1
        try:
            xn = locmax[0][1]+1
        except:
            xn = xg+20
                
        '''        
                
        x01 = centers[xg]
        x02 = centers[xn]
        
        a1 = ratioHist[xg]
        a2 = ratioHist[xn]
        
        locmaxvec[i][0] = x01
        locmaxvec[i][1] = x02
        locmaxamp[i][0] = a1
        locmaxamp[i][1] = a2

        print(str(upperErg/100) + ' ' + str(x01) + ' ' +str(x02))
        
        print(str(locmax[0]))
        '''
 
        # Guesses
        a1 = max(ratioHist)
        a2 = a1 /11
        x01 = centers[np.argmax(ratioHist)]
        x02 = 0.8891
        '''
        sigma1 = 0.0025
        sigma2 = 0.0025
        
        # Fit double gaussian
        popt, pcov = curve_fit(_2gaussian, centers, ratioHist,   \
                                p0 =[a1, x01, sigma1,  a2, x02, sigma2], maxfev=10000)
        
        gauss_2 = _2gaussian(xSmooth, *popt)
        # Separate into 2 gaussians
        pars_g = popt[0:3]
        pars_n = popt[3:6]
        
        # Add parameters to output
        parsmat[i,0,:] = pars_g
        parsmat[i,1,:] = pars_g
        
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
            plt.xlim(x01-0.07,x01+0.07)
            plt.savefig('PSD_'+str(lowerErg)+'_'+str(upperErg), dpi=500)

    return ergvec, fomvec, parsmat
        
def _gaussian(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))
    
def _2gaussian(x, a1, x01, sigma1,  a2, x02, sigma2):
    return a1*np.exp(-(x-x01)**2/(2*sigma1**2)) +\
           a2*np.exp(-(x-x02)**2/(2*sigma2**2)) 
  
def DiscLine(x, a, b,c,):
    return a*np.exp( -b*x) + c
'''
def DiscLine(x, a, b, c, d):
    return a*x**3 + b*x**2 + c*x + d
 '''
    
    
    
    
    