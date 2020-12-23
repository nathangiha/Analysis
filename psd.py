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
from scipy.signal import find_peaks
from scipy.stats import zscore

#data = X14

#ph = data[0][:,0]
#conv_factor = 0.18835

#data = B1cf

def barsPSD(barList, _binnum=300, _ergLimit=1000, _discLine=0, _binsize=50):
    numFits = len(barList)
    print(numFits)
    
    psdData = []
    fitParams = []
    

    for i in range(numFits):
        barPSD = CalNClip(barList[i][0,:], barList[i][3,:])
        psdData.append( barPSD )
        fitParams.append(PSD_hist(barPSD[0], barPSD[1], binnum = _binnum, ergLimit = _ergLimit, discLine = _discLine))
        
    return psdData, fitParams

# Calibrate pulse height/integral to light output. clip at specified upper limit
def CalNClip(intg, ratio, conv_factor = 1, ergLimit = 1500):
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

# Generates 2D PSD histogram
def PSD_hist(pi, ratio, parType = 'n', binnum=300, ergLimit=1500, discLine=0, binsize=50):
    plt.figure()
    
    
    
    
    intp = pi[(ratio >=0) & (ratio <=1) ]
    ratiop = ratio[(ratio >=0) & (ratio <=1) ]
    
#    psd = plt.hist2d(ph,ratio, bins=(binnum,50*binnum ), cmap=plt.cm.jet,
 #                    norm=LogNorm())
    psd = plt.hist2d(intp, ratiop, bins = (binnum, binnum * 10),\
        range = [[0, ergLimit], [0, 1]], cmap = plt.cm.jet, norm = LogNorm())
    
    maxInd = np.where(psd[0] == np.amax(psd[0]))
    distLoc = psd[2][maxInd[1][0]]
    plt.clf()
    plt.xlim(0, ergLimit)
    
    # Set ratio bounds
    lowerRatio = distLoc - 0.2
    upperRatio = distLoc + 0.2
    plt.ylim(lowerRatio, upperRatio)
    plt.close()


    plt.figure(figsize = (7, 4))
    
    # Inltialize fit params
    fitParams = []
    # If we are plotting a discrimination line (where the value of discLine is
    # the number of standard deviations away from the gamma Gaussian)
    if discLine > 0:
        slices = PSD_ergslice(pi, ratio, ergbin = binsize, maxerg = ergLimit, plot = 0)
        erg = slices[2]
        parsMat = slices[4]
        gPars = parsMat[:,0,:]
        nPars = parsMat[:,1,:]
        
        plotSigma = 6
        lowerRatio = gPars[2,1] - plotSigma * gPars[2,2]
        upperRatio = nPars[2,1] + plotSigma * nPars[2,2]
        
        if parType == 'g':
            discPoints = nPars[:,1] - discLine * nPars[:,2] # Mean minus some stdevs
        else:            
            discPoints = gPars[:,1] + discLine * gPars[:,2] # Mean plus some stdevs
        #plt.autoscale(False)
        #plt.scatter(erg, discPoints)
        
        # Remove outliers prior to fitting
        zScoreMin = -5
        zScoreMax = 5
        zScores = np.array(zscore( discPoints ).flatten() )
        
        # print(zScores)
        
        upperCut = np.nonzero(zScores < zScoreMax)[0]
        lowerCut = np.nonzero(zScores > zScoreMin)[0]
        usableIndices = np.intersect1d(upperCut, lowerCut)
        
        ergUsing = np.array(erg)[usableIndices]
        discPointsUsing = np.array(discPoints)[usableIndices]
        
        deriv = np.diff(discPointsUsing)
        '''
        print(deriv)
        print(max(deriv))
        '''
        '''
        derivProhibitedIndices = np.argmax(np.abs(deriv) > 0.01)
        
        ergUsing = np.delete(ergUsing,derivProhibitedIndices)
        discPointsUsing = np.delete(discPointsUsing, derivProhibitedIndices)
        '''
        
        if len(discPointsUsing) < len(discPoints):
            print('Warning: Threw out ' + str( len(discPoints) - \
                                              len(discPointsUsing) ) + ' slice(s)')
        
        plt.scatter(ergUsing, discPointsUsing, zorder = 1, s = 15, c='k')
        plt.errorbar(ergUsing, discPointsUsing, xerr = binsize/2, \
                     capsize = 5,fmt = 'none', zorder = 1, c='k')
        
        numPointsToPlot = len(discPointsUsing)       
        fitWeights = np.ones(numPointsToPlot)
        
        
        for i in range(numPointsToPlot):
            fitWeights[i] = 0.5 + np.exp(- i)
            
        fitWeights[0] = np.max(fitWeights) * 100
        # print(fitWeights)
        # print(len(fitWeights))
        
        
        # Perform discrimination line fit
        if parType == 'g':
            guessIntercept = nPars[1,1] - discLine * nPars[1,2]
            # print('Intercept guess: ' + str(guessIntercept))

            popt, pcov = curve_fit(DiscLineN, ergUsing, discPointsUsing,   \
                                    p0 = [-0.15, 0, -0.15, 0.016, guessIntercept], \
                                    sigma = fitWeights, maxfev=100000)
            # print(popt)
            eSmooth = np.arange(0,np.max(ergUsing),1)
            discFit = DiscLineN(eSmooth, *popt)
            fitParams = popt
            
        else:
            guessIntercept = gPars[1,1] + discLine * gPars[1,2]

            popt, pcov = curve_fit(DiscLineN, ergUsing, discPointsUsing,   \
                                    p0 = [0.1, 0.015, 0, 0, guessIntercept], \
                                    sigma = fitWeights, maxfev=100000)
            
            eSmooth = np.arange(0,np.max(ergUsing),1)
            discFit = DiscLineN(eSmooth, *popt)
            fitParams = popt
        
        plt.plot(eSmooth, discFit, 'k', label = str(discLine)+ r'$\sigma$')

    
    
    psd = plt.hist2d(intp,ratiop, bins=(binnum, binnum),\
        range = [[0, ergLimit],[lowerRatio, upperRatio]], cmap=plt.cm.jet,\
        zorder = 0, norm=LogNorm())
    '''
    
    psd = plt.hist2d(intp,ratiop, bins=(binnum,binnum ),\
        range = [[0,ergLimit],[0, 1]], cmap=plt.cm.jet,\
        zorder = 0, norm=LogNorm())
    '''
    plt.colorbar()
    plt.xlabel(r'Light output (keVee)')
    plt.ylabel(r'$r_{PSD}')
    plt.legend()
    plt.tight_layout()
    return fitParams
 
# Plots PSD figure-of-merit as a function of light output
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

    
# Fits double gaussians PSD energy slices
def PSD_ergslice( pi, ratio, ergbin=50, minerg = 0, maxerg = 1000, cal=1, plot=True ):
    
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
        binNum = 500
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
        
        '''
        xg = np.argmax(smoothedratioHist)+3
        
        locmaxvec[i] = locmax[0][0]
        
        try:
            xn = locmax[0][np.argmin(  np.abs(locmax[0] - xg ) ) +1]+3
        except:
            xn = xg+13
        '''
        
        
        
        xg = locmax[0]+1
        try:
            xn = locmax[1]+1
        except:
            xn = xg + 10
        
        # While loop acts as high-pass filter in histogram, filtering out
        # low-amplitude noise before choosing maxima
        
        '''
        while (ratioHist[xg]< 15):
            xg = locmax[j]+1
            j = j+1
        try:
            xn = locmax[0][j+1]+1
        except:
            xn = xg+20
            
            
        # Band-aid solution to the argrelmax function sucking sometimes
        if xn > xg+50:
            xn = xg+20
    
            
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



        # Print energy and max. found location
        '''
        print(str(upperErg/100) + ' ' + str(x01) + ' ' +str(x02))    
        print(str(locmax[0]))
        '''
        
        
        
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
        
        # Go to next iteration if gamma peak is smaller than neutron peak
        if pars_g[0] < pars_n[0]:
            continue
    
        # Add parameters to output
        parsmat[i,0,:] = pars_g
        parsmat[i,1,:] = pars_n
        
        gauss_g = _gaussian(xSmooth, *pars_g)
        gauss_n = _gaussian(xSmooth, *pars_n)
        
        fom = (pars_n[1] - pars_g[1]) / (2 * np.sqrt(2 * np.log(2)) *(abs(pars_g[2])+abs(pars_n[2])))     
        fomvec[i] = fom
        ergvec[i] = (lowerErg + upperErg)/2
        
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
            plt.xlim(x01-0.1,x01+0.3)
            plt.savefig('PSD_'+str(lowerErg)+'_'+str(upperErg), dpi=500)

    return [centers, ratioHist, ergvec, fomvec, parsmat]
        
def _gaussian(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))
    
def _2gaussian(x, a1, x01, sigma1,  a2, x02, sigma2):
    return a1*np.exp(-(x-x01)**2/(2*sigma1**2)) +\
           a2*np.exp(-(x-x02)**2/(2*sigma2**2)) 

def DiscLineN(x, a, b, a2, b2, c):
    return a*np.exp( -b*x) + a2*np.exp(-b2*x) + c
'''
def DiscLineG(x, a, b,c,):
    return -a*np.exp( -b*x) + c
'''
def DiscLineG(x, a, b, c, d):
    return a*x**3 + b*x**2 + c*x + d
 
    
    
    
    
    