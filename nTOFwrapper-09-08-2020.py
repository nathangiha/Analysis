# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 12:59:54 2020

@author: giha
Wrapper script for calling processing functions to process nTOF data
"""

import numpy as np
from loader import DataLoad
from bar import Bar, Bars
from edgeCal import EdgeCal
from psd import CalNClip, MovingAvg, PSD_hist, FOM_plot, PSD_ergslice, DiscLine
from timeHist import TimeHist
import matplotlib.pyplot as plt
from coincFinder import CoincFind
from scipy.optimize import curve_fit

##############################################################################
##############################################################################
##############################################################################
##############################################################################

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

def TOFtoEnergy(tof, pathlength = 36.2 ): #tof in ns, pathlength in cm
    massNeutron = 1.008665 # amu
    
    convFactor = 1.03642697 # MeV per amu * (cm/ns)^2
    
    energy = 0.5 * massNeutron * (pathlength / tof)**2   * convFactor
    
    return energy
    
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


def _powerLaw(E, a, b):
    return a* np.power(E, b)

    
def FitLightOutput(energy, lightOutput, aGuess, bGuess):
        popt, pcov = curve_fit(_powerLaw, energy, lightOutput,   \
                                p0 =[aGuess, bGuess], maxfev=10000)
        return popt
    
def _gaussian(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))
    
    


##############################################################################
##############################################################################
##############################################################################
##############################################################################


# Declare gamma-misclass rate, effectively
psdSigma = 4

results_dir = 'C:\\Users\\giha\\Documents\\Repos\\Analysis\\figures\\'


flightTime = 1.207 #1.356 #ns

barNum = 4
plotSlices = False

##############################################################################
##############################################################################
##############################################################################
##############################################################################


'''
X27 = DataLoad('D:\X27data.p')
X28 = DataLoad('D:\X28data.p')
X29 = DataLoad('D:\X29data.p')
X30 = DataLoad('D:\X30data.p')
X31 = DataLoad('D:\X31data.p')
X32 = DataLoad('D:\X32data.p')
'''

plt.close('all')

readIn = False
if readIn:
    X27 = DataLoad('D:\X27data.p')
    X30 = DataLoad('D:\X30data.p')
    X29 = DataLoad('D:\X29data.p')
    X28 = DataLoad('D:\X28data.p')
    csCalData = Bars(X27)
    nTOFData = Bars(X30)    
    tagCalData = X28


    print('Done reading')
currentFile = X30
fname = 'X30'

### Calibration and PSD for bars


calibrate = False

if calibrate:
    csCalsTof = []
    for i in range(len(csCalData)):
        csCalsTof.append( EdgeCal( csCalData[i][0,:], histLabel = 'Bar '+str(i), xCal=0, integral = True ) )
        plt.savefig(results_dir + fname + 'csCals.png')
    
    
    
    psdDataLO = []
    fitLOParams = []
    for i in range(len(nTOFData)):
        psdDataLO.append( CalNClip(nTOFData[i][0,:], nTOFData[i][3,:], csCalsTof[i][0] ) )
        fitLOParams.append( PSD_hist( psdDataLO[i][0], psdDataLO[i][1] , \
                                   ergLimit=1400, discLine=psdSigma  )    \
                        )
        plt.savefig(results_dir + fname+ 'bar' + str(i)+'PSD.png')
   
    
    ### Calibration and PSD for tag detector    
    tagCalData = X28[8]
    plt.figure()
    tagCal = EdgeCal( tagCalData[:,0], histLabel = 'Tag', xCal=0, integral=True)
    
    tagData = X30[8]
    '''
    tagPI = tagData[:,0]
    tagTails = tagData[:,2]
    tagTotals = tagData[:,3]
    tagRatio = tagTails / tagTotals
    tagPSD = CalNClip( tagPI, tagRatio, tagCal[0])
    PSD_hist( tagPSD[0], tagPSD[1], ergLimit=1000, discLine=4  )
    '''




if fname == 'X29':
    offset = []
    fwhm = []

if fname == 'X30':
    lightOutputFit = []

for i in range(barNum):
    
    
    # Histogram intra-bar coincidences
    barBot = currentFile[i*2]
    barTop = currentFile[i*2+1]
    tag    = currentFile[8]


    barBotTime = barBot[:,4] + barBot[:,5]
    barTopTime = barTop[:,4] + barTop[:,5]


    b1Coinc = CoincFind(barBotTime, barTopTime, 5)
    botInd = b1Coinc[:,0]
    topInd = b1Coinc[:,1]

    bot = barBot[botInd]
    top = barTop[topInd]

    dtBar = (bot[:,4] - top[:,4]) + (bot[:,5] - top[:,5])
    
    mean = np.mean(dtBar)

    plt.figure()
    plt.hist(dtBar[ np.abs(dtBar) < 5], bins=100, label= r'$\mu=$ ' + str(mean))
    plt.xlabel(r'$\Delta$t (ns)')
    plt.title('Bar ' + str(i))
    plt.legend()
    plt.tight_layout()
    plt.savefig(results_dir + fname+ 'bar' + str(i)+'intra.png')






    # Histogram bar-tag coincidences
    barAllEvents = Bar(bot,top)
    
    # Now, get rid of gamma events
    barLO       = barAllEvents[0,:]*csCalsTof[i][0]
    barPSDRatio = barAllEvents[3,:]
    
    # Fit parameters
    a = fitLOParams[i][0]
    b = fitLOParams[i][1]
    c = fitLOParams[i][2]
    
    neutronInds = np.argwhere(  barPSDRatio > DiscLine(barLO, a, b, c  ) )[:,0]
 
    tagTime = tag[:,4] + tag[:,5]
    
    
    
    
    # Set coincidence window between bar and tag detector
    barTagWindow = 60 #ns
    
    
    # Use all data if it's time res. measurement
    if fname == 'X29':
        barTime = barAllEvents[4,:] + barAllEvents[5,:]
        tagC = CoincFind(barTime, tagTime, barTagWindow)
        barInd = tagC[:,0]
        tagInd = tagC[:,1]
        
        barCoinc = barAllEvents[:,barInd]
        botCoinc = bot[barInd,:]
        topCoinc = top[barInd,:]
        tagCoinc = tag[tagInd, :]
        
        
        
    # Use only neutron data if it's TOF   
    if fname =='X30':
        barNEvents = barAllEvents [ :,  neutronInds     ]
    
        PSD_hist( barNEvents[0,:]*csCalsTof[i][0], barNEvents[3,:]   )
        plt.savefig(results_dir + fname+ 'bar' + str(i)+'PSDcut.png')

        barTime = barNEvents[4,:] + barNEvents[5,:]    
    
        tagC   = CoincFind(barTime, tagTime, barTagWindow)
        barInd = tagC[:,0]
        tagInd = tagC[:,1]
        
        barCoinc = barNEvents[:,barInd]
        botCoinc = (bot[neutronInds,:])[barInd,:]
        topCoinc = (top[neutronInds,:])[barInd,:]
        tagCoinc = tag[tagInd, :]

    
    
    
    
    ''' 
    ttt = ( botCoinc[:,4] + topCoinc[:,5] ) / 2.0 - tagCoinc[:,4]
    cfd = ( botCoinc[:,4] + topCoinc[:,5] ) / 2.0 - tagCoinc[:,5]
    '''   
    cfd = barCoinc[4,:] - tagCoinc[:,4]
    ttt = barCoinc[5,:] - tagCoinc[:,5]
    
    dt = ttt + cfd
    meanDt = np.mean(dt)
    
    
    
    dtp = dt[dt>0]
    plt.figure()
    
    if fname == 'X29':
        binedges = np.linspace(-5+meanDt,5+meanDt,301)
        centers = (binedges[:-1] + binedges[1:]) / 2
        hist = plt.hist(dt, bins=binedges)[0]
        plt.close()
        plt.figure()
        # Gaussian guesses
        a = np.max(hist)
        x0 = meanDt
        sigma = 0.5
        
        # Fit double gaussian
        popt, pcov = curve_fit(_gaussian, centers, hist,   \
                                p0 =[a, x0, sigma], maxfev=10000)
    
        gauss_g = _gaussian(centers, *popt)
        fitMean = popt[1]
        hist = plt.hist(dt, bins=binedges, label = r'$\mu =$' + str(np.round(fitMean,5)))[0]
        
        plt.plot(centers, gauss_g, color = 'k', label = 'Fit')
        plt.xlim(-5 + fitMean, 5 + fitMean)
        plt.title('Bar '+ str(i))
        plt.legend()
        plt.xlabel(r'$\Delta$t (ns)')
        plt.ylabel('Counts')
        plt.tight_layout()
        plt.savefig(results_dir + fname+ 'bar' + str(i)+'TOF.png')
    
        offset.append(fitMean)
        fwhm.append( FWHM(centers, hist )   )
      
    if fname == 'X30':
        # Plot tdif histogram
        timeBinEdges = np.arange( 0, barTagWindow, fwhm[i] )
        timeBinCenters = (timeBinEdges[:-1] + timeBinEdges[1:]) / 2
        
        dtCorrected = dt - offset[i] + flightTime
        
        dtHist = plt.hist( dtCorrected, bins=timeBinEdges )[0]
        #plt.hist(dt , bins=500)
        #plt.xlim(-30, 100)
        maxTime = timeBinCenters[ np.argmax( dtHist ) ]
        maxErg = TOFtoEnergy(maxTime)
        plt.title('Bar '+ str(i) + r', maxErg = ' + str(maxErg))
        plt.xlabel(r'$\Delta$t (ns)')
        plt.ylabel('Counts')
        
        plt.tight_layout()
        plt.savefig(results_dir + fname+ 'bar' + str(i)+'TOF.png')
        
        # Determine minimum time bin we're concerned with (corresponds to max. energy)
        # 10 MeV neutron has flight time of ~9 ns, let's set the min. to 10
        
        binEdgesRevised = np.arange(10, barTagWindow, fwhm[i])
        timeBinCenters =  (binEdgesRevised[:-1] + binEdgesRevised[1:]) / 2 
        
        energies = TOFtoEnergy( timeBinCenters )
        lightOutputVec = np.zeros( len(energies) )
        
        lowerEnergy = 1 # MeV
        upperEnergy = 5 # MeV
        
        upperTimeIndex = np.argmax(energies < lowerEnergy )
        lowerTimeIndex = len(energies) - np.argmax(energies[::-1] > upperEnergy)
        
        
        
        
        # Find LOs of interactions in each time bin 
        for j in range( lowerTimeIndex, upperTimeIndex + 1 ): #len( binEdgesRevised ) - 1 ):
            
            # Find upper and lower bound of bins
            lowerBound = binEdgesRevised[j]
            upperBound = binEdgesRevised[j+1]
            lowerEnergy = TOFtoEnergy(tof = upperBound)
            upperEnergy = TOFtoEnergy(tof = lowerBound)
            
            # Find data from bin
            indicesInBin = np.argwhere(np.logical_and( dtCorrected > lowerBound, dtCorrected < upperBound     ) )
            lightOutputs = np.sort( barCoinc[0,indicesInBin]*csCalsTof[i][0] )
            

            # Set LO bins
            # Choose max bin based on percentage of data to include
            sortedLO = np.sort(lightOutputs, axis=None)
            maxBin = sortedLO[ int(np.floor( len(sortedLO) * 0.99 )) ]
            lightOutputBinEdges = np.linspace(0, maxBin + 500, 41)
            lightOutputCenters = (lightOutputBinEdges[:-1] + lightOutputBinEdges[1:]) / 2
          
            # Plot LO hist for time bin
            plt.figure()  
            lightOutputHist = plt.hist(lightOutputs, bins=lightOutputBinEdges, histtype = 'step')[0]


            # Draw vertical line at some percent of total counts
            countSum = 0
            pickoffIndex = 0
            while countSum < 0.88 * np.sum(lightOutputHist ):
                countSum += lightOutputHist[pickoffIndex]
                pickoffIndex +=1
            
            lightOutputPercent = lightOutputBinEdges[pickoffIndex]
            plt.axvline(x = lightOutputPercent )


            # Smooth hist, find deriv, and pickoff edge            
            lightOutputSmoothed = MovingAvg(lightOutputHist)
            lightOutputDeriv = np.diff(lightOutputSmoothed)
            
            lightOutputDeriv = MovingAvg(lightOutputDeriv)

            
            
            minDerivLightOutput = lightOutputCenters[3:-2][ np.argmin(lightOutputDeriv)]
            
            while minDerivLightOutput < lightOutputPercent:
                lightOutputDeriv[ np.argmin(lightOutputDeriv) ] = 0
                minDerivLightOutput = lightOutputCenters[3:-2][ np.argmin(lightOutputDeriv)]
            
            # Weight actual min by 
           # minDerivBefore = lightOutputCenters[2:-1][ np.argmin(lightOutputDeriv)
            
            plt.axvline( x = minDerivLightOutput , c = 'k', label = 'min. deriv.')
            
            plt.plot(lightOutputCenters[1:-1], lightOutputSmoothed, label = 'Smoothed')
            plt.plot(lightOutputCenters[3:-2], lightOutputDeriv, label = 'Deriv.')
            
            
            plt.title( str( np.round(lowerBound, 3) ) + ' - ' + \
                      str( np.round(upperBound, 3) ) + ' ns' + '\n ' + \
                      str( np.round(lowerEnergy, 3) ) + ' - ' + \
                      str( np.round(upperEnergy, 3) ) + ' MeV'
                      )
            plt.xlabel('Light Output (keVee)')
            plt.legend()
            plt.savefig(results_dir + fname+ 'bar' + str(i)+'_' + \
                        str( round(lowerEnergy,3) ) +'-' \
                        + str(round(upperEnergy,3)) + 'mev.png')

            if ~plotSlices:
                plt.close()
            
            lightOutputVec[j] = minDerivLightOutput# lightOutput
            
            
            
            
        plt.figure()
        lightOutputMeV =  lightOutputVec [ lightOutputVec > 0] * 1e-3
        
        
        energiesToPlot = energies[lightOutputVec > 0][1:]
        lightOutputsToPlot = lightOutputMeV[1:]
        
        
        plt.scatter(energiesToPlot, lightOutputsToPlot, c = 'k', label = 'Data')
        # plt.plot([0,10], [0,10])
        plt.xlim( left = 0, right = max(energiesToPlot) * 1.2 )
        plt.ylim( bottom = 0 , top = max(lightOutputsToPlot) * 1.2)
        plt.xlabel(r'E$_n$ (MeV)')
        plt.ylabel('LO (MeVee)')
        plt.tight_layout()
        
        # Plot fit
        aFit, bFit = FitLightOutput( energiesToPlot, lightOutputsToPlot, 1, 1)
        
        lightOutputFit.append( [aFit, bFit] )
        
        ergSmooth = np.linspace(0,max(energiesToPlot) * 1.2, 1e4)
        
        plt.plot(ergSmooth, aFit * np.power(ergSmooth, bFit), c = 'k', \
                 label = 'Fit' )
        plt.savefig(results_dir + fname+ 'bar' + str(i)+'_' + 'nLO' + '.png')
        
        
        
        
if fname == 'X30':
    np.savetxt( results_dir + 'lightOutputCurves.txt', lightOutputFit)        
            
            
            
        
        








    
