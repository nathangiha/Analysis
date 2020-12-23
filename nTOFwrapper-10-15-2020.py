# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 13:30:47 2020

@author: giha

Wrapper script for calling processing functions to process nTOF data

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import copy

from loader import DataLoad
from bar import Bars
from edgeCal import EdgeCal, listEdgeCal
from psd import CalNClip, MovingAvg, PSD_hist, \
DiscLineN, DiscLineG, barsPSD
from timeHist import doTimeHist
from coincFinder import CoincFind
from calibrateGainShift import linearCalibrate
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

def TOFtoEnergy(tof, pathlength): #tof in ns, pathlength in cm
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
    
def lengthToCTime(length): # Input length in cm
    return length / 2.99792e8 * 1e7 # Output time in ns

def identity(matrix):
    return matrix
##############################################################################
##############################################################################
##############################################################################
##############################################################################


# Declare gamma-misclass rate, effectively
psdSigma = 4

results_dir = 'C:\\Users\\giha\\Documents\\Repos\\Analysis\\figures\\'

# Declare flight paths
flightLengthFront = 48.65 # cm
flightLengthBack = 49.88 # cm

barDepth = 0.6 # cm
barUncert = barDepth / np.sqrt(12) # cm

flightLengthTag = 2.54 # cm

tagDepth = 2.54 # cm
tagUncert = tagDepth / np.sqrt(12) # cm

timeToFront = lengthToCTime(flightLengthFront)
timeToBack = lengthToCTime(flightLengthBack)
timeToTag = lengthToCTime(flightLengthTag)





barNum = 6
plotSlices = False

##############################################################################
##############################################################################
##############################################################################
##############################################################################


plt.close('all')

readIn = True
if readIn:
    X71 = DataLoad('D:\X71data.p')
    
    lower = 74
    upper = 85
    
    dataTof = []
      
    for i in range(lower, upper + 1):
        cache = DataLoad('D:\X' + str(i) + 'data.p')
        dataTof.append(cache)

    print('Done reading')
    '''
    X74 = DataLoad('D:\X74data.p')
    X75 = DataLoad('D:\X75data.p')
    X76 = DataLoad('D:\X76data.p')
    
    X77 = DataLoad('D:\X77data.p')
    X78 = DataLoad('D:\X78data.p')
    X79 = DataLoad('D:\X79data.p')
    
    X80 = DataLoad('D:\X80data.p')
    X81 = DataLoad('D:\X81data.p')
    X82 = DataLoad('D:\X82data.p')
    
    X83 = DataLoad('D:\X83data.p')
    X84 = DataLoad('D:\X84data.p')
    X85 = DataLoad('D:\X85data.p')
    '''
    tagCalData = X71[12]
    # tof1Tag = X75[12]

transform = True
if transform:
    calBefTofRaw = []
    calAftTofRaw = []
    nTofDataRaw = []
    tofTagRaw = []
    
    for i in range(4):
        calBefTofRaw.append(Bars(dataTof[i * 3]))
        calAftTofRaw.append(Bars(dataTof[i * 3 + 2]))
        nTofDataRaw.append(Bars(dataTof[i * 3 + 1]))
        tofTagRaw.append(dataTof[i * 3 + 1][12])
    

    # Data protection
    calBefTof = copy.deepcopy(calBefTofRaw)
    calAftTof = copy.deepcopy(calAftTofRaw)
    nTofData = copy.deepcopy(nTofDataRaw)
    tofTag = copy.deepcopy(tofTagRaw)
        
# currentFile = X75
fname = 'X75'

### Calibration and PSD for bars


calibrate = True

if calibrate:

    tagCal = EdgeCal(tagCalData[:, 0], histLabel = 'Tag',
                     xCal = 0, integral = True)[0]

    # Calibrate bars and tag in tof measurements
    for i in range(4):
        thisCalBefTof = listEdgeCal(calBefTof[i])
        thisCalAftTof = listEdgeCal(calAftTof[i])
        
        for j in range(barNum):
            pi = nTofData[i][j][0, :]
            nTofData[i][j][0, :] = linearCalibrate(pi,
                    thisCalBefTof[j], thisCalAftTof[j])
        
    
        tofTag[i][:, 0] *= tagCal
     
# Smash measurements together
smash = True
if smash:
    print('Smashing...')
    nTofDataCopy = copy.deepcopy(nTofData)
    tofTagCopy = copy.deepcopy(tofTag)
    
    timeTagSeparator = np.float64(5e15)
    plt.figure()
    for i in range(4):
        for j in range(barNum):
            nTofDataCopy[i][j][5, :] += timeTagSeparator * i
            
        tofTagCopy[i][:, 5] += timeTagSeparator * i
        plt.hist(nTofDataCopy[i][0][5, :], bins = 500)
    
    nTofDataSmashed = []
    for i in range(barNum):
        barStack = np.hstack((nTofDataCopy[0][i], nTofDataCopy[1][i], nTofDataCopy[2][i], nTofDataCopy[3][i]))
        nTofDataSmashed.append(barStack)
      
    plt.figure()
    plt.hist(nTofDataSmashed[0][5, :], bins = 500)
    
    tofTagSmashed = np.vstack((tofTagCopy[0], tofTagCopy[1], tofTagCopy[2], tofTagCopy[3]))
    plt.hist(tofTagSmashed[:, 5], bins = 500)
    
    print('Smashed.')
        
    
    '''
    tof1Bars = Bars(X75)   
    tagData = identity(tof1Tag)

    nTofData = tof1Bars


    # Calibration for bars and tag detector
    calBefTof1 = listEdgeCal(csTof1Bef)
    calAftTof1 = listEdgeCal(csTof1Aft)

    tagCalData = X71[12]
    plt.figure()
    tagCal = EdgeCal(tagCalData[:, 0], histLabel = 'Tag', xCal=0, integral=True)[0]
    
    tagData[:, 0] = tagData[:, 0] * tagCal
    
    # Calibrate bar light output
    for i in range(len(nTofData)):
        pi = nTofData[i][0, :]
        nTofData[i][0, :] = linearCalibrate(pi, calBefTof1[i], calAftTof1[i])
    '''        
        
        
    '''
    
    psdDataLO = []
    fitLOParams = []
    
    for i in range(len(nTOFData)):
        psdDataLO.append(CalNClip(nTOFData[i][0,:], nTOFData[i][3,:],
                                   csCalsTof[i]))
        fitLOParams.append(PSD_hist( psdDataLO[i][0], psdDataLO[i][1] , binnum = 300,\
                                   ergLimit=1000, discLine=psdSigma  )    \
                        )
        # plt.savefig(results_dir + fname+ 'bar' + str(i)+'PSD.png')
    '''
    

psd = True 
if psd:
    print('Applying PSD...')
    # PSD for bars and tag detector        
    nTofPsdData, nTofPsdFits = barsPSD(nTofDataSmashed, _discLine = psdSigma)
    del(nTofPsdData)
    
    tagPI = tofTagSmashed[:, 0]
    tagTails = tofTagSmashed[:, 2]
    tagTotals = tofTagSmashed[:, 3]
    tagRatio = tagTails / tagTotals
    tagPSD = CalNClip(tagPI, tagRatio)
    tagParams = PSD_hist(tagPSD[0], tagPSD[1], parType = 'g', ergLimit = 1000, discLine = psdSigma)
    '''
    tagEnergies = tagPI
    # Cut neutron events out of tag
    tagGammaInds = np.argwhere(tagRatio < DiscLineN(tagEnergies, *tagParams))[:, 0]
    
    tagGamma = tofTagSmashed[tagGammaInds, :]
    PSD_hist(tagGamma[:, 0], (tagGamma[:, 2] / tagGamma[:, 3]))
    '''    
    print('PSD applied.')
    
doLO = False
if doLO:
    lightOutputFit = []
    
    for i in range(barNum):
        barEvents = nTofDataSmashed[i]
        tagEvents = tofTagSmashed
        
        # Now, get rid of gamma events in bars, neutron events in tag
        barEnergy   = barEvents[0, :]
        barPSDRatio = barEvents[3, :]
        
        tagEnergy = tagEvents[:, 0]
        tagPSDRatio = tagEvents[:, 2] / tagEvents[:, 3]
            
        barNeutronInds = np.argwhere(barPSDRatio > DiscLineN(barEnergy, *nTofPsdFits[i]))[:, 0]
        tagGammaInds = np.argwhere(tagPSDRatio < DiscLineN(tagEnergy, *tagParams))[:, 0]
        
        tagGamma = tagEvents[tagGammaInds, :]
     
        
        
        tagTime = tagGamma[:, 4] + tagGamma[:, 5]
            
        # Use only neutron data in bars and only gamma data in tag
        barNEvents = barEvents[:, barNeutronInds]
    
        PSD_hist(barNEvents[0, :], barNEvents[3, :])
        plt.savefig(results_dir + fname + 'bar' + str(i) + 'PSDcut.png')
    
        barTime = barNEvents[4, :] + barNEvents[5, :]    
        
        # Set coincidence window between bar and tag detector
        barTagWindow = 80 #ns
    
        tagC   = CoincFind(barTime, tagTime, barTagWindow)
        barInd = tagC[:, 0]
        tagInd = tagC[: ,1]
        
        barCoinc = barNEvents[: ,barInd]
        tagCoinc = tagGamma[tagInd, :]
    
    
        cfd = barCoinc[4,:] - tagCoinc[:,4]
        ttt = barCoinc[5,:] - tagCoinc[:,5]
        
        dt = ttt + cfd
        meanDt = np.mean(dt)
        
        
          
        # Plot tdif histogram
        timeBinEdges = np.arange(-10, barTagWindow, fwhm[i])
        timeBinCenters = (timeBinEdges[:-1] + timeBinEdges[1:]) / 2
        
        # Correct dt for time offsets
        if i < 4:
            flightPath = flightLengthFront
        else:
            flightPath = flightLengthBack
        
        neutronFlightTime = dt - offset[i] + lengthToCTime(flightPath)
        
        fig, ax = plt.subplots(figsize = (7, 4))
        dtHist = plt.hist(neutronFlightTime, bins = timeBinEdges)[0]
        #plt.hist(dt , bins=500)
        #plt.xlim(-30, 100)
        maxTime = timeBinCenters[np.argmax(dtHist)]
        maxErg = TOFtoEnergy(maxTime, flightPath)
        plt.title('Bar '+ str(i) + r', maxErg = ' + str(maxErg))
        plt.xlabel(r'$\Delta$t (ns)')
        plt.ylabel('Counts')
        
        plt.tight_layout()
        plt.savefig(results_dir + fname + 'bar' + str(i) + 'TOF.png')
    
            
        # Determine minimum time bin we're concerned with (corresponds to max. energy)
        # 10 MeV neutron has flight time of ~11 ns, let's set the min. to 10
        
        binEdgesRevised = np.arange(10, barTagWindow, fwhm[i])
        timeBinCenters =  (binEdgesRevised[:-1] + binEdgesRevised[1:]) / 2 
                
        energies = TOFtoEnergy(timeBinCenters, flightPath)
        lightOutputVec = np.zeros( len(energies) )
        
        lowerEnergy = 1 # MeV
        upperEnergy = 5 # MeV
        
        upperTimeIndex = np.argmax(energies < lowerEnergy)
        lowerTimeIndex = len(energies) - np.argmax(energies[::-1] > upperEnergy)
        
        
        
        
        # Find LOs of interactions in each time bin 
        for j in range( lowerTimeIndex, upperTimeIndex + 1 ): #len( binEdgesRevised ) - 1 ):
            
            # Find upper and lower bound of bins
            lowerBound = binEdgesRevised[j]
            upperBound = binEdgesRevised[j + 1]
            lowerEnergy = TOFtoEnergy(upperBound, flightPath)
            upperEnergy = TOFtoEnergy(lowerBound, flightPath)
            
            # Find data from bin
            indicesInBin = np.argwhere(np.logical_and(neutronFlightTime > lowerBound, neutronFlightTime < upperBound))
            lightOutputs = np.sort(barCoinc[0, indicesInBin])
            
    
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
            
            if not plotSlices:
                plt.close()
            
            lightOutputVec[j] = minDerivLightOutput# lightOutput
            
            
            
            
        plt.figure()
        
        
        lightOutputMeV =  lightOutputVec[lightOutputVec > 0] * 1e-3
        
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
        aFit, bFit = FitLightOutput(energiesToPlot, lightOutputsToPlot, 1, 1)
        
        lightOutputFit.append( [aFit, bFit] )
        
        ergSmooth = np.linspace(0,max(energiesToPlot) * 1.2, 1e4)
        
        plt.plot(ergSmooth, aFit * np.power(ergSmooth, bFit), c = 'k', \
                 label = 'Fit' )
        plt.savefig(results_dir + fname+ 'bar' + str(i)+'_' + 'nLO' + '.png')
            
            
        
            
    np.savetxt(results_dir + 'lightOutputCurves.txt', lightOutputFit)


