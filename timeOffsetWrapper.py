# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 15:53:37 2020

@author: giha
Wrapper script for finding timing offsets of bars.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from loader import DataLoad
from bar import Bar, Bars
from edgeCal import EdgeCal, listEdgeCal
from psd import CalNClip, MovingAvg, PSD_hist, FOM_plot, PSD_ergslice, \
DiscLineN, DiscLineG, barsPSD
from timeHist import TimeHist
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
    
def lengthToCTime(length): # Input length in cm
    return length / 2.99792e8 * 1e7 # Output time in ns


##############################################################################
##############################################################################
##############################################################################
##############################################################################


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
fname = 'X73'

##############################################################################
##############################################################################
##############################################################################
##############################################################################


plt.close('all')

readIn = True
if readIn:
    X71 = DataLoad('D:\X71data.p')
    X72 = DataLoad('D:\X72data.p')
    X73 = DataLoad('D:\X73data.p')
    X74 = DataLoad('D:\X74data.p')

    tagCalData = X71[12]
    csOffsetBef = Bars(X72)
    naTimeRes = X73
    naBars = Bars(X73)
    naTag = X73[12]
    csOffsetAft = Bars(X74)

    print('Done reading')
    
### Calibration and PSD for bars

calibrate = True

if calibrate:
    # Calibration for bars and tag detector
    calOffsetBef = listEdgeCal(csOffsetBef)
    calOffsetAft = listEdgeCal(csOffsetAft)
    
    plt.figure()
    tagCal = EdgeCal( tagCalData[:,0], histLabel = 'Tag', xCal=0, integral=True)

    
    # Calibrate bar light output, assuming linear gain shift
    for i in range(len(naBars)):
        pi = naBars[i][0,:]
        naBars[i][0,:] = linearCalibrate(pi, calOffsetBef[i], calOffsetAft[i] )
        
        
        
        
initVars = True
if initVars:
    offset = []
    fwhm = []


for i in range(barNum):
    
    
    # Histogram intra-bar coincidences
    barBot = naTimeRes[i * 2]
    barTop = naTimeRes[i * 2 + 1]
    tag    = naTimeRes[12]


    barBotTime = barBot[:, 4] + barBot[:, 5]
    barTopTime = barTop[:, 4] + barTop[:, 5]


    intraBarCoinc = CoincFind(barBotTime, barTopTime, 5)
    
    botInd = intraBarCoinc[:, 0]
    topInd = intraBarCoinc[:, 1]

    bot = barBot[botInd]
    top = barTop[topInd]

    dtBar = (bot[:, 4] - top[:, 4]) + (bot[:, 5] - top[:, 5])
    
    mean = np.mean(dtBar)

    plt.figure()
    plt.hist(dtBar[np.abs(dtBar) < 5], bins =  100, label= r'$\mu=$ ' + str(mean))
    plt.xlabel(r'$\Delta$t (ns)')
    plt.title('Bar ' + str(i))
    plt.legend()
    plt.tight_layout()
    plt.savefig(results_dir + fname + 'bar' + str(i) + 'intra.png')


    # Get bar data
    bar = naBars[i]
    
    # Set coincidence window between bar and tag detector
    barTagWindow = 60 #ns

    # Find coincidences between bars and tag det.
    tagTime = tag[:, 4] + tag[:, 5]
    barTime = bar[4, :] + bar[5, :]
    tagC = CoincFind(barTime, tagTime, barTagWindow)
    
    barInd = tagC[:,0]
    tagInd = tagC[:,1]
    
    barCoinc = bar[:, barInd]
    tagCoinc = tag[tagInd, :]
           
    
    cfd = barCoinc[4,:] - tagCoinc[:,4]
    ttt = barCoinc[5,:] - tagCoinc[:,5]
    
    dt = ttt + cfd
    meanDt = np.mean(dt)
    '''
    # Correct time difference for photon flight paths
    if i < 4:
        timeToBar = timeToFront
    else:
        timeToBar = timeToBack
    
    dtCorrected = dt - timeToBar + tagTime
    '''
    
    dtPos = dt[dt > 0]
    plt.figure()
    
    binEdges = np.arange(-10 + meanDt, 10 + meanDt, 0.1)
    centers = (binEdges[:-1] + binEdges[1:]) / 2
    hist = plt.hist(dt, bins = binEdges)[0]
    plt.close()
    plt.figure()
    
    # Gaussian guesses
    a = np.max(hist)
    x0 = meanDt
    sigma = 0.5
    
    # Fit gaussian
    popt, pcov = curve_fit(_gaussian, centers, hist,   \
                            p0 =[a, x0, sigma], maxfev=10000)

    xSmooth = np.arange(centers[0], centers[-1], 0.001)
    gauss_g = _gaussian(xSmooth, *popt)
    fitMean = popt[1]
    hist = plt.hist(dt, bins = binEdges, label = r'$\mu =$' + str(np.round(fitMean, 5)) + ' ns')[0]
    
    thisFWHM = FWHM(centers, hist)
    fitFWHM = FWHM(xSmooth, gauss_g)
    
    plt.plot(xSmooth, gauss_g, color = 'k', label = 'Gaussian fit')
    plt.xlim(-5 + fitMean, 5 + fitMean)
    plt.title('Bar '+ str(i) + '\nFWHM = ' + str(round(thisFWHM, 4)) + ' ns')
    plt.legend()
    plt.xlabel(r'$\Delta$t (ns)')
    plt.ylabel('Counts')
    plt.tight_layout()
    plt.savefig(results_dir + fname+ 'bar' + str(i)+'TOF.png')

    offset.append(fitMean)
    fwhm.append(thisFWHM)
    
offset = np.array(offset)
fwhm = np.array(fwhm)      
            
