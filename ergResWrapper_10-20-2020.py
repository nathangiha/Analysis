# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 13:23:23 2020

@author: giha

Wrapper for calculating energy resolution of 6-pillar glass H2DPI system
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from loader import DataLoad
from bar import Bar, Bars
from edgeCal import EdgeCal, listEdgeCal, peakCal, fitBackscatterPeak
from psd import CalNClip, MovingAvg, PSD_hist, FOM_plot, PSD_ergslice, \
DiscLineN, DiscLineG, barsPSD
from timeHist import doTimeHist
from coincFinder import CoincFind
from calibrateGainShift import linearCalibrate, listBarLinCal
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

def _gaussian(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))
    
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


plt.close('all')


results_dir = 'C:\\Users\\giha\\Documents\\Repos\\Analysis\\figures\\'
barNum = 6
coincWindow = 70
plotLinewidth = 1
figureSize = (7, 4)
binEdgesE = np.arange(0, 900, 10)
centersE = (binEdgesE[:-1] + binEdgesE[1:]) / 2

readIn = False
if readIn:
    lower = 87
    upper = 100
    
    data = []
      
    for i in range(lower, upper + 1):
        cache = DataLoad('D:\X' + str(i) + 'data.p')
        data.append(cache)

transform = False
if transform:
    datNaiCals = []
    datBefCals = []
    datAftCals = []
    datErgRes = []
    
    for i in range(4):
        datNaiCals.append(data[i * 4][12])
    
    for i in range(3):
        datBefCals.append(Bars(data[i * 4 + 1]))
        datAftCals.append(Bars(data[i * 4 + 3]))
        datErgRes.append(data[i * 4 + 2])


plotNaI = True
if plotNaI:
    fig, ax = plt.subplots()
    piThreshold = 5

    for i in range(len(datNaiCals)):
        dataToHist = datNaiCals[i][:, 0]
        dataToHist = dataToHist[dataToHist > piThreshold]
        
        hist = plt.hist(dataToHist,
                        bins = 500,
                        label = 'Cal. ' + str(i),
                        density = True,
                        histtype = 'step')
        
    
    plt.legend()
    ax.set_xlabel('Pulse Integral')
    ax.set_ylabel('Counts (norm.)')
    plt.tight_layout()
    plt.show()

'''
cBefWorking = datBefCals[1]
ergWorking = datErgRes[1]
cAftWorking = datAftCals[1]

befCals = listEdgeCal(cBefWorking)
aftCals = listEdgeCal(cAftWorking)

ergWorkingBar = Bars(ergWorking)
calInts = listBarLinCal(ergWorkingBar, befCals, aftCals)
'''

processBars = True

if processBars:
    barBackscatterEvents = []
        
    for meas in range(1, 3):
        cBefWorking = datBefCals[meas]
        ergWorking = datErgRes[meas]
        cAftWorking = datAftCals[meas]
    
        befCals = listEdgeCal(cBefWorking)
        aftCals = listEdgeCal(cAftWorking)
        
        ergWorkingBar = Bars(ergWorking)
        calInts = listBarLinCal(ergWorkingBar, befCals, aftCals)
    
        
        
        for bar in range(barNum):    
        
            # Setup for first bar
            b0Bef = befCals[bar]
            b0Aft = aftCals[bar]
            b0Erg = ergWorkingBar[bar]
            nai = ergWorking[12]
            
            # Time-dependent calibration
            b0Erg[0, :] = calInts[bar]
            
            # Get coincidences
            barTime = b0Erg[4, :] + b0Erg[5, :]
            naiTime = nai[:, 4] + nai[:, 5]
            coincInds = CoincFind(barTime, naiTime, coincWindow)
            
            barInd = coincInds[:,0]
            naiInd = coincInds[:,1]
                    
            barCoinc = b0Erg[:, barInd]
            naiCoinc = nai[naiInd, :]
            
            
            cfd = naiCoinc[:, 4] - barCoinc[4, :]
            ttt = naiCoinc[:, 5] - barCoinc[5, :]
            
            dt = cfd + ttt
            
            plt.figure(figsize = figureSize)
            binEdges = np.linspace(-coincWindow, coincWindow, 150)
            plt.hist(dt, bins = binEdges, histtype = 'step')
            
            # Cut waveforms outside of coincidence peak
            lowerTime = 40 # ns
            upperTime = 55 # ns
            windowIndices = np.logical_and(dt > lowerTime, dt < upperTime)
            
            barWindow = barCoinc[:, windowIndices]
            naiWindow = naiCoinc[windowIndices, :]
        
            cfdW = naiWindow[:, 4] - barWindow[4, :]
            tttW = naiWindow[:, 5] - barWindow[5, :]
            dtW = cfdW + tttW
            plt.hist(dtW, bins = binEdges)
            
            
            chanceLargerFactor = 4
            windowSize = upperTime - lowerTime
            
            
            upperTimeChance = lowerTime - windowSize
            lowerTimeChance = upperTimeChance - windowSize * chanceLargerFactor
            
            chanceIndices = np.logical_and(dt > lowerTimeChance,
                                           dt < upperTimeChance)
            
            
            
            barChance = barCoinc[:, chanceIndices]
            naiChance = naiCoinc[chanceIndices, :]
            
            cfdC = naiChance[:, 4] - barChance[4, :]
            tttC = naiChance[:, 5] - barChance[5, :]
            dtC = cfdC + tttC
            plt.hist(dtC, bins = binEdges)
            
            
            
            plt.xlim(left = -coincWindow, right = coincWindow)
            plt.xlabel(r'$\Delta$t')
            plt.ylabel(r'Counts/$\Delta$t')
            plt.tight_layout()
            plt.savefig(results_dir + 'meas' + str(meas) + 'bar' + str(bar) + 'backscatterTimeDif.png')
            
            
            naiCal = peakCal(nai[:,0], plot = False)
            naiEnergy = naiWindow[:,0] * naiCal[0]
            
            
            barEnergy = barWindow[0, :]
            barEnergyHist = np.histogram(barEnergy, bins = binEdgesE)[0]
    
            barEnergyChance = barChance[0, :]
            barEnergyChanceHist = np.histogram(barEnergyChance, bins = binEdgesE)[0]
                        
            backscatterPeak = 661.7 / (1 + 661.7 / 511 * 2 )
            comptonEdge = 661.7 - backscatterPeak
            
            resolution = 0.15
            backscatterUpper = backscatterPeak * (1 + resolution)
            backscatterLower = backscatterPeak * (1 - resolution)
            
            naiBackscatterInds = np.logical_and(naiEnergy > backscatterLower,
                                                naiEnergy < backscatterUpper)
            
            naiBackscatter = naiEnergy[naiBackscatterInds]
            barBackscatter = barEnergy[naiBackscatterInds]
            
            barBackscatterHist = np.histogram(barBackscatter,
                                              bins = binEdgesE)[0]

            # barBackscatterHistChanceSubt = (barBackscatterHist / np.sum(barBackscatterHist)  -
            #                 barEnergyChanceHist / np.sum(barEnergyChanceHist))
        
            barBackscatterHistChanceSubt = barBackscatterHist - barEnergyChanceHist
        
            
            
            # Plot NaI spectra
            fig, ax = plt.subplots(figsize = figureSize)
            plt.hist(naiEnergy,
                     bins = binEdgesE,
                     histtype = 'step')
            
            plt.hist(naiBackscatter,
                     bins = binEdgesE)
            ax.set_xlabel('Energy (keV)')
            ax.set_ylabel('Counts / keV')
            plt.tight_layout()
            plt.savefig(results_dir + 'meas' + str(meas) + 'bar' + str(bar) + 'backscatterNaI.png')
            
            
            # Plot bar spectra
            fig, ax = plt.subplots(figsize = figureSize)
            plt.hist(barEnergyChance,
                     bins = binEdgesE,
                     histtype = 'step',
                     label = 'Chance coinc.')
            
            plt.hist(barEnergy,
                     bins = binEdgesE,
                     histtype = 'step',
                     label = 'Coinc.')
            
            plt.hist(barBackscatter,
                     bins = binEdgesE,
                     histtype = 'step',
                     label = 'Gated')
            
            
            ax.set_xlabel('Energy (keV)')
            ax.set_ylabel('Counts / keV')
            ax.legend()
            
            fig, ax = plt.subplots(figsize = figureSize)
            plt.step(centersE,
                     barBackscatterHistChanceSubt,
                     linewidth = plotLinewidth,
                     c = 'k',
                     where = 'mid',
                     label = 'Chance subt.')

            
            peakCal(barBackscatter,
                    userBinEdges = binEdgesE,
                    userPeakErg = comptonEdge)
            
            plt.tight_layout()
            plt.savefig(results_dir + 'meas' + str(meas) + 'bar' + str(bar) + 'backscatterBar.png')
            
            barBackscatterEvents.append(barBackscatter)
            
            test = fitBackscatterPeak(barBackscatterHist,
                                      barEnergyChanceHist,
                                      binEdges = binEdgesE,
                                      unevenFactor = chanceLargerFactor)
            
            '''
            totEnergy = naiEnergy[naiBackscatterInds] + barEnergy[naiBackscatterInds]
            
            fig, ax = plt.subplots(figsize = figureSize)
            plt.hist(totEnergy,
                     bins = binEdgesE,
                     density = True)
            
            peakCal(totEnergy,
                    percentMaxThreshold = 0.5,
                    userBinEdges = binEdgesE)
            
            ax.set_xlabel('Energy (keV)')
            ax.set_ylabel('Counts / keV')
            
            plt.tight_layout()
            plt.savefig(results_dir + 'meas' + str(meas) + 'backscatterSum.png')
            '''

smash = True
if smash:
    barBackscatterSmash = []
    
    for i in range(barNum):
        thisBar = np.hstack((barBackscatterEvents[i],
                             barBackscatterEvents[i + 6]))
        barBackscatterSmash.append(thisBar)


plot = True
if plot:
    ergRes = np.zeros(barNum)
    for i in range(barNum):
        thisCal = peakCal(barBackscatterSmash[i],
                          userBinEdges = binEdgesE,
                          userPeakErg = comptonEdge)
        
        thisErgRes = thisCal[1]
        ergRes[i] = thisErgRes
        
        plt.savefig('C:\\Users\\giha\\Documents\\Papers\\glassBars\\' + 'bar' + str(i) + 'EnergyRes.png',
                    dpi = 300)

np.savetxt('C:\\Users\\giha\\Documents\\Repos\\Analysis\\results\\barsErgRes.txt',
           ergRes)

def thrownDataset():
    ergThrow = datErgRes[0]
    nai = ergThrow[12][:,0]
    numEvents = np.size(nai)
    
    eventBins = np.arange(0, numEvents, 10000)
    
    plt.figure()
    binEdges = np.arange(5, 600, 3)
    
    peakVec = np.zeros(len(eventBins)-1)
    
    
    for i in range(len(eventBins) - 1):
        dataToPlot = nai[eventBins[i]:eventBins[i + 1]]
        '''
        hist = plt.hist(dataToPlot,
                 bins = binEdges,
                 label = 'segment ' + str(i),
                 histtype = 'step')[0]
        '''
        if i % 1 == 0:
            calibFactor, fwhm, peakLoc = peakCal(dataToPlot)
            peakVec[i] = peakLoc
        
        
        
    plt.legend()
    
    plt.figure()
    plt.plot(peakVec[np.nonzero(peakVec)])
    plt.xlabel('Time (a.u.)')
    plt.ylabel(r'Peak location $(V\cdot ns)$')
    
    return peakVec


























