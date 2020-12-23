# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 19:16:15 2020

@author: giha
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import curve_fit
import copy

from loader import DataLoad
from bar import Bar, Bars
from edgeCal import EdgeCal, listEdgeCal
from calibrateGainShift import linearCalibrate

plt.close('all')

results_dir = 'C:\\Users\\giha\\Documents\\Repos\\Analysis\\figures\\'


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

def _3rdOrder(x, a, b, c, d):
    return a*x**3 + b*x**2 + c*x + d

def _3rdOrderDeriv(x, a, b, c, d):
    return 3*a*x**2 + 2*b*x + c

def _1stOrder(x, a, b):
    return a*x + b

def _1stOrderDeriv(x, a, b):
    return a

def _expConst(x, a, b, c):
    return a * np.exp(-b * x) + c

posIndices = np.array([1, 3, 6, 8, 11, 14, 17, 20, 23, 25, 28])
posNum = 11
barNum = 1

read = False
if read:
    lower = 35
    upper = 67
    
    dataCalBef = []
    dataCalAft = []
    dataPos = []
    
    
    for i in range(upper - lower + 1):
        cache = DataLoad('D:\X' + str(lower + i)  + 'data.p')
        if (i in posIndices):
            dataPos.append(Bars(cache))
            
        if (i in posIndices - 1):
            dataCalBef.append(Bars(cache))
            
        if (i in posIndices + 1):
            dataCalAft.append(Bars(cache))
            
    '''
    pos0 = Bars(dataPos[1])
    pos1 = Bars(dataPos[3])
    pos2 = Bars(dataPos[6])
    pos3 = Bars(dataPos[8])
    pos4 = Bars(dataPos[11])
    pos5 = Bars(dataPos[14])
    pos6 = Bars(dataPos[17])
    pos7 = Bars(dataPos[20])
    pos8 = Bars(dataPos[23])
    pos9 = Bars(dataPos[25])
    pos10 = Bars(dataPos[28])
    bkg = Bars(dataPos[31])
    '''


calibrate = False
if calibrate:
    # Data protection
    arrCalBef = copy.deepcopy(dataCalBef)
    arrCalAft = copy.deepcopy(dataCalAft)
    arrPos = copy.deepcopy(dataPos)

    calFactorsBef = np.zeros((posNum, barNum))
    calFactorsAft = np.zeros((posNum, barNum))

    for i in range(posNum):
        calFactorsBef[i, :] = listEdgeCal(arrCalBef[i])[:barNum]
        calFactorsAft[i, :] = listEdgeCal(arrCalAft[i])[:barNum]
        
        for j in range(barNum):
            pi = arrPos[i][j][0, :]
            arrPos[i][j][0, :] = linearCalibrate(pi,
                                                 calFactorsBef[i, j],
                                                 calFactorsAft[i, j])

            

process = True
if process:
    pos = np.array([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50]) + 2.5 - 25
    ergBinEdges = np.arange(25, 1050, 100)
    ergCenters = (ergBinEdges[:-1] + ergBinEdges[1:]) / 2
    ergBinNum = len(ergCenters)
    
    binedges = np.arange(0, 1, 0.005)
    center = (binedges[:-1] + binedges[1:]) / 2
    times = np.array([144, 248, 288, 247, 206, 213, 192, 178, 277, 280, 201, 519/3]) * 300
    
    bkgNum = 10
        
    xLower = 0.2
    xUpper = 0.8
    
    maxes = np.zeros((barNum, posNum, ergBinNum))
    stdevs = np.zeros((barNum, posNum, ergBinNum))
    fwhms = np.zeros((barNum, posNum, ergBinNum))
    
    resolutions = np.zeros((barNum, ergBinNum))
    errResolutions = np.zeros((barNum, ergBinNum))
    
    fits = []   
    
    color = plt.cm.rainbow(np.linspace(0, 1, posNum))
    
    for k in range(ergBinNum):
        ergUpper = ergBinEdges[k + 1]
        ergLower = ergBinEdges[k]

    
        for i in range(barNum):
            thisBkg = arrPos[bkgNum][i]
            
            fig = plt.figure(figsize = (11, 6))
            ax = fig.add_subplot(111, projection='3d')
    
            ergBkgInds = np.logical_and(thisBkg[0] > ergLower,
                                        thisBkg[0] < ergUpper)
            
            bkgHist, temp = np.histogram(thisBkg[6, ergBkgInds], bins = binedges)
            bkgHistNormed = bkgHist / times[bkgNum]
      
            #for i in range
            for j in range(posNum - 1):
                col = color[j]
                
                thisPos = arrPos[j][i]
                
                    
                ergBinInds = np.logical_and(thisPos[0] > ergLower,
                                            thisPos[0] < ergUpper)
                
                hist, temp = np.histogram(thisPos[6, ergBinInds], bins = binedges)
                histNormedSubt = (hist / times[j]) - bkgHistNormed
                histNormedSubt[histNormedSubt < 0] = 0
                
    
                # Integral normalized
                # histNormedSubt /= np.sum(histNormedSubt)
                
                
                ax.bar(center, histNormedSubt,
                       zs = pos[j],
                       zdir = 'y',
                       width = binedges[1] - binedges[0] ,
                       label = str(pos[j]) + 'mm',
                       edgecolor = None,
                       #color = col,
                       alpha = 0.9)
                
                '''
                plt.bar(center, histNormedSubt,
                        width = binedges[1] - binedges[0],
                        label = str(pos[j])+ 'mm',
                        color = col, alpha = 0.3)
                '''
        
                popt, pcov = curve_fit(_gaussian,
                   center,
                   histNormedSubt,
                   p0 = [np.max(histNormedSubt),
                         center[np.argmax(histNormedSubt)], 0.1],
                   maxfev = 10000)
    
                xSmooth = np.arange(xLower, xUpper, 0.001)
                ax.plot(xSmooth,
                        _gaussian(xSmooth, *popt),
                        zs = pos[j],
                        linewidth = 1,
                        zdir = 'y',
                        c = 'k')
                
    
                maxes[i, j, k] = popt[1]
                stdevs[i, j, k] = np.abs(popt[2])
                #fwhms[i, j, k] = FWHM(center, histNormedSubt)
            
            
            
            ax.set_xlabel('Light collection ratio')
            ax.set_ylabel('Collimator position (mm)')
            ax.set_zlabel(r'Net count rate $(s)^{-1}$')
            ax.view_init(25, 45)
    
            #plt.xlabel('Bottom/Total')
            plt.xlim((xLower, xUpper))
            # plt.tight_layout()
            # plt.title(str(ergLower) + '-' + str(ergUpper) + ' keVee')
            plt.savefig( results_dir + 'Bar'+ str(i) + 'posHists.png')
            
            
            
            plt.figure()
    
            posToPlot = pos[:posNum - 1]
            maxesToPlot = maxes[i,:posNum - 1, k]
            errToPlot = np.array(stdevs[i])[:posNum - 1, k]
            #np.savetxt('positions.txt', )    
            
            plt.errorbar(posToPlot, maxesToPlot, xerr = 0.5, yerr = errToPlot, c = 'k',
                         ls='none', marker = '.', elinewidth = 1, capsize = 5)
            
            popt = np.polyfit(maxesToPlot, posToPlot, 1)
            x = np.linspace(0.2, 0.8, 1001)
            plt.plot(_1stOrder(x, *popt), x, c = 'k' )
            
    
            
            
            plt.xlabel('Height (mm)')
            plt.ylabel('Bottom/Total')
            plt.xlim((-30, 30))
            plt.tight_layout()
            plt.savefig(results_dir + 'Bar ' + str(i) + 'posFit.png')
            
            fits.append(popt)
            
            
            fig, ax = plt.subplots(figsize = (11, 6))
            
            plt.scatter(posToPlot,
                        errToPlot)
            plt.xlabel('Collimator position (mm)')
            plt.ylabel(r'$\sigma$ (ratio)')
            
            fig, ax = plt.subplots(figsize = (11, 6))         
            relErr = errToPlot * maxesToPlot
            
            resolutions[i, k] = np.mean(relErr) * _1stOrderDeriv(maxesToPlot, *popt)* -1
            errResolutions[i, k] = np.std(relErr) * _1stOrderDeriv(maxesToPlot, *popt)* -1
            #plt.scatter(posToPlot, np.mean(relErr) * _1stOrderDeriv(maxesToPlot, *popt)* -1)


plot = True
if plot:
    
    fig, ax = plt.subplots(figsize = (11, 6))
    
    color = plt.cm.rainbow(np.linspace(0, 1, 6))
    
    for i in range(barNum):
        maxesToPlot = maxes[i,:posNum - 1, 0]
        col = color[i]
        plt.errorbar(posToPlot, maxesToPlot, xerr = 0.5, yerr = errToPlot,
                     ls='none', marker = '.', elinewidth = 1, 
                     capsize = 5, label = 'Bar ' + str(i), c = col)
        
        plt.plot(x, _1stOrder(x, *fits[i]),
                 c = col)
        
    plt.xlabel('Beam height (mm)')
    plt.ylabel('Bottom/Total')
    plt.legend()
    plt.tight_layout()
    plt.savefig(results_dir + 'allBarsPosFit.png') 
        
    
    np.savetxt(results_dir + 'posFit.txt', fits)
    
    
    for i in range(barNum):
        fig, ax = plt.subplots(figsize = (11, 6))
        
        plotLinewidth = 1
        thisResolutions = resolutions[i, :]
        thisErr = errResolutions[i, :]
        
        plt.scatter(ergCenters,
                    thisResolutions,
                    s = 10,
                    c = 'k')
        
        
        plt.errorbar(ergCenters,
                     thisResolutions,
                     yerr = thisErr,
                     c = 'k',
                     elinewidth = plotLinewidth,
                     capsize = 2,
                     fmt = 'none')
    
        
        poptR, pcovR = curve_fit(_expConst,
           ergCenters,
           thisResolutions,
           p0 = [7, 0.015, 1.25],
           maxfev = 10000)
        
        xSmooth = np.linspace(ergBinEdges[0],
                              ergBinEdges[-1],
                              10000)
        
        plt.plot(xSmooth,
                 _expConst(xSmooth, *poptR),
                 'k--',
                 linewidth = plotLinewidth)
        
        
        plt.ylim(bottom = 0, top = np.floor(thisResolutions[0] + thisErr[0]) + 1)
        
        plt.xlabel('Light output (keVee)')
        plt.ylabel('Resolution (mm)')
        
        plt.savefig(results_dir + 'bar'+ str(i) + 'resPlot.png')

