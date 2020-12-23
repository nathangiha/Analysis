# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 17:05:38 2020

@author: giha
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from loader import DataLoad
from bar import Bar, Bars
from edgeCal import EdgeCal
from psd import CalNClip, MovingAvg, PSD_hist, FOM_plot,PSD_ergslice
from timeHist import TimeHist

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

read = False
if read:
    lower = 35
    upper = 67
    
    
    dataPos = []
    
    
    for i in range(upper - lower + 1):
        cache =  DataLoad('D:\X' + str(lower + i)  + 'data.p')
        dataPos.append(cache)

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

pos0 = Bars(DataLoad('D:\X36data.p'))
pos1 = Bars(DataLoad('D:\X38data.p'))
pos2 = Bars(DataLoad('D:\X41data.p'))
pos3 = Bars(DataLoad('D:\X43data.p'))
pos4 = Bars(DataLoad('D:\X46data.p'))
pos5 = Bars(DataLoad('D:\X49data.p'))
pos6 = Bars(DataLoad('D:\X52data.p'))
pos7 = Bars(DataLoad('D:\X55data.p'))
pos8 = Bars(DataLoad('D:\X58data.p'))
pos9 = Bars(DataLoad('D:\X60data.p'))
pos10 = Bars(DataLoad('D:\X63data.p'))
bkg = Bars(DataLoad('D:\X66data.p'))
'''
allpos = [pos0, pos1, pos2, pos3, pos4, pos5, pos6, pos7, pos8, pos9, pos10, bkg]
pos = np.array([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, -1]) + 2.5
binedges = np.linspace(0, 1, 201)
center = (binedges[:-1] + binedges[1:]) / 2
times = [144, 248, 288, 247, 206, 213, 192, 178, 277, 280, 201, 519/3]

barNum = len(pos0)

maxes = np.zeros(( barNum, len(allpos) ))
stdevs = np.zeros(( barNum, len(allpos) ))
fwhms = np.zeros(( barNum, len(allpos) ))
fits = []



color = plt.cm.rainbow(np.linspace(0, 1, 10))

for i in range( barNum ):
    plt.figure()
    bkgHist, temp = np.histogram(pos10[i][6], bins = binedges)
    bkgHistNormed = bkgHist / times[10]

    #for i in range
    for j in range(10):
        col = color[j]

        # plt.hist(allpos[j][2][6], bins=500, density = True, label = str(j*5)+ 'mm' , alpha = 0.5 )
        hist, temp = np.histogram(allpos[j][i][6], bins = binedges)
        histNormedSubt = hist / times[j] - bkgHistNormed
        histNormedSubt[histNormedSubt < 0] = 0
        histNormedSubt /= np.max(histNormedSubt)
        plt.bar(center, histNormedSubt,
                width = binedges[1] - binedges[0],
                label = str(pos[j])+ 'mm',
                color = col, alpha = 0.3)
        
        popt, pcov = curve_fit(_gaussian,
                       center,
                       histNormedSubt,
                       p0 = [np.max(histNormedSubt),
                             center[np.argmax(histNormedSubt)], 0.1],
                       maxfev = 10000)

        
        xSmooth = np.arange(0, 1, 0.001)
        plt.plot(xSmooth, _gaussian(xSmooth, *popt), c = col)
        
        maxes[i, j] = popt[1]
        stdevs[i, j] = popt[2]
        fwhms[i, j] = (FWHM (center, histNormedSubt))
        
        
        
        
    plt.xlabel('Bottom/Total')
    plt.xlim((0.2, 0.8))
    plt.legend()
    plt.tight_layout()
    plt.savefig( results_dir + 'Bar'+ str(i) + 'posHists.png')
    
    
    
    plt.figure()
    posToPlot = pos[:-2]
    maxesToPlot = maxes[i,:-2]
    errToPlot = np.array(stdevs[i] / 2.355)[:-2]
    #np.savetxt('positions.txt', )    
    
    plt.errorbar(posToPlot, maxesToPlot, xerr = 1, yerr = errToPlot, c = 'k',
                 ls='none', marker = '.', elinewidth = 1, capsize = 5)
    
    p = np.polyfit(posToPlot, maxesToPlot, 3)
    x = np.linspace(0,50,1001)
    plt.plot(x, _3rdOrder(x, *p), c = 'k' )
    
    
    plt.xlabel('Height (mm)')
    plt.ylabel('Bottom/Total')
    plt.xlim((-5,55))
    plt.tight_layout()
    plt.savefig(results_dir + 'Bar ' + str(i) + 'posFit.png')
    
    fits.append(p)


fig, ax = plt.subplots(figsize = (11, 6))

color = plt.cm.rainbow(np.linspace(0, 1, 6))

for i in range(6):
    maxesToPlot = maxes[i,:-2]
    col = color[i]
    plt.errorbar(posToPlot, maxesToPlot, xerr = 1, yerr = errToPlot,
                 ls='none', marker = '.', elinewidth = 1, 
                 capsize = 5, label = 'Bar ' + str(i), c = col)
    
    plt.plot(x, _3rdOrder(x, *fits[i]),
             c = col)
    
plt.xlabel('Beam height (mm)')
plt.ylabel('Bottom/Total')
plt.legend()
plt.tight_layout()
plt.savefig(results_dir + 'allBarsPosFit.png') 
    
    
    

np.savetxt(results_dir + 'posFit.txt', fits)