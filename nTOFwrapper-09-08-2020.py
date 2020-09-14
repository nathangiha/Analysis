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

'''
X27 = DataLoad('D:\X27data.p')
X28 = DataLoad('D:\X28data.p')
X29 = DataLoad('D:\X29data.p')
X30 = DataLoad('D:\X30data.p')
X31 = DataLoad('D:\X31data.p')
X32 = DataLoad('D:\X32data.p')
'''

plt.close('all')
'''
X27 = DataLoad('D:\X27data.p')
X30 = DataLoad('D:\X30data.p')
X29 = DataLoad('D:\X29data.p')
X28 = DataLoad('D:\X28data.p')

print('Done reading')
'''

### Calibration and PSD for bars
'''
csCalData = Bars(X27)
tagCalData = X28
csCals = []
for i in range(len(csCalData)):
    csCals.append( EdgeCal( csCalData[i][0,:], histLabel = 'Bar '+str(i), xCal=0, integral = True ) )



nTOFData = Bars(X30)    
psdData = []
fitParams = []
for i in range(len(csCalData)):
    psdData.append( CalNClip(nTOFData[i][0,:], nTOFData[i][3,:], csCals[i][0] ) )
    fitParams.append( PSD_hist( psdData[i][0], psdData[i][1] , \
                               ergLimit=1000, discLine=3  )    \
                    )

### Calibration and PSD for tag detector    
tagCalData = X28[8]
plt.figure()
tagCal = EdgeCal( tagCalData[:,0], histLabel = 'Tag', xCal=0, integral=True)

tagData = X30[8]
tagPI = tagData[:,0]
tagTails = tagData[:,2]
tagTotals = tagData[:,3]
tagRatio = tagTails / tagTotals
tagPSD = CalNClip( tagPI, tagRatio, tagCal[0])
PSD_hist( tagPSD[0], tagPSD[1], ergLimit=1000, discLine=3  )
'''

results_dir = 'C:\\Users\\giha\\Documents\\Repos\\Analysis\\figures\\'

currentFile = X30
fname = 'X30'

flightTime = 1.356 #ns


if fname == 'X29':
    offset = []


barNum = 1
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
    
    neutronInds = 
    
    barNEvents = barAllEvents [                       ]
    
    
    
    
    barTime = barNEvents[4,:] + barNEvents[5,:]
    tagTime = tag[:,4] + tag[:,5]
    
    tagC = CoincFind(barTime, tagTime, 100)
    barInd = tagC[:,0]
    tagInd = tagC[:,1]
    
    barCoinc = barNEvents[:,barInd]
    
    botCoinc = bot[barInd,:]
    topCoinc = top[barInd,:]   
    tagCoinc = tag[tagInd, :]
    
    
    
    
    
    ttt = ( botCoinc[:,4] + topCoinc[:,5] ) / 2.0 - tagCoinc[:,4]
    cfd = ( botCoinc[:,4] + topCoinc[:,5] ) / 2.0 - tagCoinc[:,5]
    
    '''
    ttt = barCoinc[4,:] - tagCoinc[:,4]
    cfd = barCoinc[5,:] - tagCoinc[:,5]
    '''
    dt = ttt + cfd
    meanDt = np.mean(dt)
    
    
    
    
    dtp = dt[dt>0]
    plt.figure()
    
    if fname == 'X29':
        plt.hist(dt, bins=500, label = r'$\mu =$' + str(np.round(meanDt,5)))
    
        plt.xlim(2, 10)
        plt.title('Bar '+ str(i))
        plt.legend()
    
        offset.append(meanDt)
        
    if fname == 'X30':
        #plt.hist(dt - offset[i] + flightTime, bins=1000)
        plt.hist(dt , bins=500)
        plt.xlim(-30, 100)
        plt.title('Bar '+ str(i))

    plt.xlabel(r'$\Delta$t (ns)')
    plt.tight_layout()
    plt.savefig(results_dir + fname+ 'bar' + str(i)+'TOF.png')








