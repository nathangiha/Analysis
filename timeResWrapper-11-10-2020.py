# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 10:42:40 2020

@author: giha

Examine time resolution of glass bars, perhaps performing energy cuts as well
"""


from loader import DataLoad
from bar import Bar
from edgeCal import EdgeCal
from timeHist import doTimeHist
import matplotlib.pyplot as plt
import numpy as np
from coincFinder import CoincFind



plt.close('all')


resultPath = 'C:\\Users\\giha\\Documents\\Papers\\glassBars\\'


loadData = False
if loadData:
    X23 = DataLoad('D:\\X23data.p')
    X24 = DataLoad('D:\\X24data.p')

doBars = True
if doBars:    
    b1cs = Bar(X23[0],X23[1])
    b2cs = Bar(X23[2],X23[3])

    b1na = Bar(X24[0],X24[1])
    b2na = Bar(X24[2],X24[3])


calibrate = True
if calibrate:

    b1Cal = EdgeCal(b1cs[0, :], histLabel = 'b1', xCal = 0, integral = True)
    b2Cal = EdgeCal(b2cs[0, :], histLabel = 'b2', xCal = 0, integral = True)

    b1na[0, :] *= b1Cal[0]
    b2na[0, :] *= b2Cal[0]

ergCuts = True
if ergCuts:
    
    # Apply 50 keVee threshold
    totalThresh = 50
    b1totInds = b1na[0, :] > totalThresh
    b2totInds = b2na[0, :] > totalThresh
    b1naThr = b1na[:, b1totInds]
    b2naThr = b2na[:, b2totInds]


    # Gate on Compton edge
    lowerThresh = 290
    upperThresh = 390
    
    b1ThresholdInds = np.logical_and(b1na[0, :] > lowerThresh,
                                          b1na[0, :] < upperThresh)
    b2ThresholdInds = np.logical_and(b2na[0, :] > lowerThresh,
                                          b2na[0, :] < upperThresh)
    
    b1naCut = b1na[:, b1ThresholdInds]
    b2naCut = b2na[:, b2ThresholdInds]
    
    fig, ax = plt.subplots()
    binEdges = np.arange(0, 600, 10)
    b1Spec = plt.hist(b1na[0, :], bins = binEdges,
                      histtype = 'step',
                      label = 'Bar 1')
    b2Spec = plt.hist(b2na[0, :], bins = binEdges,
                      histtype = 'step',
                      label = 'Bar 2')

    b1SpecCut = plt.hist(b1naCut[0, :], bins = binEdges,
                         alpha = 0.5,
                         label = 'Cut 1')
    b2SpecCut = plt.hist(b2naCut[0, :], bins = binEdges,
                         alpha = 0.5,
                         label = 'Cut 2')
    
    ax.set_xlabel('Light collected (keVee)')
    ax.set_ylabel(r'Counts/$\Delta$keVee')
    ax.legend()
    fig.tight_layout()
    
    
    
timeHist = True
if timeHist:
    # 50 keVee threshold
    b1Time = b1naThr[4, :] + b1naThr[5, :]
    b2Time = b2naThr[4, :] + b2naThr[5, :]
    
    barCoinc = CoincFind(b1Time, b2Time, 5) # 5ns window
    b1Ind = barCoinc[: ,0]
    b2Ind = barCoinc[:, 1]

    b1Coinc = b1naThr[:, b1Ind]
    b2Coinc = b2naThr[:, b2Ind]

    cfd = b1Coinc[4, :] - b2Coinc[4, :]
    ttt = b1Coinc[5, :] - b2Coinc[5 ,:]
    tdif = cfd + ttt
    
    doTimeHist(tdif)
    plt.savefig(resultPath + 'timeResGlass50keVThreshold.png',
                dpi = 300)
    
    # With energy cuts
    b1Time = b1naCut[4, :] + b1naCut[5, :]
    b2Time = b2naCut[4, :] + b2naCut[5, :]
    
    barCoinc = CoincFind(b1Time, b2Time, 5) # 5ns window
    b1Ind = barCoinc[: ,0]
    b2Ind = barCoinc[:, 1]

    b1Coinc = b1naCut[:, b1Ind]
    b2Coinc = b2naCut[:, b2Ind]



    cfd = b1Coinc[4, :] - b2Coinc[4, :]
    ttt = b1Coinc[5, :] - b2Coinc[5 ,:]
    tdif = cfd + ttt
    
    doTimeHist(tdif)
    plt.savefig(resultPath + 'timeResGlass290-390keVee.png',
            dpi = 300)

