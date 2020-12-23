# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 11:00:24 2020

@author: giha

Functions for applying linear calibration over time, to account for gain shift

- Function for calculating calibration values at start and end of measurement
  based on calibration time (important for pos. measurements)
- Function for converting pulse integral from V-ns to keVee

"""

import numpy as np
from loader import DataLoad
#from psd import CalNClip, MovingAvg, PSD_hist, FOM_plot, PSD_ergslice, DiscLineN, DiscLineG
#from timeHist import TimeHist
import matplotlib.pyplot as plt
from coincFinder import CoincFind
from scipy.optimize import curve_fit

from bar import Bar, Bars
from edgeCal import EdgeCal

plt.close('all')

def listBarLinCal(datList, startCalList, endCalList):
    numCals = len(datList)
    calibratedData = []
    for i in range(numCals):
        calibratedData.append(linearCalibrate(datList[i][0,:],
                                              startCalList[i],
                                              endCalList[i]))
    return calibratedData
        

def linearCalibrate(pi, startCal, endCal):
    numPulses = len(pi)
    calibrationVec = np.linspace(startCal, endCal, numPulses)
    return pi * calibrationVec


'''
loadLinCal = False

if loadLinCal:
    X74 = DataLoad('D:\X74data.p')
    X75 = DataLoad('D:\X75data.p')
    X76 = DataLoad('D:\X76data.p')
    
    before = Bars(X74)
    resBars = Bars(X75)
    after  = Bars(X76)


calibrate = True
if calibrate:
    csCalsBeforebefore = []
    csCalsBefore = []
    csCalsAfter = []
    for i in range(len(before)):
        plt.figure()
        csCalsBeforebefore.append( EdgeCal( beforebefore[i][0,:], histLabel = 'bb Bar '+str(i), xCal=0, integral = True ) )
        csCalsBefore.append( EdgeCal( before[i][0,:], histLabel = 'b Bar '+str(i), xCal=0, integral = True ) )
        csCalsAfter.append(EdgeCal( after[i][0,:], histLabel = 'a Bar '+str(i), xCal=0, integral = True ) )
        # plt.savefig(results_dir + fname + 'csCals.png')


'''  
