#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 13:37:04 2019

@author: giha
"""

from loader import DataLoad
from bar import Bar
from edgeCal import EdgeCal
from psd import CalNClip, MovingAvg, PSD_hist, FOM_plot,PSD_ergslice
from timeHist import TimeHist
import matplotlib.pyplot as plt

'''
X23 = DataLoad('/media/giha/DATA/X23data.p')
X24 = DataLoad('/media/giha/DATA/X24data.p')
'''
B1cs = Bar(X23[0],X23[1])
B2cs = Bar(X23[2],X23[3])

B1na = Bar(X24[0],X24[1])
B2na = Bar(X24[2],X24[3])

plt.close('all')
B1 = EdgeCal(B1cs[0,:], histLabel='B1', xCal=0, integral = True)
B2 = EdgeCal(B2cs[0,:], histLabel='B2', xCal=0, integral = True)




cfd = B1na[4,:] - B2na[4,:]
ttt = B1na[5,:] - B2na[5,:]
tdif = cfd + ttt

test = TimeHist(tdif)