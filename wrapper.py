#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 11:00:35 2019

@author: giha
"""

from loader import DataLoad
from bar import Bar
from edgeCal import EdgeCal
from psd import CalNClip, MovingAvg, PSD_hist, FOM_plot,PSD_ergslice
import matplotlib.pyplot as plt
'''
#X17 = DataLoad('D:\\X17data.p')
X18 = DataLoad('D:\\X18data-0-20-370.p')
#X19 = DataLoad('D:\\X19data.p')
'''
'''
X20 = DataLoad('/media/giha/DATA/X20data.p')
X21 = DataLoad('/media/giha/DATA/X21data.p')
X22 = DataLoad('/media/giha/DATA/X22data.p')
'''
# Perform bar analysis
B1cs1 = Bar(X20[0],X20[1])
B2cs1 = Bar(X20[2],X20[3])

B1cf = Bar(X21[0],X21[1])
B2cf = Bar(X21[2],X21[3])

B1cs2 = Bar(X22[0],X22[1])
B2cs2 = Bar(X22[2],X22[3])

plt.close('all')
# Energy calibration
B11 = EdgeCal(B1cs1[0,:], histLabel='B11', xCal=0, integral = True)
B12 = EdgeCal(B1cs2[0,:], histLabel='B12', xCal=0, integral = True)
B21 = EdgeCal(B2cs1[0,:], histLabel='B21', xCal=0, integral = True)
B22 = EdgeCal(B2cs2[0,:], histLabel='B22', xCal=0, integral = True)

plt.savefig('/home/giha/Figures/CsCal.png',dpi=200)
# Evaluate PSD
B1psd = CalNClip(B1cf[0,:],B1cf[3,:], (B11[0]+B12[0])/2 )
B1fom = FOM_plot(B1psd[0], B1psd[1], binsize=25 )
plt.savefig('/home/giha/Figures/B1psd.png',dpi =200 )


B2psd = CalNClip(B2cf[0,:],B2cf[3,:], (B21[0]+B22[0])/2 )
B2fom = FOM_plot(B2psd[0], B2psd[1], binsize=25 )
plt.savefig('/home/giha/Figures/B2psd.png',dpi =200 )
