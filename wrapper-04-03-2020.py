# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 11:08:52 2020

@author: giha
"""


import numpy as np
from loader import DataLoad
from bar import Bar
from edgeCal import EdgeCal
from psd import CalNClip, MovingAvg, PSD_hist, FOM_plot,PSD_ergslice
from timeHist import TimeHist
import matplotlib.pyplot as plt

'''
X27 = DataLoad('D:\X27data.p')
X28 = DataLoad('D:\X28data.p')
X29 = DataLoad('D:\X29data.p')
X30 = DataLoad('D:\X30data.p')
X31 = DataLoad('D:\X31data.p')
X32 = DataLoad('D:\X32data.p')
'''

plt.close('all')

b1b = X29[0]
b1t = X29[1]
tag = X29[8]


b1bt = b1b[:,4] + b1b[:,5]
b1tt = b1t[:,4] + b1t[:,5]


b1c = CoincFind(b1bt, b1tt, 5)
botInd = b1c[:,0]
topInd = b1c[:,1]

bot = b1b[botInd]
top = b1t[topInd]

B1 = Bar(bot,top)

B1t = B1[4,:] + B1[5,:]
Tagt = tag[:,4] + tag[:,5]

tagC = CoincFind(B1t, Tagt, 5)
barInd = tagC[:,0]
tagInd = tagC[:,1]

barCoinc = B1[:,barInd]
tagCoinc = tag[tagInd, :]
ttt = barCoinc[4,:] - tagCoinc[:,4]
cfd = barCoinc[5,:] - tagCoinc[:,5]

dt = ttt + cfd
dtp = dt[dt>0]
plt.hist(dt, bins=100)







