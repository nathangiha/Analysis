#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 10:33:12 2020

@author: giha
"""

from loader import DataLoad
from bar import Bar
from edgeCal import EdgeCal
from psd import CalNClip, MovingAvg, PSD_hist, FOM_plot,PSD_ergslice
from timeHist import TimeHist
import matplotlib.pyplot as plt

'''
X25 = DataLoad('/media/giha/DATA/X25data.p')
X26 = DataLoad('/media/giha/DATA/X26data.p')
'''
X25 = DataLoad('D:\X25data.p')
X26 = DataLoad('D:\X26data.p')


#B1cf = Bar(X25[0],X25[1])
B2cf = Bar(X25[2],X25[3])
B3cf = Bar(X25[4],X25[5])
B4cf = Bar(X25[6],X25[7])

B1cs = Bar(X26[0],X26[1])
B2cs = Bar(X26[2],X26[3])
B3cs = Bar(X26[4],X26[5])
B4cs = Bar(X26[6],X26[7])


plt.close('all')
#B1 = EdgeCal(B1cs[0,:], histLabel='B1', xCal=0, integral = True)
B2 = EdgeCal(B2cs[0,:], histLabel='B2', xCal=0, integral = True)
B3 = EdgeCal(B3cs[0,:], histLabel='B3', xCal=0, integral = True)
B4 = EdgeCal(B4cs[0,:], histLabel='B4', xCal=0, integral = True)



B2psd = CalNClip(B2cf[0,:],B2cf[3,:], B2[0] )
B3psd = CalNClip(B3cf[0,:],B3cf[3,:], B3[0] )
B4psd = CalNClip(B4cf[0,:],B4cf[3,:], B4[0] )

fig, ax = plt.subplots()



B4fom = FOM_plot(B4psd[0], B4psd[1], binsize=25 )

plt.savefig('/home/giha/Figures/fom.png',dpi =200 )
