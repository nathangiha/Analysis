# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 12:14:17 2020

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
X33 = DataLoad('D:\X33data.p')
X34 = DataLoad('D:\X33data.p')
'''

plt.close('all')



B1cs1 = Bar(X27[0],X27[1])
B2cs1 = Bar(X27[2],X27[3])
B3cs1 = Bar(X27[4],X27[5])
B4cs1 = Bar(X27[6],X27[7])

B1cs2 = Bar(X33[0],X33[1])
B2cs2 = Bar(X33[2],X33[3])
B3cs2 = Bar(X33[4],X33[5])
B4cs2 = Bar(X33[6],X33[7])
B5cs2 = Bar(X33[8],X33[9])
B6cs2 = Bar(X33[10],X33[11])

B1i = EdgeCal(B1cs1[0,:], histLabel='B1i', xCal=0, integral = True)
B2i = EdgeCal(B2cs1[0,:], histLabel='B2i', xCal=0, integral = True)
B3i = EdgeCal(B3cs1[0,:], histLabel='B3i', xCal=0, integral = True)
B4i = EdgeCal(B4cs1[0,:], histLabel='B4i', xCal=0, integral = True)



B1f = EdgeCal(B1cs2[0,:], histLabel='B1f', xCal=0, integral = True)
B2f = EdgeCal(B2cs2[0,:], histLabel='B2f', xCal=0, integral = True)
B3f = EdgeCal(B3cs2[0,:], histLabel='B3f', xCal=0, integral = True)
B4f = EdgeCal(B4cs2[0,:], histLabel='B4f', xCal=0, integral = True)
B5f = EdgeCal(B5cs2[0,:], histLabel='B5f', xCal=0, integral = True)
B6f = EdgeCal(B6cs2[0,:], histLabel='B6f', xCal=0, integral = True)



plt.tight_layout()

'''
plt.figure()
plt.xlabel('Bar #')
plt.ylabel('keVee/V-ns')
plt.scatter(np.arange(1,7), stil, label = 'Stilbene')
plt.scatter(np.arange(1,5), glass, label='Glass')
plt.axhline(y=savg, label= 'Stilbene avg',color = 'blue')
plt.axhline(y=gavg, label = 'Glass avg', color = 'orange')


plt.legend()
plt.tight_layout()
'''

ini = np.array(  [B1i[0], B2i[0], B3i[0], B4i[0] ] )
fin = np.array( [ B1f[0], B2f[0], B3f[0], B4f[0], B5f[0], B6f[0] ])
plt.figure()
plt.xlabel('Bar #')
plt.ylabel('Compton Edge Calibration (keVee/V-ns)')
plt.scatter(np.arange(1,5), ini, label = 'March', c='k')
plt.scatter(np.arange(1,7), fin, label='July', c='b')
plt.axhline(y= np.mean(ini), label= 'March avg',color = 'black')
plt.axhline(y= np.mean(fin), label = 'July avg', color = 'blue')
plt.legend()
plt.tight_layout()



B1cf = Bar(X34[2],X34[3])
B2cf = Bar(X34[2],X34[3])
B3cf = Bar(X34[4],X34[5])
B4cf = Bar(X34[6],X34[7])
B5cf = Bar(X34[8],X34[9])
B6cf = Bar(X34[10],X34[11])









