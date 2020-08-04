# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 10:17:05 2020

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
B1cs1 = Bar(X27[0],X27[1])
B2cs1 = Bar(X27[2],X27[3])
B3cs1 = Bar(X27[4],X27[5])
B4cs1 = Bar(X27[6],X27[7])
Scs1 = X28[8]

B1cs2 = Bar(X31[0],X31[1])
B2cs2 = Bar(X31[2],X31[3])
B3cs2 = Bar(X31[4],X31[5])
B4cs2 = Bar(X31[6],X31[7])
Scs2 = X32[8]


B1i = EdgeCal(B1cs1[0,:], histLabel='B1i', xCal=0, integral = True)
B2i = EdgeCal(B2cs1[0,:], histLabel='B2i', xCal=0, integral = True)
B3i = EdgeCal(B3cs1[0,:], histLabel='B3i', xCal=0, integral = True)
B4i = EdgeCal(B4cs1[0,:], histLabel='B4i', xCal=0, integral = True)


'''
B1f = EdgeCal(B1cs2[0,:], histLabel='B1f', xCal=0, integral = True)
B2f = EdgeCal(B2cs2[0,:], histLabel='B2f', xCal=0, integral = True)
B3f = EdgeCal(B3cs2[0,:], histLabel='B3f', xCal=0, integral = True)
B4f = EdgeCal(B4cs2[0,:], histLabel='B4f', xCal=0, integral = True)
'''
'''
B1i = EdgeCal(B1cs1[0,:], histLabel='B1i', xCal=0, integral = True)
B1f = EdgeCal(B1cs2[0,:], histLabel='B1f', xCal=0, integral = True)
'''

print( B1i[0], B2i[0], B3i[0], B4i[0], B1f[0], B2f[0], B3f[0], B4f[0])


plt.tight_layout()


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
'''
