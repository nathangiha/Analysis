# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 14:31:18 2019

@author: giha

Do PSD
"""
import matplotlib.pyplot as plt
import numpy as np

data = X14

ph = data[0][:,1]

clipped = np.where(ph < 8000)[0]

tail = data[0][clipped,2]
total = data[0][clipped,3]



ratio = np.divide(tail,total)

binnum=200
plt.close()
psd = plt.hist2d(ph[clipped],ratio, bins=(binnum,500*binnum ), cmap=plt.cm.jet)
plt.colorbar()
plt.xlim(0, 8191)
plt.ylim(0, 0.4) 