# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 11:46:38 2019

@author: giha

Convert 2-channel bar data to single channel output
"""

import numpy as np
from coincFinder import CoincFind





# This will be a function when I'm done, but I have it like this for easy debugging

def Bar(bottom, top):
#bottom = X18[0]
#top = X18[1]

    barBotTime = bottom[:,4] + bottom[:,5]
    barTopTime = top[:,4] + top[:,5]


    barCoinc = CoincFind(barBotTime, barTopTime, 5)
    botInd = barCoinc[:,0]
    topInd = barCoinc[:,1]

    bottom = bottom[botInd]
    top = top[topInd]

    
    
    ph = bottom[:,1] + top[:,1]
    pi = bottom[:,0] + top[:,0]
    zratio = bottom[:,0] / pi
    cfd = (bottom[:,4] + top[:,4]) /2
    ttt = (bottom[:,5] + top[:,5]) /2
    # tstart = cfd + ttt
    #tails = bottom[:,2] + top[:,2]
    #totals = bottom[:,3] + top[:,3]
    tails = np.sqrt( bottom[:,2]**2 + top[:,2]**2  )
    totals = np.sqrt( bottom[:,3]**2 + top[:,3]**2  )
    PSDratio = tails/totals
    
    bar = np.vstack((pi, tails, totals, PSDratio, cfd, ttt, zratio))
    
    
    return bar


def Bars(filedata):
    print('Running Bars...')
    barNum = int(np.floor(  len(filedata )/2 ))
    bardata = []
    for i in range(barNum):
        temp = Bar( filedata[2*i], filedata[2*i+1]     )
        bardata.append(temp)
        
    return bardata






























    