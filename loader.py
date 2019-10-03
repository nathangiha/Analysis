# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 12:22:38 2019

@author: giha

Function for loading .p file produced by PreProcess.py package
"""

# Libraries
import numpy as np


def DataLoad(pFile):
    
    # Open file safely
    with open(pFile) as f:
        # Load channel count
        chCount = np.loadtxt( f,  comments='#', max_rows = 1)
        chCount = [chCount.astype(int)]
        chNum = len(chCount)
        
        # Load data
        data = np.loadtxt(f)

        # Separate data into channels using list
        a = []
        lineNum = int(0)
        for i in range(chNum):
            if chCount[i] ==0:
                a.append([])
            else:
                a.append(data[ lineNum:(lineNum + chCount[i]) ,:])
                lineNum = lineNum + chCount[i]
    return a




