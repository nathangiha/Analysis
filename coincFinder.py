# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 11:22:45 2020

@author: giha

Function for finding coincidences within a specified window (ns)

Inputs: Time vec1 of length M, time vec2 of length N, window size
Output: Vectors of indices corresponding to coinc for vec1, vec 2 
"""
import numpy as np

def CoincFind(v1, v2, window):
    
    coincInd = np.full( (max( [len(v1),len(v2)] ), 2), np.nan)

    i = 0
    for n in range(len(v1)):
        while( v1[n] + window > v2[i]):
            if( v1[n] - window < v2[i] and v1[n] + window > v2[i]):
                coincInd[n,:] = [n,i]
            if i < len(v2)-1:
                i+=1
            else:
                break
    
    coincInd = coincInd[~np.isnan(coincInd).any(axis=1)].astype(int)
    return coincInd
    
    
    
    
    
