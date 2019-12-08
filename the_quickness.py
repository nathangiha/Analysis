# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 17:13:16 2019

@author: giha

Short script that reads my data in for debugging my bar PSD stuff
"""
from loader import DataLoad
from bar import Bar
from edgeCal import EdgeCal
'''
#X17 = DataLoad('D:\\X17data.p')
X18 = DataLoad('D:\\X18data-0-20-370.p')
#X19 = DataLoad('D:\\X19data.p')
'''
'''
X17 = DataLoad('/media/giha/DATA/X17data.p')
X18 = DataLoad('/media/giha/DATA/X18data-0-20-370.p')
X19 = DataLoad('/media/giha/DATA/X19data.p')
'''
X18 = DataLoad('/media/giha/DATA/X18data-0-40-370.p')

B1cs1 = Bar(X17[0],X17[1])
B2cs1 = Bar(X17[2],X17[3])

B1cf = Bar(X18[0],X18[1])
B2cf = Bar(X18[2],X18[3])

B1cs2 = Bar(X19[0],X19[1])
B2cs2 = Bar(X19[2],X19[3])

B11 = EdgeCal(B1cs1[0,:], histLabel='B11', xCal=0, integral = True)
B12 = EdgeCal(B1cs2[0,:], histLabel='B12', xCal=0, integral = True)
B21 = EdgeCal(B2cs1[0,:], histLabel='B21', xCal=0, integral = True)
B22 = EdgeCal(B2cs2[0,:], histLabel='B22', xCal=0, integral = True)
