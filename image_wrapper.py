# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 16:20:26 2020

@author: giha

A script for creating an image file from glass H2DPI calibrations
"""


import numpy as np
from loader import DataLoad
from bar import Bar, Bars
from edgeCal import EdgeCal
from psd import CalNClip, MovingAvg, PSD_hist, FOM_plot, PSD_ergslice, DiscLine
from timeHist import TimeHist
import matplotlib.pyplot as plt
from coincFinder import CoincFind
from scipy.optimize import curve_fit
from scipy import interpolate
from mpl_toolkits.mplot3d import Axes3D





##############################################################################
##############################################################################
##############################################################################
##############################################################################




def DiscLine(x, a, b,c,):
    return a*np.exp( -b*x) + c

# Define functions for using calibration data
def _InvPowerLaw(L, a, b):
    return (L / a)**(1/b)

def _InvZPos(ratio, p):
    heights = np.linspace(-25,75,1e4)
    ratios = p[0]*heights**3 + p[1]*heights**2 + p[2]*heights + p[3]
    findRatio = interpolate.interp1d(ratios, heights)
    return findRatio(ratio)

def _ZPos(height,p):
    return p[0]*height**3 + p[1]*height**2 + p[2]*height + p[3]

# Function for finding alpha
def GetLeverAlpha(d1,d2,tof,e1):
    # Switch order
    if tof < 0:
        temp = d1
        d1 = d2
        d2 = temp
    tail = d2
    # Lever arm
    vectof = (np.asarray(d1)-np.asarray(d2))
    
    # Calculate alpha
    dtof = np.linalg.norm(vectof) # cm
    tof = tof # ns
    '''
    if broaden == True:
        tof = tof + np.random.normal(0,timeRes)
    '''
    
    nmass = 1.008665 # amu
    
    etof = 0.5 * nmass * (dtof/tof)**2
    etof = etof * 1.03642697  # amu (cm/ns)^2 to MeV
    alpha = etof/ (etof + e1)
    
    
    return tail, vectof, alpha, etof, tof



##############################################################################
##############################################################################
##############################################################################
##############################################################################





sipm = 0.633

# Make dictionary for x,y positions of bars
det ={
0 : [-3*sipm, 0],
1 : [-sipm,0],
2 : [sipm,0],
3 : [3*sipm,0],
4 : [-2*sipm,-2*sipm],
5 : [2*sipm,-2*sipm]
}

outfile = 'results\\glassCones.p'



plt.close('all')

load = False
if load:
    # Load in data and calibrations
    lightOutputCal = np.loadtxt('figures\\lightOutputCurves.txt')
    zPosCal = np.loadtxt('figures\\posFit.txt')
    X69 = DataLoad('d:\X69data.p')
    barData = Bars(X69)


# Convert to bar data 

# Extract pulse times from bar data
barNum = len( barData )
eventNums = np.zeros( barNum )


# Calibrate and apply PSD cuts
loadCal = False
if loadCal:
    csCalData = Bars( DataLoad('D:\X68data.p') )

psdSigma = 4

calibrate = True
if calibrate:
    csCals = []
    for i in range(len(csCalData)):
        csCals.append( EdgeCal( csCalData[i][0,:], histLabel = 'Bar '+str(i), xCal=0, integral = True ) )
    
    
    
    psdData = []
    fitParams = []
    for i in range(len(barData)):
        print('PSD on bar ' + str(i))
        psdData.append( CalNClip(barData[i][0,:], barData[i][3,:], csCals[i][0] ) )
        fitParams.append( PSD_hist( psdData[i][0], psdData[i][1] , \
                                   ergLimit=1000, discLine=psdSigma  )  )
        plt.title('bar ' + str(i))


applyPSD = True
if applyPSD:
    nBarData = []

    for i in range(barNum):
        # Get PSD line fit params
        a = fitParams[i][0]
        b = fitParams[i][1]
        c = fitParams[i][2]
            
        # Get LO and tail/total
        barLO       = barData[i][0,:]*csCals[i][0]
        barPSDRatio = barData[i][3,:]

        # Find neutron indices
        neutronInds = np.nonzero(barPSDRatio > DiscLine(barLO, a, b, c))[0]
        
        # Construct nBarData
        nBarData.append( barData[i][:,neutronInds] )
        PSD_hist(nBarData[i][0,:]*csCals[i][0], nBarData[i][3,:])
        plt.title('Bar ' + str(i))
        plt.savefig('results\\bar' + str(i)+'psd.png')

applyEnergyThreshold = True
if applyEnergyThreshold:
    loThreshold = 100 # keVee
    for i in range(barNum):
        barLO = nBarData[i][0,:] * csCals[i][0]
        lowerCut = np.nonzero(barLO > loThreshold)[0]
        nBarData[i] = nBarData[i][:, lowerCut]


applyZCut = True
if applyZCut:
    for i in range(barNum):
        zMax = _ZPos(0, zPosCal[i])
        zMin = _ZPos(50, zPosCal[i])
        zRatios = nBarData[i][6,:]
        upperCut = np.nonzero(zRatios < zMax)[0]
        lowerCut = np.nonzero(zRatios > zMin)[0]
        passingInds = np.intersect1d(upperCut, lowerCut)
        nBarData[i] = nBarData[i][:,passingInds]


createTimeLists = True
if createTimeLists:
    timesList = [ ]
    zList = [ ]
    
    # Make lists of times and z-ratios
    for i in range( barNum ):
        eventNums[i] = np.size( nBarData[i], axis = 1 )
        timesList.append( (nBarData[i][4,:] + nBarData[i][5,:]))
        zList.append( (nBarData[i][6]))



'''
coincidenceTest = [CoincFind(timesList[b1], timesList[b2], 3) \
                for b1 in range(barNum) for b2 in range(b1 + 1, barNum) ]
'''

findCoinc = True
if findCoinc:
    numDoubles = 0
    coincidences = []
    channels = []
    for i in range(barNum):
        for j in range(i+1, barNum):
            print(i,j)
            channels.append((i,j))
            coince = CoincFind(timesList[i], timesList[j],10)
            numDoubles += np.size(coince, axis=0)
            coincidences.append( coince  )
            if coince.size == 0:
                print('Empty')

print('Event Count is ' + str(numDoubles))

# Calculate z-positions and light outputs
'''
calcZandE1 = True
if calcZandE1:
    zPos = []
    E1 = []
    for i in range(barNum):
        barRatios = nBarData[i][6,:]
        barZs = _InvZPos(barRatios, zPosCal[i])
        zPos.append(barZs)
        
        barLOs = nBarData[i][0,:] * csCals[i][0]
        BarE1s = _InvPowerLaw(barLOs, a = lightOutputCal[0,0], b = lightOutputCal[0,1])
        E1.append(barE1s)
'''


# Write cone parameters to script
                
cones = True

if cones:
    coneMat = np.zeros((numDoubles, 8))
    flightPaths = np.zeros((numDoubles,7))
    
    coneMat[:,7]= 0.005   
    eventsWritten = 0
    
    numPairs = len(coincidences)
    for i in range(numPairs):
        print('Analyzing pair ' + str(i))
        # Get number of coincidences in pair
        numPairDoubles = np.size(coincidences[i], axis=0)

        # Get x,y positions
        b1, b2 = channels[i]
        xy1 = det[b1]
        xy2 = det[b2]
        
        # Get coincidence indices for relevant bars
        b1Inds = coincidences[i][:,0]
        b2Inds = coincidences[i][:,1]
        
        # Get flight times
        b1CFDs = nBarData[b1][4,b1Inds]
        b1TTTs = nBarData[b1][5,b1Inds]
        b2CFDs = nBarData[b2][4,b2Inds]
        b2TTTs = nBarData[b2][5,b2Inds]
        
        cfds = b2CFDs - b1CFDs
        ttts = b2TTTs - b1TTTs
        
        dts = ttts + cfds
        
        # Get E1
        b1LOs = nBarData[b1][0,b1Inds] * csCals[b1][0]
        B1Es = _InvPowerLaw(b1LOs, a = lightOutputCal[0,0], b = lightOutputCal[0,1])
        # Get z's
        b1zratios = nBarData[b1][6,b1Inds]
        b2zratios = nBarData[b2][6,b2Inds]
                
        b1Zs = _InvZPos(b1zratios, zPosCal[b1]) / 10 # mm to cm
        b2Zs = _InvZPos(b2zratios, zPosCal[b2]) / 10 # mm to cm
        
        # Construct x,y,z interaction locations
        xyz1 = np.zeros((numPairDoubles,3))
        xyz2 = np.zeros((numPairDoubles,3))
        
        xyz1[:,:2] = xy1
        xyz2[:,:2] = xy2
        xyz1[:,2] = b1Zs
        xyz2[:,2] = b2Zs
        
        # Write positions to matrix
        coneMat[eventsWritten:eventsWritten+numPairDoubles,0:3] = xyz1
        coneMat[eventsWritten:eventsWritten+numPairDoubles,3:6] = xyz2
        
        '''
        print (dts[:10])
        plt.figure()
        plt.hist(dts,bins=100)
        '''

        for j in range(numPairDoubles):
            indexOverall = eventsWritten + j
            tof = dts[j]
            pos1 = xyz1[j]
            pos2 = xyz2[j]
            if tof < 0:
                temp = pos1
                pos1 = pos2
                pos2 = temp
            
            [tail, lever, alpha, Etof, TOF] = GetLeverAlpha(pos1, pos2, tof, B1Es[j]/1000 )
            coneMat[indexOverall, 6] = alpha
            flightPaths[indexOverall, :3] = tail
            flightPaths[indexOverall, 3:6] = lever
            flightPaths[indexOverall, 6] = alpha
            
        
        tempfile = 'results\\glassConesTemp.p'
        np.savetxt(tempfile, coneMat[eventsWritten:eventsWritten+j])
        SBP(tempfile)
        pair = str(b1) + '-' + str(b2) 
        plt.title(pair)
        plt.savefig('results\\' + pair + 'image.png' )           
            
        eventsWritten += numPairDoubles
        
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
flightPathsPlot = flightPaths[np.floor(np.linspace(0,numDoubles-1,500)).astype(int),:  ]

plt.quiver(flightPathsPlot[:,0], flightPathsPlot[:,1], flightPathsPlot[:,2], \
           flightPathsPlot[:,3], flightPathsPlot[:,4], flightPathsPlot[:,5], flightPathsPlot[:,6])
ax.set_xlim([-3, 3])
ax.set_ylim([-2, 1])
ax.set_zlim([-1, 6])
plt.savefig('results\\quiver.png')

      
write = True
if write:
    f = open(outfile, 'w')
    np.savetxt(f, coneMat)        
f.close()    

SBP(outfile)
plt.savefig('results\\image.png')











