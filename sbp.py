# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 01:11:09 2020

@author: giha
"""


#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as col
import matplotlib.patches as patches
from matplotlib import ticker, cm
import numpy.ma as ma

import os.path
import glob
import time



################################################################
############### user-specified variables #######################
################################################################
#
# Initialize IO path


# Source location
srcloc = (-0.378, 2.98)
sw = 20
sh= 20



# Specify bin numbers
azBin = 360
altBin = 180

plt.close('all')


################################################################
############### user-specified variables #######################
################################################################

#def srcCircleInd(srcloc, sw, sh):
alt, az = np.ogrid[-89.5:90.5,-179.5:180.5]
dist_from_center = np.sqrt((az - srcloc[0])**2 + (alt - srcloc[1])**2)
mask = dist_from_center <= sw/2
    
   # return mask
    
    


def sph2cart(r, theta, phi):
    return np.array([r * np.sin(theta)*np.cos(phi),
         r * np.cos(theta)*np.cos(phi),
         r * np.sin(phi)
    ])

        

#######
def SBPPlot(mat):
    # Plot ring
    # plt.close()
    
    
    fig, ax = plt.subplots()
    
    cs = plt.imshow(mat, cmap = 'jet',  aspect = 'auto', origin = 'lower')
    norm= col.LogNorm(vmin= np.max(mat)*1e-12, vmax=np.max(mat))
    
    #norm= col.LogNorm(vmin=1e-3, vmax=cs.cvalues.max())
    sm = plt.cm.ScalarMappable(norm=norm, cmap = cs.cmap)
    sm.set_array([])
    #fig.colorbar(sm, ticks=np.linspace(0,1.1,12))
    fig.colorbar(sm)
    
    
    #fig.tight_layout()
    plt.xlabel(r'Azimuth ($\theta$)')
    plt.ylabel(r'Altitude ($\phi$)')
    xlabels = np.round(np.linspace(-180,180, 9)).astype(int)
    ylabels = np.round(np.linspace(-90, 90, 7)).astype(int)
    plt.xticks(xlabels+179.5, xlabels)
    plt.yticks(ylabels + 89.5, ylabels)
    
    
    fig.set_size_inches(9,4)
    plt.tight_layout()

    plt.show()
    
    
    # Create a Rectangle patch
    ell = patches.Ellipse((179.5 + srcloc[0], 89.5 + srcloc[1]),sw, sh,linewidth=1,edgecolor='w',facecolor='none')

    # Add the patch to the Axes
    ax.add_patch(ell)
    plt.tight_layout()
    
    return
    



def SBP(conefile):
    
    # Construct array of bin vectors
    bVecs = np.zeros([3,altBin,azBin])
    
    # Bin Edges
    azEdge = np.linspace(-np.pi, np.pi, azBin+1)
    altEdge = np.linspace(-np.pi/2, np.pi/2, altBin+1)
    
    # Convert to centers
    azC = (azEdge[1:] + azEdge[:-1] ) / 2
    altC = (altEdge[1:] + altEdge[:-1] ) / 2
    
    # Convert to Cartesian for computational ease
    for i in range(altBin):
        for j in range(azBin):
            bVecs[:,i,j] = sph2cart(30,azC[j],altC[i])
    
    simdata = np.loadtxt(conefile)
    if not simdata.size:
        return

    try:
        v1 = simdata[:, 0:3]
        v2 = simdata[:, 3:6]
        alphas = simdata[:,6]
        varAlphas = simdata[:,7]
    except:
        print('Threw out pair with 1 cone')
        return        
    del simdata
    
    sigB = 0.001
    #varA = 0.05
    
    # Find beta for each bin
    bMat = np.zeros([altBin,azBin])
    cMat = np.zeros([altBin,azBin])
    totMat = np.zeros([altBin,azBin])
    
    
    # pointervec = [0,1,0]
    
    histnum = len(v1)
    #histnum=4
    for k in range(histnum):
        #bMat = np.empty([azBin,altBin])
        #cMat = np.empty([azBin,altBin])
        alpha = alphas[k]
        lever =  (v1[k,:] - v2[k,:]) / np.linalg.norm(v1[k,:] - v2[k,:])
        v1Vecs = (np.ones((3,altBin,azBin)).T * v1[k,:].T).T
        
        aVecs = bVecs - v1Vecs
        aVecsNormed = aVecs / np.linalg.norm(aVecs, axis=0)
        
        dots = np.sum(aVecsNormed.T * lever, axis=2).T
        
        bMat = np.square(dots)
        
        # Get variance from file
        varA = varAlphas[k]
        
        cMat = np.exp(-(bMat - alpha)**2  / (2*(sigB**2 + varA) ) )
        cMat[ dots<0] = 0
        
        '''
        for i in range(altBin):
            for j in range(azBin):
                aVec = ( bVecs[:,i,j]-v1[k,:] )
                aVecNormed = aVec / np.linalg.norm(aVec)
                dot = np.dot(lever, aVecNormed)        
                bMat[i][j] = dot**2
                # dot2 = np.dot(lever, pointervec)
                if dot > 0:            
                    cMat[i][j] = np.exp(-(bMat[i][j]-alpha)**2 / (2*(sigB**2 + sigA**2) ) ) 
                else:
                    cMat[i][j] = 0
        '''                
        # Normalize if nonzero
        if not np.all(cMat==0):           
            cMat = cMat / np.sum(cMat)    
            
        totMat = totMat + cMat
        # time.sleep(.1)
        
    
    
    #totMat = ma.array(totMat, mask = totMat==0)
    # Define hotspot and background images based on mask
        
        
    SBPPlot(totMat)
    '''
    # Find rotation angle from filename
    inpfile = os.path.splitext(   conefile     )
    phi = inpfile[0][-2:]
    ext = inpfile[1][1:]
    
    plt.title(str(histnum) + " cones " + phi)
    plt.tight_layout()

    plt.savefig(ext+'alt'+phi ,dpi=200 )
    
    '''
    
    #plt.close()
    return totMat#, phi


# Takes image matrix and plots hotspot/background, outputs hotspot fraction and background matrix
def ImageQuality(mat):
    hotspot = ma.masked_array(mat, mask = ~mask)
    SBPPlot(hotspot)
    
    background = ma.masked_array(mat, mask = mask)
    SBPPlot(background)
    
    fracHotspot = np.sum(hotspot)/np.sum(mat)
    
    return fracHotspot, background

'''
    
conefiles = glob.glob('*.pt')
#conefiles.extend(glob.glob('*.pt'))

#SBP(conefiles[-1])

summedMat =  np.zeros([altBin,azBin])

quals = []
angles = []

for f in conefiles[:]:
    mat, angle = SBP(f)
    summedMat = summedMat + mat
    
    angles.append(angle)
    quals.append(ImageQuality(mat)) 
    #time.sleep(.1)

angles = np.array(angles)
indangles = np.argsort(angles)
angles.sort()
quals = np.array(quals)[indangles]
fig, ax = plt.subplots()
plt.scatter(angles, quals[:,0]*100, s=10, c='k')
plt.xlabel(r'Tilt ($\phi$)')
plt.ylabel('Hotspot fraction')
plt.xlim(-0.1,6.1)
plt.tight_layout()
plt.savefig('hotspotfrac.png')



'''








    
    #plt.savefig(ext+'alt'+phi ,dpi=200 )
'''
summedMirror = np.flip(summedMat, axis=0)
SBPPlot(summedMat)
SBPPlot(summedMirror)
SBPPlot(summedMat+summedMirror)
'''