# Silke Henkes 11.07.2013: Accelerated version of pebble game and rigid cluster code
#!/usr/bin/python

import random
import sys
import os
import glob
# Note: changed to simple pickle on move to python 3. Compatibility problems with loading old pickle files expected!
#import cPickle as pickle
import pickle as pickle
import copy as cp
import numpy as np
import scipy as sp
import matplotlib
import matplotlib.pyplot as plt

import Configuration as CF
import Pebbles as PB
import Hessian as HS
import Analysis as AN
import Tiling as TY

# ========================= Sample execution script. Will be removed in bulk version. ===========

#foldername = '/home/melle/Documents/Code/RigidLibrary/DataSimulation/conf1_N64_mu10_phi804/'
foldername = '/home/melle/Documents/Code/RigidLibrary/DataSimulation/conf1/'
form = 'simulation'
start = 420 #410 is partially rigid for the large system
stop = 2000
step = stop-start

outfolder = foldername

if not os.path.exists(outfolder):
    os.makedirs(outfolder)

#Hardcoded constants of the simulation
N = 1024
mu = 10

# Colors for the clusters. They are random, but I have redefined the first couple of ones to have preset colors, starting with black.
random_read = 0.8*np.random.rand(3*N, 3)
random_read[0, :] = [0, 0, 0]
random_read[1, :] = [0.8, 0, 0.8]
random_read[2, :] = [0, 0.8, 0]
random_read[3, :] = [0, 0, 0.8]
random_read[4, :] = [0.8, 0, 0]

#dummy= plt.figure(figsize=(15,7.5))
# dummy, (ax1, ax2) = plt.subplots(nrows=1,ncols1,figsize=(15,7.5))
clusout2 = np.zeros(int((stop-start)/step+1))
clusterall2 = [[] for _ in range(int((stop-start)/step+1))]

#Create config
ThisConf = CF.Configuration(foldername, form, mu, step)

for k in range(start, stop, step):
    #Read sim data
    ThisConf.readSimdata(k, False)
    
    ########## Setting up and playing the pebble game
    ThisPebble = PB.Pebbles(ThisConf,3,3,'nothing',False)
    
    # play game
    ThisPebble.play_game()
    # compute rigid clusters
    cidx, clusterall, clusterallBonds, clusteridx, BigCluster=ThisPebble.rigid_cluster()
    

    ########### Setting up the dynamical matrix and getting eigenmodes
    # This itself does very little, just creates an empty Hessian class
    # __init__(self,conf0):
    ThisHessian = HS.Hessian(ThisConf)
    ThisTiling = TY.Tiling(ThisConf)
    
    #Make the Maxwell-Cremona tiling
    ThisTiling.tile()
    
    ########## Have a look at some analysis functions of the rigid clusters
    #def __init__(self,conf0,pebbles0,hessian0,verbose=False):
    ThisAnalysis=AN.Analysis(ThisConf,ThisPebble,ThisHessian,ThisTiling,0.01,False)
    # stress statistics
    zav,nm,pres,fxbal,fybal,torbal,mobin,mohist,sxx,syy,sxy,syx=ThisAnalysis.getStressStat()
    # cluster statistics
    frac,fracmax,lenx,leny=ThisAnalysis.clusterStatistics()
    #def plotStresses(self,plotCir,plotVel,plotCon,plotF,plotStress,**kwargs):
    #fig1 = ThisAnalysis.plotStresses(True,False,False,True,False)
    #def plotPebbles(self,plotCir,plotPeb,plotPebCon,plotClus,plotOver,**kwargs):
    #ThisAnalysis.plotPebbles(True,True,True,False,False)
    #fig2 = ThisAnalysis.plotPebbles(True,False,False,True,False)
    
    #Plotting the contact network
    #fig3 = ThisAnalysis.contactnetwork()
    
    #Plotting Maxwell-Cremona tiling
    #Colorscheme options filled = False: cluster, force, colorblind, random
    #Colorscheme options filled = True: colorblind, random
    
    #fig5 = ThisAnalysis.tileplotter(colorscheme='cluster', filled=False)
    #fig6 = ThisAnalysis.tileplotter(colorscheme='random', filled=True)
    
    ThisAnalysis.tiling_statistics()
    """
    fig1.set_size_inches(15,15)
    fig1.savefig('data_1024_rigid.pdf', dpi=100)
    
    fig2.set_size_inches(15,15)
    fig2.savefig('pebble_1024_floppy.pdf', dpi=100)
    
    fig5.set_size_inches(15,15)
    fig5.savefig('tiles_1024_rigid.pdf', dpi=100)
    
    fig6.set_size_inches(15,15)
    fig6.savefig('tiles_1024_rigid_alt.pdf', dpi=100)
    """
    
    #plt.show()