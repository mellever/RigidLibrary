import random
import sys, os, glob


# Location of the library files: insert into python path
pebblepath='/home/sh18581/Documents/Friction/RigidLibrary/'
sys.path.insert(1,pebblepath)

# Load pebble game and dynamical matrix libraries
import Configuration as CF
import Pebbles as PB
import Hessian as HS
import Analysis as AN

import matplotlib.pyplot as plt
import numpy as np

#topdir='/directory/where/experimental/data/is/located/'
topdir='./Data/'

# experimental friction coefficient
mu=0.3

#Change this if multiple experiments were used and use this to locate the correct data per experiment
experiment_nums=['1', '2']





# Loop over experiment
for experiment in experiment_nums:
        # Extract some information from the data set:
        try:
            data = np.loadtxt(topdir+experiment+'/Adjacency_list.txt', delimiter=',')
            #nsteps is the maximal frame number
            nsteps = np.max(data[:,0]).astype(int)
        except:
            print('No valid data found for experiment ' + str(experiment))
            #if no data is found for an experiment number go to the next experiment
            continue
        
        # Loop over strain steps for a given experiment
        # Start at 1 since steps start at 1. Ends at nsteps.
        for u in range(1, nsteps+1):
                #Creating configuration
                ThisConf = CF.Configuration(topdir+experiment,'experiment',mu, u)
                #Reading in the data
                ThisConf.ReadExpdata(verbose=False)
                #Adding boundary contacts
                #ThisConf.AddBoundaryContacts()
               
                #Setting up and playing the pebble game
                ThisPebble = PB.Pebbles(ThisConf,3,3,'nothing',False)
                #Play pebble game
                ThisPebble.play_game()
                # compute rigid clusters
                cidx, clusterall, clusterallBonds, clusteridx, BigCluster=ThisPebble.rigid_cluster()
                

                ########### Setting up the dynamical matrix and getting eigenmodes
                # This itself does very little, just creates an empty Hessian class
                # __init__(self,conf0):
                ThisHessian = HS.Hessian(ThisConf)
                
                ########## Have a look at some analysis functions of the rigid clusters
                #def __init__(self,conf0,pebbles0,hessian0,verbose=False):
                ThisAnalysis=AN.Analysis(ThisConf,ThisPebble,ThisHessian,0.01,True)
                # stress statistics
                #zav,nm,pres,fxbal,fybal,torbal,mobin,mohist,sxx,syy,sxy,syx=ThisAnalysis.getStressStat()
                # cluster statistics
                frac,fracmax,lenx,leny=ThisAnalysis.clusterStatistics()
                #def plotStresses(self,plotCir,plotVel,plotCon,plotF,plotStress,**kwargs):
                #fig1 = ThisAnalysis.plotStresses(True,False,False,True,False)
                #def plotPebbles(self,plotCir,plotPeb,plotPebCon,plotClus,plotOver,**kwargs):
                #ThisAnalysis.plotPebbles(True,True,True,False,False)
                
                #Plot pebbles has the following arguments: plotCir,plotPeb,plotPebCon,plotClus,plotOver
                fig2 = ThisAnalysis.plotPebbles(True,True,False,True,False)
                #fig2 = ThisAnalysis.plotPebbles(True,True,True,True,False)
                
                
                ######### continuing with the Hessian now 
                # constructing the matrix
                #  makeHessian(self,frictional,recomputeFnor,stabilise,verbose=False):
                #ThisHessian.makeHessian(True,False,0,False)
                # diagonalising the matrix
                # def getModes(self,debug=False):
                #ThisHessian.getModes(False)
                
                ##### a couple of checks on the modes (optional)
                #plotModes(self,usepts):
                #usepts=[3*ThisConf.N-5,3*ThisConf.N-4,3*ThisConf.N-3,3*ThisConf.N-2,3*ThisConf.N-1]
                #ThisHessian.plotModes(usepts)
                #plotZeroModes(self,thresh=2e-8,simple=True):
                #ThisHessian.plotZeroModes()
                
                ############ Now look for the cross-correlations
                # what is rigid according to modes with a given threshold:
                #fig3 = ThisAnalysis.ModeClusters('translations',2e-4)
                # how this cross-correlates to rigidy according to pebbles:
                #P_eig_if_pebble,P_pebble_if_eig,fig5 = ThisAnalysis.RigidModesCorrelate(2e-4)
                # These are the conditional probabilities of being rigid by mode while being rigid by pebble and the reverse
                #print (P_eig_if_pebble,P_pebble_if_eig)
                # if there is a next data set
                # Needs revision in any case, don't use for now
                #if ThisConf.Nnext>0:
                #    P_disp_if_pebble,P_pebble_if_disp, fig6 = ThisAnalysis.RigidDisplacementsCorrelate(2e-4)
                #    # Conditional probabilities of being rigid by displacement while being rigid by pebble
                #    print (P_disp_if_pebble,P_pebble_if_disp)
                #    # D2_min, needs assessment
                #    fig7 = ThisAnalysis.DisplacementCorrelateD2min(True)    
                plt.show()                              
