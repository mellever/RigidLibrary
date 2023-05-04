# Silke Henkes, 29.08.18:
# Created Configuration class to read both experimental and simulation data
# Detangled code parts, contains:
# - Read-in for all circumstances
# - Analysis functions which depend on particle positions only, called by Analysis and others

# Silke Henkes 11.07.2013: Accelerated version of pebble game and rigid cluster code
#!/usr/bin/python

import sys
import os
import glob
import numpy as np
import pandas
import math 
import matplotlib.pyplot as plt
import random



class Configuration:

        #=================== Create a new configruation ===============
        # Current choices: either simulation or experimental through datatype. Parameters are either read (simulation) or passed to the configuration via optional arguments (experiment)
        def __init__(self, folder, datatype, mu0, strainstep):
            self.folder = folder
            self.datatype = datatype
            # These get set to true if add boudary functions are called
            self.addBoundarySquare = False
            self.addBoundaryAnnulus = False
            if (datatype == 'simulation'):
                print("Reading simulation data!")
                self.getParameters(folder)
                # Simulation data has periodic boundary conditions and also angles as output
                self.periodic = True
                self.hasAngles = True
                # basis for mass computations
                self.density = 1.0
                # Prefactor of the stiffness coefficients
                self.stiffness = 1.0
                # radius conversion factor
                self.rconversion = 1.0
                self.height = 1.0
                self.width = 1.0
                self.test = False
            elif (datatype == 'simulation_test'):
                print("Reading simulation data!")
                # Simulation data has periodic boundary conditions and also angles as output
                self.periodic = True
                self.hasAngles = True
                # basis for mass computations
                self.density = 1.0
                # Prefactor of the stiffness coefficients
                self.stiffness = 1.0
                # radius conversion factor
                self.rconversion = 1.0
                self.height = 1.0
                self.width = 1.0
                self.test = True
            elif (datatype == 'experiment_square'):
                self.step = strainstep
                self.mu = mu0
                # Experimental data does not have periodic boundary conditions and there is no angle data either
                self.periodic = False
                self.hasAngles = False
                # density of the material in kg/m^3
                self.density = 1060.0
                # Prefactor of the stiffness coefficient, k_n = \frac{\pi}{6} E h = 6490 kg /s^2
                self.stiffness = 6490.0
                # radius conversion factor from pixel to m.
                self.rconversion = 2.7e-4
                # height of the disks in m
                self.height = 3.1e-3
                # boundary width
                self.width = 20e-3
                # Prefixes of the data name
                self.prefix1='DSC'
                self.prefix2='_solved_Tracke_'
            elif (datatype == 'experiment_annulus'):
                #Radius of the inner ring in experiment
                self.R1 = 1548
                #Radius of the outer ring in experiment
                self.R2 = 2995
                #Coordinates of midpoint
                self.mid_x = 3134
                self.mid_y = 3028
                #Size of boundary texture
                self.texture = 30
                #Strainstep of the experiment
                self.step = strainstep
                #Set friction coefficient
                self.mu = mu0
                # Experimental data does not have periodic boundary conditions and there is no angle data either
                self.periodic = False
                self.hasAngles = False
                # density of the material in kg/m^3
                self.density = 1060.0
                # Prefactor of the stiffness coefficient, k_n = \frac{\pi}{6} E h = 6490 kg /s^2
                self.stiffness = 6490.0
                # radius conversion factor from pixel to m.
                self.rconversion = 2.7e-4
                # height of the disks in m
                self.height = 3.1e-3
                # boundary width
                self.width = 20e-3
                self.experiment = True
            elif (datatype == 'test'):
                # Use this for analysing artificial packings
                #Strainstep of the experiment
                self.step = strainstep
                #Set friction coefficient
                self.mu = mu0
                # Experimental data does not have periodic boundary conditions and there is no angle data either
                self.periodic = False
                self.hasAngles = False
                # density of the material in kg/m^3
                self.density = 1060.0
                # Prefactor of the stiffness coefficient, k_n = \frac{\pi}{6} E h = 6490 kg /s^2
                self.stiffness = 6490.0
                # radius conversion factor from pixel to m.
                self.rconversion = 2.7e-4
                # height of the disks in m
                self.height = 3.1e-3
                # boundary width
                self.width = 20e-3
                self.experiment = True
            else:
                print("Error: Unknown type of data! Doing nothing.")

        #======= Simulation parameter read-in =====================
        # use Parameter read-in as function called by the constructor
        def getParameters(self, folder):
                parname = folder + 'Parameters_final.dat'
                if os.path.exists(parname):    column = 1
                else:
                        parname = folder + 'Parameters_LE.par'
                        if os.path.exists(parname):
                            column = 2
                        else:
                            print('Error: Couldn\'t find parameter file.')
                parfile = open(parname, 'r')
                for line in parfile:
                        s = line.split()
                        if s == []: s = ['']
                        if s[0] == 'np': self.N = eval(
                            s[column])        # number of particles
                        # spring constant
                        elif s[0] == 'k_rep': self.krep = eval(s[column])
                        # system size
                        elif s[0] == 'boxVec1x': self.L = eval(s[column])
                        # friction coefficient
                        elif s[0] == 'mu': self.mu = eval(s[column])
                        # viscous damping coefficient
                        elif s[0] == 'xivisc': self.xi = eval(s[column])
                        # packing fraction
                        elif s[0] == 'phi': self.phi = eval(s[column])
                        # packing fraction
                        elif s[0] == 'strain_rate': self.gammadot = eval(s[column])
                parfile.close()
                self.Lx = self.L
                self.Ly = self.L
                self.rad = np.loadtxt(folder + 'posinit.dat', usecols=(3,))

        #######========== Simulation data read-in ==================
        # snap is the label, and distSnapshot tells me how many time units they are apart (in actual time, not steps)
        def readSimdata(self, snap, verbose0, distSnapshot0=1000):
                self.verbose = verbose0
                self.distSnapshot = distSnapshot0
                if not self.test:
                    self.strain = self.gammadot*(1.0*snap)*self.distSnapshot

                if (self.verbose):
                        print('L=' + str(self.L))
                        print('Np=' + str(self.N))
                        print('mu=' + str(self.mu))
                        print('xi=' + str(self.xi))
                        print('phi=' + str(self.phi))
                        print('strain=' + str(self.strain))
                filename = self.folder + '/PosVel' + str(snap) + '.dat'
                if (self.verbose):
                        print(filename)
                isPosdata = False
                # Only if: 1) The file exists 2) is not empty
                try:
                        data = np.loadtxt(filename)
                        if (len(data.shape) > 1):
                                isPosdata = True
                except:
                        isPosdata = False
                

                if (not isPosdata):
                        if self.verbose:
                            print('error: no position data at snapshot' + str(snap))
                        self.x = 0.0
                        self.y = 0.0
                        self.dx = 0.0
                        self.dy = 0.0
                        self.alpha = 0.0
                        self.N = 0
                else:   
                        self.x = data[:,0]
                        self.y = data[:,1]
                        self.alpha = data[:,2]
                        self.dx = data[:,3]
                        self.dy = data[:,4]
                        self.dalpha = data[:,5]
                        del data

                if self.test: 
                    self.Lx = np.amax(self.x)-np.amin(self.x)
                    self.Ly = np.amax(self.y)-np.amin(self.y)
                    self.L = np.amax([self.Lx, self.Ly])
                    self.strain = 0
                
                filename = self.folder + '/Contact' + str(snap) + '.dat'
                if (self.verbose):
                        print(filename)
                isCondata = False
                # Only if: 1) The file exists 2) is not empty 3) is longer than 1 contact (garbage avoiding device .. hope for the best. Some of these got cut off mid-writing)
                try:
                        data = np.loadtxt(filename)
                        if (len(data.shape) > 1):
                                isCondata = True
                except:
                        isCondata = False
                if ((not isCondata) or (not isPosdata)):
                        if self.verbose:
                            print('error: no position data at snapshot' + str(snap))
                        self.noConf = True
                        self.I = -1
                        self.J = -1
                        self.fullmobi = 0
                        self.Ifull = -1
                        self.Jfull = -1
                else:
                        self.noConf = False
                        self.I = data[:,0].astype(int)
                        self.J = data[:,1].astype(int)
                        self.nx = data[:,2]
                        self.ny = data[:,3]
                        self.fullmobi = data[:,4].astype(int)
                        self.fnor = data[:,5]
                        self.ftan = data[:,6]+data[:,7]
                        del data
                        if not self.test:
                            self.x -= self.L*np.round(self.x/self.L)
                            self.y -= self.L*np.round(self.y/self.L)
                        self.ncon = len(self.I)
        #========== Experimental data read-in for annulus ==================
        def ReadExpdataAnnulusPandas(self, verbose):
                prefix = self.folder +'/particle_positions.txt'
                try:
                    #read in data of following format:
                    #framenumber, particle id, x, y, radius, boundary
                    pos_data = pandas.read_csv(prefix, names=['n','id', 'x', 'y', 'r', 'b'])
                    pos_data = pos_data[pos_data['n'] == self.step]
                    self.id = pos_data.loc[:,'id'].to_numpy()
                    self.x = pos_data.loc[:,'x'].to_numpy()
                    self.y = pos_data.loc[:,'y'].to_numpy()
                    self.rad = pos_data.loc[:,'r'].to_numpy()
                    self.boundary = pos_data.loc[:,'b'].to_numpy()
                    self.N = len(self.rad)
                    self.Lx = np.amax(self.x)-np.amin(self.x)
                    self.Ly = np.amax(self.y)-np.amin(self.y)
                except:
                    self.x = 0
                    self.y = 0
                    self.rad = 1  # purely so that we can divide by it ...
                    self.N = 0
                    print("Error: there is no position data here")
                    return 1
                
                prefix = self.folder +'/Adjacency_list.txt'
                try:
                    #read in data of following format:
                    #framenumber, id1, id2, ft, fn
                    con_data = pandas.read_csv(prefix, names=['n', 'id1', 'id2', 'ft', 'fn'])
                    #length of previous framenumer
                    if self.step == 1:
                        len_prev = 0
                    else:
                        len_prev = len(con_data[con_data['n'] <= self.step-1]['id1'])
                    #data for current framenumber
                    con_data= con_data[con_data['n'] == self.step]
                    
                except:
                    print("Error: there is no contact data here")
                    return 1
                            
                #Temporary
                con_frame = con_data
                
                #Create lists
                self.I=[]
                self.J=[]
                fn0=[]
                ft0=[]
                fm0=[]                   

                #Filter entries and determine if slipping or not
                for k in range(len_prev, len_prev + len(con_frame['id1'])):    
                    try:
                        #Add particle id's, identified by indices, to list.
                        i = con_frame['id1'][k]
                        j = con_frame['id2'][k]
                        argi = np.argwhere(self.id == i).flatten()[0]
                        argj = np.argwhere(self.id == j).flatten()[0]
                        self.I.append(argi)
                        self.J.append(argj)
                        
                        # Search other force and take the mean, might not be relevant.
                        norm2 = con_frame[ (con_frame['id1']==con_frame['id2'][k]) & (con_frame['id2']==con_frame['id1'][k])]['fn'].iloc[0]
                        tan2 = con_frame[ (con_frame['id1']==con_frame['id2'][k]) & (con_frame['id2']==con_frame['id1'][k])]['ft'].iloc[0]
                        norm = (con_frame['fn'][k] + norm2)/2
                        tan = (con_frame['ft'][k] + tan2)/2
                        
                        fn0.append(norm)
                        ft0.append(tan)
                        if (abs(tan)/norm>self.mu):
                            fm0.append(1)
                        else:
                            fm0.append(0)
                        #Now delete corresponding entry to remove duplicates. 
                        #Check: len(self.I) = (#True in dups_bool)/2 + (#False in dups_bool)
                        con_frame.drop(con_frame[(con_frame['id1']==con_frame['id2'][k]) & (con_frame['id2']==con_frame['id1'][k])].index, inplace=True)
                    except:
                        if verbose:
                            print('dropped duplicate')
                                  
                #Final arrays
                self.ncon=len(fm0)
                self.I = np.array(self.I)
                self.J = np.array(self.J)
                self.fnor=np.array(fn0)
                self.ftan=np.array(ft0)
                self.fullmobi=np.array(fm0)

                self.nx=np.zeros(self.ncon)
                self.ny=np.zeros(self.ncon)
                for k in range(self.ncon):
                        x1=self.x[self.I[k]]
                        y1=self.y[self.I[k]]
                        x2=self.x[self.J[k]]
                        y2=self.y[self.J[k]]
                        rij=np.sqrt((x1-x2)**2+(y1-y2)**2)
                        self.nx[k]=(x2-x1)/rij
                        self.ny[k]=(y2-y1)/rij
                
                print("Config frame #" +str(self.step)+ " created")      
                return 0
        
        def ReadExpdataAnnulusNumpy(self, verbose):    
                prefix = self.folder +'/particle_positions.txt'
                try:
                    coords=np.loadtxt(prefix, delimiter=',')
                    coords=coords[coords[:,0] == self.step]
                    self.id=coords[:,1]
                    self.x=coords[:,2]
                    self.y=coords[:,3]
                    self.rad=coords[:,4]
                    self.boundary=coords[:,5]
                    self.N=len(self.rad)
                    self.Lx=np.amax(self.x)-np.amin(self.x)
                    self.Ly=np.amax(self.y)-np.amin(self.y)
                    del coords
                except:
                    isPosdata=False
                    self.x=0
                    self.y=0
                    self.rad=1 # purely so that we can divide by it ...
                    self.N=0
                    print("Error: there is no position data here")
                #Load in contact data
                prefix = self.folder +'/Adjacency_list.txt'
                try:
                    condata=np.loadtxt(prefix, delimiter=',')
                    condata=condata[condata[:,0] == self.step]
                except:
                    isCondata=False
                    print ("Error: there is no contact data here")
                    return 1
                
                #Create empty lists
                #Lists with particle ids
                self.I=[]
                self.J=[]
                
                #Lists with particle indices, ommit for now.
                #self.argI=[]
                #self.argJ=[]
                
                #Lists with forces
                fn0=[]
                ft0=[]
                
                #Frictional or sliding list
                fm0=[]
                
                #Drop duplicates
                for k in range(len(condata[:,0])):
                    if condata[:,1][k] > condata[:,2][k]:
                        #Add contact id's
                        i = condata[:,1][k].astype(int)
                        j = condata[:,2][k].astype(int)
                        argi = np.argwhere(self.id == i).flatten()[0]
                        argj = np.argwhere(self.id == j).flatten()[0]
                        
                        """
                        #self.I.append(i)
                        #self.J.append(j)
                        self.argI.append(argi)
                        self.argJ.append(argj)
                        """
                        
                        #For now we identify particle id's with indices in the list. Not ideal for debugging, but it works.
                        #Can be fixed if necessary, see above code for start. 
                        self.I.append(argi)
                        self.J.append(argj)
                        
                        #Extract forces
                        fn = condata[:,4][k]
                        ft = condata[:,3][k]
                        fn0.append(fn)
                        ft0.append(ft)
                        
                        #Determine frictional or sliding
                        if (abs(ft)/fn>self.mu):
                            fm0.append(1)
                        else:
                            fm0.append(0)
                
                #Final arrays
                self.ncon=len(fm0)
                self.I = np.array(self.I)
                self.J = np.array(self.J)
                self.fnor=np.array(fn0)
                self.ftan=np.array(ft0)
                self.fullmobi=np.array(fm0)
                
                self.nx=np.zeros(self.ncon)
                self.ny=np.zeros(self.ncon)
                for k in range(self.ncon):
                        x1=self.x[self.I[k]]
                        y1=self.y[self.I[k]]
                        x2=self.x[self.J[k]]
                        y2=self.y[self.J[k]]
                        rij=np.sqrt((x1-x2)**2+(y1-y2)**2)
                        self.nx[k]=(x2-x1)/rij
                        self.ny[k]=(y2-y1)/rij

                print("Config frame #" +str(self.step)+ " created")      
                return 0
        
        
        
        #========== Maxwell Cremona tiling -> port this to its own class ==================
        def Tiling(self, graph, tiling):
            #Plotting the graph
            if graph:
                for i in self.I:
                    for k in self.J[self.I==i]:
                        x0, x1, y0, y1 = getConPos2(i, k)
                        plt.plot([x0, x1], [y0, y1])
                    break
            plt.show()
            
            #Plotting the Maxwell Cremona tiling
            if tiling:
                color = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)]) for i in range(np.max(self.I)+1)]
                for i in self.I:
                    #Scaling parameter, remove later
                    scale = 1000
                    #Origin coordinates
                    xor = self.x[i]
                    yor = self.y[i]
                    #Determine contacts
                    x = np.argwhere(self.J==i).flatten()
                    y = np.argwhere(self.I==i).flatten()
                    print(y)
                    z = np.concatenate([x,y])
                    xhat = self.J[y]
                    yhat = self.I[x]
                    zhat = np.concatenate([xhat,yhat])
                    #Plot
                    for k in range(len(z)):
                        theta = np.arctan2(self.y[i]-self.y[zhat[k]], self.x[i]-self.x[zhat[k]])
                        fx = scale*self.fnor[z[k]]*np.cos(theta)
                        fy = scale*self.fnor[z[k]]*np.sin(theta)
                        print(i,k,theta, self.fnor[z[k]])
                        plt.arrow(xor, yor, fx, fy, width=.08, facecolor=color[i])
                        #plt.plot([xor, xor+fx], [yor, yor+fy])
                        xor += fx
                        yor += fy
                plt.show()                



         #========== Experimental data read-in for square ==================
        def ReadExpdataSquare(self,numlabel):    
                prefix = self.prefix1 + numlabel + self.prefix2
                print ("Starting analysis of experimental step " + prefix)
                # Let's get the positions first
                isPosdata=True
                print (self.folder+prefix + 'ParticleData.dlm')
                try:
                    coords=np.loadtxt(self.folder+prefix + 'ParticleData.dlm',delimiter=',')
                    self.x=coords[:,1]
                    self.y=coords[:,2]
                    self.rad=coords[:,3]
                    self.N=len(self.rad)
                    self.Lx=np.amax(self.x)-np.amin(self.x)
                    self.Ly=np.amax(self.y)-np.amin(self.y)
                    del coords
                except:
                    isPosdata=False
                    self.x=0
                    self.y=0
                    self.rad=1 # purely so that we can divide by it ...
                    self.N=0
                    print("Error: there is no position data here")
                
                # The contact data lives in the adjacency matrices
                # As non-sparse N by N matrices (!)
                isCondata=True
                try:
                    contacts=np.loadtxt(self.folder+prefix+'BinaryAdjacencyMatrix.dlm',delimiter=',')
                    fnor0=np.loadtxt(self.folder+prefix+'NormWeightedAdjacencyMatrix.dlm',delimiter=',')
                    ftan0=np.loadtxt(self.folder+prefix+'TanWeightedAdjacencyMatrix.dlm',delimiter=',')
                except:
                    isCondata=False
                    print ("Error: there is no contact data here")
                    return 1
                    
                if isCondata:
                    self.I=[]
                    self.J=[]
                    fn0=[]
                    ft0=[]
                    fm0=[]
                    for k in range(self.N):
                        thisline=contacts[k,:].astype(int)
                        iscontact=np.nonzero(thisline)
                        for idx in iscontact[0]:
                            isnew=True
                            # Now check if this contact has been seen before
                            # Usually yes, except for some quite weak contacts, which show up in only one direction (and boundary?!)
                            if idx<k:
                                isidx=np.nonzero(self.I==idx)
                                if len(isidx[0])>0:
                                    for idx2 in isidx[0]:
                                        # if k itself shows up in this list it's not a new contact
                                        if self.J[idx2]==k:
                                            isnew=False
                                            #print self.J[idx2]
                            # if it's new, add it to the list
                            if isnew:
                                self.I.append(k)
                                self.J.append(idx)
                                fn0.append(fnor0[k,idx])
                                ft0.append(ftan0[k,idx])
                                if (abs(ftan0[k,idx])/fnor0[k,idx]>self.mu):
                                    fm0.append(1)
                                else:
                                    fm0.append(0)
                            # otherwise merge it with its twin
                            else: 
                                # take the average of the forces
                                fn0[idx2]=0.5*(fn0[idx2]+fnor0[k,idx])
                                ft0[idx2]=0.5*(ft0[idx2]+ftan0[k,idx])
                                # and revisit the mobilisation issue
                                if (abs(ft0[idx2])/fn0[idx2]>self.mu):
                                    fm0[idx2]=1
                                else:
                                    fm0[idx2]=0
                    self.ncon=len(fm0)
                    self.fnor=np.array(fn0)
                    self.ftan=np.array(ft0)
                    self.fullmobi=np.array(fm0)
                    # clean up, these are large matrices
                    del contacts
                    del fnor0
                    del ftan0
                else:
                    self.I=-1
                    self.J=-1
                    self.fullmobi=0
                    self.Ifull=-1
                    self.Jfull=-1
                    self.ncon=1
                print ("Read in a configuration of " + str(self.N) + " particles with " + str(self.ncon) + " unique contacts.")
                
                # create the local normal coordinate system. Note that experiment never has periodic boundary conditions
                self.nx=np.zeros(self.ncon)
                self.ny=np.zeros(self.ncon)
                for k in range(self.ncon):
                        x1=self.x[self.I[k]]
                        y1=self.y[self.I[k]]
                        x2=self.x[self.J[k]]
                        y2=self.y[self.J[k]]
                        rij=np.sqrt((x1-x2)**2+(y1-y2)**2)
                        self.nx[k]=(x2-x1)/rij
                        self.ny[k]=(y2-y1)/rij
                        
                return 0

        # same kind of procedure, but get only the next positions and don't redefine things
        # This is for comparison of the displacements compared to the the next dataset to pebbles and modes
        def ReadExpdataNext(self,numlabel,scale=False):
                prefix = self.prefix1 + numlabel + self.prefix2
                print ("Reading experimental step " + prefix + " as next data set.")
                # Let's get the positions first
                isPosdata=True
                try:
                    coords=np.loadtxt(self.folder+prefix + 'ParticleData.dlm',delimiter=',')
                    self.xnext=coords[:,1]
                    self.ynext=coords[:,2]
                    self.radnext=coords[:,3]
                    self.Nnext=len(self.rad)
                    self.Lxnext=np.amax(self.x)-np.amin(self.x)
                    self.Lynext=np.amax(self.y)-np.amin(self.y)
                    del coords
                except:
                    isPosdata=False
                    self.xnext=0
                    self.ynext=0
                    self.radnext=1 # purely so that we can divide by it ...
                    self.Nnext=0
                # Slight overkill, but then this enforces unified syntax with simulation, making our job below easier
                if len(self.xnext)==len(self.x):
                    self.dx=self.xnext-self.x
                    self.dy=self.ynext-self.y
         
        #### ======================== Boundary integration =======================================================
        def AddBoundaryContactsSquare(self,threshold=20,Brad=20.0):
            self.addBoundarySquare=True
            # Threshold to check if a particle is close enough to walls.
            upidx=np.argmax(self.y)
            downidx=np.argmin(self.y)
            leftidx=np.argmin(self.x)
            rightidx=np.argmax(self.x)
            
            # Boundary posiitons:
            # coordinates of virtual boundary particles: in the middle, one Brad off from the edge of the outermost particle
            up=self.y[upidx]
            yup = up+self.rad[upidx]
            down=self.y[downidx]
            ydown = down-self.rad[downidx]
            left=self.x[leftidx]
            xleft=left-self.rad[leftidx]
            right=self.x[rightidx]
            xright=right+self.rad[rightidx]
            
            # coordinates of virtual boundary particles: in the middle, one Brad off from the edge of the outermost particle
            Boundaries=np.zeros((4,3)) # Four boundary particles with their x,y and rad
            Boundaries[0,:]=[(left+right)*0.5,yup+Brad,Brad]
            Boundaries[1,:]=[(left+right)*0.5,ydown-Brad,Brad]
            Boundaries[2,:]=[xleft-Brad,(up+down)*0.5,Brad]
            Boundaries[3,:]=[xright+Brad,(up+down)*0.5,Brad]
            
            # Find the particles in contact with the boundary, and label correctly
            self.bindices=[self.N,self.N+1,self.N+2,self.N+3]
            padd=[]
            labels=[]
            pup =  np.nonzero(np.abs(self.y+self.rad-yup)<threshold)[0]
            padd.extend(pup)
            labels.extend([0 for k in range(len(pup))])
            pdown =  np.nonzero(np.abs(self.y-self.rad-ydown)<threshold)[0]
            padd.extend(pdown)
            labels.extend([1 for k in range(len(pdown))])
            pleft = np.nonzero(np.abs(self.x-self.rad-xleft)<threshold)[0]
            padd.extend(pleft)
            labels.extend([2 for k in range(len(pleft))])
            pright = np.nonzero(np.abs(self.x+self.rad-xright)<threshold)[0]
            padd.extend(pright)
            labels.extend([3 for k in range(len(pright))])
            
            fullmobi_add=[]
            fnor_add=[]
            ftan_add=[]
            nx_add=[]
            ny_add=[]
            for k in range(len(padd)):
                # does this guy have neighbours?
                neii=np.nonzero(self.I[:self.ncon]==padd[k])[0]
                neij=np.nonzero(self.J[:self.ncon]==padd[k])[0]
                # if yes add the boundary contacts
                if (len(neii)>0 or len(neij)>0):
                    self.I.append(self.bindices[labels[k]])
                    self.J.append(padd[k])
                    if (labels[k])==0:
                        nx0=0
                        ny0=-1
                    elif (labels[k]==1):
                        nx0=0
                        ny0=1
                    elif (labels[k]==2):
                        nx0=1
                        ny0=0
                    else:
                        nx0=-1
                        ny0=0
                    # compute force on this contact by force balance
                    # two minus signs on second part cancel out
                    ftotx=np.sum(self.fnor[neii]*self.nx[neii]-self.ftan[neii]*self.ny[neii])-np.sum(self.fnor[neij]*self.nx[neij]-self.ftan[neij]*self.ny[neij])
                    ftoty=np.sum(self.fnor[neii]*self.ny[neii]+self.ftan[neii]*self.nx[neii])-np.sum(self.fnor[neij]*self.ny[neij]+self.ftan[neij]*self.nx[neij])
                    # (fx*nx+fy*ny)
                    fnor0=ftotx*nx0+ftoty*ny0
                    # (fx*(-ny)+fy*nx)
                    ftan0=ftotx*(-ny0)+ftoty*nx0
                    #print (ftan0)
                    if (abs(ftan0)/fnor0>self.mu):
                        fullmobi_add.append(1)
                    else:
                        fullmobi_add.append(0)
                    fnor_add.append(fnor0)
                    ftan_add.append(ftan0)
                    nx_add.append(nx0)
                    ny_add.append(ny0)
            # Finally stick it at the end of the existing data
            self.x=np.concatenate((self.x,Boundaries[:,0]))
            self.y=np.concatenate((self.y,Boundaries[:,1]))
            self.rad=np.concatenate((self.rad,Boundaries[:,2]))
            self.fnor=np.concatenate((self.fnor,np.array(fnor_add)))
            self.ftan=np.concatenate((self.ftan,np.array(ftan_add)))
            self.fullmobi=np.concatenate((self.fullmobi,np.array(fullmobi_add)))
            self.nx=np.concatenate((self.nx,np.array(nx_add)))
            self.ny=np.concatenate((self.ny,np.array(ny_add)))
            self.ncon=len(self.I)
            self.N+=4
            print ("Added boundaries!")
        
        #### ======================== Boundary integration =======================================================
        def AddBoundaryContactsAnnulus2(self,threshold=30):
            # For getting positions
            self.addBoundaryAnnulus=True

            #Set radius of boundary particles
            #This needs some tweaking, since we are working with innies instead of outies. 
            Brad = self.texture
            
            # Boundary positions:
            # Coordinates of virtual boundary particles: at mid_x and y one Brad off from radi
            b1 = self.mid_y + self.R1 - Brad
            b2 = self.mid_y + self.R2 + Brad
                        
            # Coordinates of virtual boundary particles: in the middle, one Brad off from the edge of the outermost particle
            Boundaries=np.zeros((2,3)) # Two boundary particles with their x,y and rad
            Boundaries[0,:]=[self.mid_x,b1,Brad]
            Boundaries[1,:]=[self.mid_x,b2,Brad]
            
            # Find the particles in contact with the boundary, and label correctly
            self.bindices=[self.N,self.N+1]
            padd=[]
            labels=[]
            
            #Distance of points from the center of the annulus
            r = np.sqrt(np.square(self.x- self.mid_x)+np.square(self.y-self.mid_y))
            p1 =  np.nonzero(np.abs(r - self.R1 - self.rad)<threshold)[0]
            padd.extend(p1)
            labels.extend([0 for k in range(len(p1))])
            p2 =  np.nonzero(np.abs(self.R2 - r - self.rad)<threshold)[0]
            padd.extend(p2)
            labels.extend([1 for k in range(len(p2))])
            
            #Some arrays for the information of the newly introduced contacts
            fullmobi_add=[]
            fnor_add=[]
            ftan_add=[]
            nx_add=[]
            ny_add=[]
            for k in range(len(padd)):
                self.I = np.append(self.I, self.bindices[labels[k]])
                self.J = np.append(self.J, padd[k])
                neii=np.nonzero(self.I[:self.ncon]==padd[k])[0]
                neij=np.nonzero(self.J[:self.ncon]==padd[k])[0]

                #Always add double bound
                fullmobi_add.append(0)
                
                #Compute forces if neighbours are present
                if (len(neii)>0 or len(neij)>0):
                    # Flip direction of force, is this correct?
                    nx0 = ny0 = -1

                    # Compute force on this contact by force balance
                    # Two minus signs on second part cancel out
                    ftotx=np.sum(self.fnor[neii]*self.nx[neii]-self.ftan[neii]*self.ny[neii])-np.sum(self.fnor[neij]*self.nx[neij]-self.ftan[neij]*self.ny[neij])
                    ftoty=np.sum(self.fnor[neii]*self.ny[neii]+self.ftan[neii]*self.nx[neii])-np.sum(self.fnor[neij]*self.ny[neij]+self.ftan[neij]*self.nx[neij])
                    fnor0=ftotx*nx0+ftoty*ny0
                    ftan0=ftotx*(-ny0)+ftoty*nx0

                else: fnor0 = ftan0 = nx0 = ny0 = 1
                
                #Add data to list
                fnor_add.append(fnor0)
                ftan_add.append(ftan0)
                nx_add.append(nx0)
                ny_add.append(ny0)
            
            #Add slipping contact between the two boundary particles
            self.I = np.append(self.I, self.bindices[0])
            self.J = np.append(self.J, self.bindices[1])
            fullmobi_add.append(1)

            #Forces between boundary particles are zero?
            fnor0 = ftan0 = nx0 = ny0 = 0

            #Add to lists
            fnor_add.append(fnor0)
            ftan_add.append(ftan0)
            nx_add.append(nx0)
            ny_add.append(ny0)

            # Finally stick it at the end of the existing data
            self.x=np.concatenate((self.x,Boundaries[:,0]))
            self.y=np.concatenate((self.y,Boundaries[:,1]))
            self.rad=np.concatenate((self.rad,Boundaries[:,2]))
            self.fnor=np.concatenate((self.fnor,np.array(fnor_add)))
            self.ftan=np.concatenate((self.ftan,np.array(ftan_add)))
            self.fullmobi=np.concatenate((self.fullmobi,np.array(fullmobi_add)))
            self.nx=np.concatenate((self.nx,np.array(nx_add)))
            self.ny=np.concatenate((self.ny,np.array(ny_add)))
            self.ncon=len(self.I)
            self.N+=2 #Since we have two boundary particles
            print ("Added boundaries!")

        #### ======================== Boundary integration =======================================================
        def AddBoundaryContactsAnnulus(self):
            # For getting positions
            self.addBoundaryAnnulus=True
            
            #Set radius of boundary particles
            #This needs some tweaking, since we are working with innies instead of outies. 
            Brad = self.texture
            Brad = 300
            
            # Boundary positions:
            # Coordinates of virtual boundary particles: at mid_x and y one Brad off from radi
            b1 = self.mid_y + self.R1 - Brad #Inner
            b2 = self.mid_y + self.R2 + Brad #Outer
                        
            # Coordinates of virtual boundary particles: in the middle, one Brad off from the edge of the outermost particle
            Boundaries=np.zeros((2,3)) # Two boundary particles with their x,y and rad
            Boundaries[0,:]=[self.mid_x,b1,Brad] #Inner
            Boundaries[1,:]=[self.mid_x,b2,Brad] #Outer
            
            #Generate two new id's for boundary particles
            self.bindices=[self.N, self.N +1]
            padd=[]
            labels=[]
            
            #If in contact with inner boundart
            p1 =  np.nonzero(self.boundary == -1)[0]
            padd.extend(p1)
            labels.extend([0 for k in range(len(p1))])
            #If in contact with outer boundary
            p2 =  np.nonzero(self.boundary == 1)[0]
            padd.extend(p2)
            labels.extend([1 for k in range(len(p2))])
            
            #Some arrays for the information of the newly introduced contacts
            fullmobi_add=[]
            fnor_add=[]
            ftan_add=[]
            nx_add=[]
            ny_add=[]
            for k in range(len(padd)):
                self.I = np.append(self.I, self.bindices[labels[k]])
                self.J = np.append(self.J, padd[k])
                neii=np.nonzero(self.I[:self.ncon]==padd[k])[0]
                neij=np.nonzero(self.J[:self.ncon]==padd[k])[0]

                #Always add double bound
                fullmobi_add.append(0)
                
                #Compute forces if neighbours are present
                if (len(neii)>0 or len(neij)>0):
                    # Flip direction of force, is this correct?
                    nx0 = ny0 = -1

                    # Compute force on this contact by force balance
                    # Two minus signs on second part cancel out
                    ftotx=np.sum(self.fnor[neii]*self.nx[neii]-self.ftan[neii]*self.ny[neii])-np.sum(self.fnor[neij]*self.nx[neij]-self.ftan[neij]*self.ny[neij])
                    ftoty=np.sum(self.fnor[neii]*self.ny[neii]+self.ftan[neii]*self.nx[neii])-np.sum(self.fnor[neij]*self.ny[neij]+self.ftan[neij]*self.nx[neij])
                    fnor0=ftotx*nx0+ftoty*ny0
                    ftan0=ftotx*(-ny0)+ftoty*nx0

                else: fnor0 = ftan0 = nx0 = ny0 = 1
                
                #Add data to list
                fnor_add.append(fnor0)
                ftan_add.append(ftan0)
                nx_add.append(nx0)
                ny_add.append(ny0)
            
            #Add slipping contact between the two boundary particles
            self.I = np.append(self.I, self.bindices[0])
            self.J = np.append(self.J, self.bindices[1])
            fullmobi_add.append(1)

            #Forces between boundary particles are zero?
            fnor0 = ftan0 = nx0 = ny0 = 0

            #Add to lists
            fnor_add.append(fnor0)
            ftan_add.append(ftan0)
            nx_add.append(nx0)
            ny_add.append(ny0)

            # Finally stick it at the end of the existing data
            self.x=np.concatenate((self.x,Boundaries[:,0]))
            self.y=np.concatenate((self.y,Boundaries[:,1]))
            self.id=np.concatenate((self.id, self.bindices))
            self.rad=np.concatenate((self.rad,Boundaries[:,2]))
            self.fnor=np.concatenate((self.fnor,np.array(fnor_add)))
            self.ftan=np.concatenate((self.ftan,np.array(ftan_add)))
            self.fullmobi=np.concatenate((self.fullmobi,np.array(fullmobi_add)))
            self.nx=np.concatenate((self.nx,np.array(nx_add)))
            self.ny=np.concatenate((self.ny,np.array(ny_add)))
            self.ncon=len(self.I)
            self.N+=2 #Since we have two boundary particles
            print ("Added boundaries!")

        def AddRandomContacts(self, percent, dmax):         
            #Create empty lists for storing data
            fullmobi_add=[]
            fnor_add=[]
            ftan_add=[]
            nx_add=[]
            ny_add=[]
            
            #Choose amount of particles random number in [1,100]
            nums = np.random.randint(1,100,len(self.x))

            #Loop over all particles
            for i in range(len(self.x)):
                #If number less than percent add bond, else continue
                if nums[i] >= percent: continue 
                
                #Loop over contacts in proximity
                for n in range(1,7):
                    skip = False
                    #Compute distance between particles
                    diffarr = np.sqrt(np.square(self.x-self.x[i])+np.square(self.y-self.y[i]))
                    #Take n'th smallest distance
                    diff = np.sort(diffarr)[n]
                    if diff >= dmax:
                        skip = True
                        break 
                    #Find corresponding index
                    arg = np.argwhere(diffarr == diff).flatten()[0] 
                    
                    #Only add new bonds
                    if arg in self.J[np.argwhere(self.I == i)].flatten():
                        skip = True
                        continue
                    if arg in self.I[np.argwhere(self.J == i)].flatten():
                        skip = True
                        continue
                
                #If distance is to large, do not add
                if skip:
                    continue
                
                #Append to contacts
                self.I = np.append(self.I, i)
                self.J = np.append(self.J, arg)
                neii=np.nonzero(self.I[:self.ncon]==arg)[0]
                neij=np.nonzero(self.J[:self.ncon]==arg)[0]

                #Always add double bound
                fullmobi_add.append(0)
                
                #Compute forces if neighbours are present
                if (len(neii)>0 or len(neij)>0):
                    # Flip direction of force, is this correct?
                    nx0 = ny0 = -1

                    # Compute force on this contact by force balance
                    # Two minus signs on second part cancel out
                    ftotx=np.sum(self.fnor[neii]*self.nx[neii]-self.ftan[neii]*self.ny[neii])-np.sum(self.fnor[neij]*self.nx[neij]-self.ftan[neij]*self.ny[neij])
                    ftoty=np.sum(self.fnor[neii]*self.ny[neii]+self.ftan[neii]*self.nx[neii])-np.sum(self.fnor[neij]*self.ny[neij]+self.ftan[neij]*self.nx[neij])
                    fnor0=ftotx*nx0+ftoty*ny0
                    ftan0=ftotx*(-ny0)+ftoty*nx0

                else: fnor0 = ftan0 = nx0 = ny0 = 1
                #Add data to list
                fnor_add.append(fnor0)
                ftan_add.append(ftan0)
                nx_add.append(nx0)
                ny_add.append(ny0)
            
            #Add everything to existing arrays
            self.fnor=np.concatenate((self.fnor,np.array(fnor_add)))
            self.ftan=np.concatenate((self.ftan,np.array(ftan_add)))
            self.fullmobi=np.concatenate((self.fullmobi,np.array(fullmobi_add)))
            self.nx=np.concatenate((self.nx,np.array(nx_add)))
            self.ny=np.concatenate((self.ny,np.array(ny_add)))
            self.ncon=len(self.I)
            print("Added Random Bonds")
                

        def AddSmartContacts(self,dmax,ang, threshold):         
            #Create empty lists for storing data
            fullmobi_add=[]
            fnor_add=[]
            ftan_add=[]
            nx_add=[]
            ny_add=[]
                    
            #Loop over all particles
            for i in self.I:
                #Determine indices of contacts and put together in one array
                argI = np.argwhere(self.I==i).flatten()
                argJ = np.argwhere(self.J==i).flatten()
                arg = np.concatenate([argI,argJ])
                #Get particle ids from indices and put together in one array
                I = self.I[argI]
                J = self.J[argJ]
                con = np.concatenate([I,J])
                
                #Define starting values forces
                fx = fy = 0
                
                #Loop over all contacts and compute total force
                for k in con:
                    #if k != k%self.ncon: print('warning')
                    k = int(k%len(self.fnor))
                                
                    #Compute components of forces
                    if k in J:
                        fx-=self.fnor[k]*self.nx[k]+self.ftan[k]*self.ny[k]
                        fy-=self.fnor[k]*self.ny[k]-self.ftan[k]*self.nx[k]
                    else:
                        fx+=self.fnor[k]*self.nx[k]+self.ftan[k]*self.ny[k]
                        fy+=self.fnor[k]*self.ny[k]-self.ftan[k]*self.nx[k]
                
                ftot = np.sqrt(np.square(fx)+ np.square(fy))

                
                if ftot > threshold:                  
                    #Compute angle needed for new contact
                    theta = (np.tan(fy/fx) + np.pi)%(2*np.pi)
                    angles = np.arctan2(self.x - self.x[i], self.y-self.y[i])
                    for k in range(len(angles)):
                        if angles[k] > theta + ang or angles[k] < theta - ang:
                            continue
                        if np.sqrt(np.square(self.x[k])+ np.square(self.y[k])) < dmax:
                            print("adding smart bond")
                            #Find correct index
                            j = int(self.id[k])
                            #Append to contacts
                            self.I = np.append(self.I, i)
                            self.J = np.append(self.J, j)

                            #Always add double bound
                            fullmobi_add.append(0)
                            
                            #Compute forces, this needs to be changed.
                            fnor0 = 1
                            ftan0 = 1
                            nx0 = 1
                            ny0 = 1
                            
                            #Add data to list
                            fnor_add.append(fnor0)
                            ftan_add.append(ftan0)
                            nx_add.append(nx0)
                            ny_add.append(ny0)
                
            
            #Add everything to existing arrays
            self.fnor=np.concatenate((self.fnor,np.array(fnor_add)))
            self.ftan=np.concatenate((self.ftan,np.array(ftan_add)))
            self.fullmobi=np.concatenate((self.fullmobi,np.array(fullmobi_add)))
            self.nx=np.concatenate((self.nx,np.array(nx_add)))
            self.ny=np.concatenate((self.ny,np.array(ny_add)))
            self.ncon=len(self.I)
            print("Added Smart Bonds")               
                    
            

        def AddNextBoundaryContacts(self,threshold=15,Brad=20.0):
            # Threshold to check if a particle is close enough to walls.
            upidx=np.argmax(self.ynext)
            downidx=np.argmin(self.ynext)
            leftidx=np.argmin(self.xnext)
            rightidx=np.argmax(self.xnext)
            
            # Boundary posiitons:
            # coordinates of virtual boundary particles: in the middle, one Brad off from the edge of the outermost particle
            up=self.ynext[upidx]
            yup = up+Brad+self.radnext[upidx]
            down=self.ynext[downidx]
            ydown = down-Brad-self.radnext[downidx]
            left=self.xnext[leftidx]
            xleft=left-Brad-self.radnext[leftidx]
            right=self.xnext[rightidx]
            xright=right+Brad+self.radnext[rightidx]
            
            # coordinates of virtual boundary particles: in the middle, one Brad off from the edge of the outermost particle
            Boundaries=np.zeros((4,3)) # Four boundary particles with their x,y and rad
            Boundaries[0,:]=[(left+right)*0.5,yup,Brad]
            Boundaries[1,:]=[(left+right)*0.5,ydown,Brad]
            Boundaries[2,:]=[xleft,(up+down)*0.5,Brad]
            Boundaries[3,:]=[xright,(up+down)*0.5,Brad]
            
            
            self.xnext=np.concatenate((self.xnext,Boundaries[:,0]))
            self.ynext=np.concatenate((self.ynext,Boundaries[:,1]))
            self.radnext=np.concatenate((self.radnext,Boundaries[:,2]))

            self.dx=self.xnext-self.x
            self.dy=self.ynext-self.y
            self.Nnext+=4
        
        #### ======================== Analysis helper functions     =================================================
        
        # computes basic contact, force, torque, and stress statistics
        def getStressStat(self):
                fsumx=np.empty((self.N,))
                fsumy=np.empty((self.N,))
                torsum=np.empty((self.N,))
                #------- those useful for plotting ----
                self.prepart=np.empty((self.N,))
                self.sxxpart=np.empty((self.N,))
                self.syypart=np.empty((self.N,))
                self.sxypart=np.empty((self.N,))
                self.syxpart=np.empty((self.N,))
                #-------------------------------------
                for u in range(self.N):
                        # problem - self.I doesn't seem to behave as an integer ...
                        # apparently this only works with arrays, not with lists??
                        c1=np.nonzero(np.array(self.I)==u)
                        c2=np.nonzero(np.array(self.J)==u)
                
                        fsumx[u]=np.sum(-self.fnor[c1[0]]*self.nx[c1[0]])+np.sum(self.fnor[c2[0]]*self.nx[c2[0]])+np.sum(self.ftan[c1[0]]*(-self.ny[c1[0]]))+np.sum(-self.ftan[c2[0]]*(-self.ny[c2[0]]))
                        fsumy[u]=np.sum(-self.fnor[c1[0]]*self.ny[c1[0]])+np.sum(self.fnor[c2[0]]*self.ny[c2[0]])+np.sum(self.ftan[c1[0]]*self.nx[c1[0]])+sum(-self.ftan[c2[0]]*self.nx[c2[0]])
                        torsum[u]=self.rad[u]*(np.sum(self.ftan[c1[0]])+np.sum(self.ftan[c2[0]]))
                        self.prepart[u]=np.sum(self.fnor[c1[0]]*self.rad[u])+np.sum(self.fnor[c2[0]]*self.rad[u])
                        self.sxxpart[u]=self.rad[u]*(np.sum(-self.fnor[c1[0]]*self.nx[c1[0]]*self.nx[c1[0]])+np.sum(self.fnor[c2[0]]*self.nx[c2[0]]*(-self.nx[c2[0]]))+np.sum(self.ftan[c1[0]]*(-self.ny[c1[0]])*(-self.ny[c1[0]]))+np.sum(-self.ftan[c2[0]]*(-self.ny[c2[0]])*(self.ny[c2[0]])))
                        self.syypart[u]=self.rad[u]*(np.sum(-self.fnor[c1[0]]*self.ny[c1[0]]*(self.ny[c1[0]]))+np.sum(self.fnor[c2[0]]*self.ny[c2[0]]*(-self.ny[c2[0]]))+np.sum(self.ftan[c1[0]]*(self.nx[c1[0]])*self.nx[c1[0]])+sum(-self.ftan[c2[0]]*self.nx[c2[0]]*(-self.nx[c2[0]])))
                        self.sxypart[u]=self.rad[u]*(np.sum(-self.fnor[c1[0]]*self.nx[c1[0]]*(-self.ny[c1[0]]))+np.sum(self.fnor[c2[0]]*self.nx[c2[0]]*(self.ny[c2[0]]))+np.sum(self.ftan[c1[0]]*(-self.ny[c1[0]])*(self.nx[c1[0]]))+np.sum(-self.ftan[c2[0]]*(-self.ny[c2[0]])*(-self.nx[c2[0]])))
                        self.syxpart[u]=self.rad[u]*(np.sum(-self.fnor[c1[0]]*self.ny[c1[0]]*(self.nx[c1[0]]))+np.sum(self.fnor[c2[0]]*self.ny[c2[0]]*(-self.nx[c2[0]]))+np.sum(self.ftan[c1[0]]*(self.nx[c1[0]])*(-self.ny[c1[0]]))+sum(-self.ftan[c2[0]]*self.nx[c2[0]]*(self.ny[c2[0]])))
                        del c1 
                        del c2
                
                # statistical output: don't do histograms here just yet; the problem of consistent and appropriate binning is too messy. 
                # return the whole thing.
                # except for the mobilization distribution which is unique
                # just get a few signature parameters here
                zav=2*len(self.I)/(1.0*self.N)
                nm=len(np.nonzero(self.fullmobi>0)[0])/(1.0*self.N)
                pres=np.sum(self.prepart)/(self.Lx*self.Ly)
                sxx=np.sum(self.sxxpart)/(self.Lx*self.Ly)
                syy=np.sum(self.syypart)/(self.Lx*self.Ly)
                sxy=np.sum(self.sxypart)/(self.Lx*self.Ly)
                syx=np.sum(self.syxpart)/(self.Lx*self.Ly)
                fxbal=np.mean(abs(fsumx))/np.mean(self.fnor)
                fybal=np.mean(abs(fsumy))/np.mean(self.fnor)
                torbal=np.mean(abs(torsum))/(np.mean(self.fnor)*np.mean(self.rad)) # correct units; do *not* normalize with ftan; consider a mu though
                mobin=np.linspace(-1,1,101)
                mohist,bin_edges=np.histogram(self.ftan/(self.mu*self.fnor),mobin)
                
                return zav,nm,pres,fxbal,fybal,torbal,mobin,mohist,sxx,syy,sxy,syx    
        
        #### =================== Computes the D2_min, amount of non-affine motion around a particle ==============
        def getD2min(self,threshold_rad):
                D2_min=np.zeros(len(self.x))
                for k in range(len(self.x)):
                    # Boundary is dubious
                    # the 10^5 is bad, but then they tend to be that value for experiment
                    # just a way to get the color scale to not f up
                    if (k in self.bindices):
                        D2_min[k]=1e5
                    else:
                        temp=[]
                        dist2=(self.x-self.x[k])**2+(self.y-self.y[k])**2
                        rad2=(self.rad[k]+np.mean(self.rad)+threshold_rad)**2
                        # yes, that includes myself, but that term drops out
                        temp=np.nonzero(dist2<rad2)
                        #for i in range(len(self.x)):
                                #if np.sqrt((self.x[k]-self.x[i])**2+(self.y[k]-self.y[i])**2)<(self.rad[k]+threshold_rad):
                                        #temp.append(i)
                        X=np.zeros((2,2))
                        Y=np.zeros((2,2))
                        epsilon=np.zeros((2,2))
                        if len(temp[0])>0:
                                
                                for neighbor in temp[0]:
                                        dx_next=self.xnext[neighbor]-self.xnext[k]
                                        dy_next=self.ynext[neighbor]-self.ynext[k]
                                        dx=self.x[neighbor]-self.x[k]
                                        dy=self.y[neighbor]-self.y[k]
                                        X[0,0]+=dx_next*dx
                                        X[0,1]+=dx_next*dy
                                        X[1,0]+=dy_next*dx
                                        X[1,1]+=dy_next*dy
                                        Y[0,0]+=dx*dx
                                        Y[0,1]+=dx*dy
                                        Y[1,0]+=dy*dx
                                        Y[1,1]+=dy*dy
                                epsilon[0,0]+=(X[0,0]/Y[0,0]+X[0,1]/Y[0,1])
                                epsilon[0,1]+=(X[0,0]/Y[1,0]+X[0,1]/Y[1,1])
                                epsilon[1,0]+=(X[1,0]/Y[0,0]+X[1,1]/Y[0,1])
                                epsilon[1,1]+=(X[1,0]/Y[1,0]+X[1,1]/Y[1,1])

                                for neighbor in temp[0]:
                                        dx_next=self.xnext[neighbor]-self.xnext[k]
                                        dy_next=self.ynext[neighbor]-self.ynext[k]
                                        dx=self.x[neighbor]-self.x[k]
                                        dy=self.y[neighbor]-self.y[k]
                                        D2_min[k]+=((dx_next- (epsilon[0,0]*dx+epsilon[0,1]*dy))**2+ (dy_next-(epsilon[1,0]*dx+epsilon[1,1]*dy))**2)
                        
                                                                    
                return D2_min
                
        
        # ====== Experimental or simulation displacements, decomposed into normal, tangential and potentiallt rotational displecements
        def Disp2Contacts(self,minThresh,debug=False):
                disp2n=np.zeros(self.ncon)
                disp2t=np.zeros(self.ncon)
                if self.hasAngles:
                    disp2r=np.zeros(self.ncon)
                    disp2gear=np.zeros(self.ncon)
                for k in range(self.ncon):
                        i=self.I[k]
                        j=self.J[k]
                        nx0=self.nx[k]
                        ny0=self.ny[k]
                        tx0=-self.ny[k]
                        ty0=self.nx[k]
                        disp2n[k] = ((self.dx[j]-self.dx[i])*nx0+(self.dy[j]-self.dy[i])*ny0)**2
                        disp2t[k] = ((self.dx[j]-self.dx[i])*tx0+(self.dy[j]-self.dy[i])*ty0)**2
                        if self.hasAngles:
                            disp2r[0]=(self.rad[j]*self.dalpha[j]-self.rad[i]*self.dalpha[i])**2
                            disp2gear[0]=((self.dx[j]-self.dx[i])*tx0+(self.dy[j]-self.dy[i])*ty0 - (self.rad[i]*self.dalpha[i]+self.rad[j]*self.dalpha[j]))**2
                thresh=np.mean(disp2n+disp2t)
                print ("Internally generated threshold of " + str(thresh))
                if thresh<minThresh:
                    thresh=minThresh
                if self.hasAngles:
                    return disp2n, disp2t, disp2r, disp2gear,thresh
                else:
                    return disp2n, disp2t, thresh
            
            
        ##================ Finally, pure helper functions: positions of both ends of a contact ===============
        # get the cleaned up, periodic boundary conditions sorted out positions corresponding to two ends of a contact. 
        # Basic plotting helper function
        def getConPos(self,k):
                x0=self.x[self.I[k]]
                x1=self.x[self.J[k]]
                y0=self.y[self.I[k]]
                y1=self.y[self.J[k]]
                if self.periodic:
                        x1=x1-self.Lx*np.round((x1-x0)/self.Lx)
                        yover=np.round((y1-y0)/self.Ly)
                        if (yover!=0):
                                y1=y1-self.Ly*yover
                                x1-=self.Lx*self.strain*yover
                                x1=x1-self.Lx*np.round((x1-x0)/self.Lx)
                if self.addBoundarySquare:
                    ival=self.I[k]
                    if ((ival==self.bindices[0]) or (ival==self.bindices[1])): #top or bottom
                        x0=x1
                    if ((ival==self.bindices[2]) or (ival==self.bindices[3])): #left or right
                        y0=y1
                if self.addBoundaryAnnulus:
                    ival=self.I[k]
                    l = 100 #Length in pixels of the contacts with the boundary
                    #If in contact with the inner boundary
                    if (ival==self.bindices[0]):
                        #Compute position vector from midpoint
                        r = np.array([x1-self.mid_x, y1-self.mid_y])
                        #Compute length of r
                        r_len = np.linalg.norm(r)
                        #Compute angle between position vector from midpoint and (1,0) vector
                        ang = np.arctan2(r[1], r[0])
                        #Compute coordinates of outward inward coordinates
                        x0 = self.mid_x + (r_len - l)*np.cos(ang)
                        y0 = self.mid_y + (r_len - l)*np.sin(ang)
                    #If in contact with the outer boundary 
                    if (ival==self.bindices[1]):
                        #Compute position vector from midpoint
                        r = np.array([x1-self.mid_x, y1-self.mid_y])
                        #Compute length of r
                        r_len = np.linalg.norm(r)
                        #Compute angle between position vector from midpoint and (1,0) vector
                        ang = np.arctan2(r[1], r[0])
                        #Compute coordinates of outward pointing coordinates
                        x0 = self.mid_x + (r_len + l)*np.cos(ang)
                        y0 = self.mid_y + (r_len + l)*np.sin(ang)
                return x0,x1,y0,y1
        
        # same, but based on existing particle labels (in case those come from elsewhere)
        def getConPos2(self,k1,k2):
            x0=self.x[k1]
            x1=self.x[k2]
            y0=self.y[k1]
            y1=self.y[k2]   
            if self.periodic:
                    x1=x1-self.Lx*np.round((x1-x0)/self.Lx)
                    yover=np.round((y1-y0)/self.Ly)
                    if (yover!=0):
                            y1=y1-self.Ly*yover
                            x1-=self.Lx*self.strain*yover
                            x1=x1-self.Lx*np.round((x1-x0)/self.Lx)
            if self.addBoundarySquare:
                if ((k1==self.bindices[0]) or (k1==self.bindices[1])): #top or bottom
                    x0=x1
                if ((k1==self.bindices[2]) or (k1==self.bindices[3])): #left or right
                    y0=y1
            if self.addBoundaryAnnulus:
                l = 100 #Length in pixels of the contacts with the boundary
                if (k1==self.bindices[0] and k2==self.bindices[1]) or (k1==self.bindices[1] and k2==self.bindices[0]):
                    return x0,x1,y0,y1
                #If in contact with the inner boundary
                if (k1==self.bindices[0]):
                    #Compute position vector from midpoint
                    r = np.array([x1-self.mid_x, y1-self.mid_y])
                    #Compute length of r
                    r_len = np.linalg.norm(r)
                    #Compute angle between position vector from midpoint and (1,0) vector
                    ang = np.arctan2(r[1], r[0])
                    #Compute coordinates of outward inward coordinates
                    x0 = self.mid_x + (r_len - l)*np.cos(ang)
                    y0 = self.mid_y + (r_len - l)*np.sin(ang)
                #If in contact with the outer boundary 
                if (k1==self.bindices[1]):
                    #Compute position vector from midpoint
                    r = np.array([x1-self.mid_x, y1-self.mid_y])
                    #Compute length of r
                    r_len = np.linalg.norm(r)
                    #Compute angle between position vector from midpoint and (1,0) vector
                    ang = np.arctan2(r[1], r[0])
                    #Compute coordinates of outward pointing coordinates
                    x0 = self.mid_x + (r_len + l)*np.cos(ang)
                    y0 = self.mid_y + (r_len + l)*np.sin(ang)
            return x0,x1,y0,y1