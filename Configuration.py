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
from scipy.spatial import ConvexHull
from scipy.spatial.distance import cdist



class Configuration:

        #=================== Create a new configruation ===============
        # Current choices: either simulation or experimental through datatype. Parameters are either read (simulation) or passed to the configuration via optional arguments (experiment)
        def __init__(self, folder, datatype, mu0, strainstep):
            self.folder = folder
            self.datatype = datatype
            self.addBoundary = False
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
            elif (datatype == 'experiment'):
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
                self.rad = sp.loadtxt(folder + 'posinit.dat', usecols=(3,))

        #######========== Simulation data read-in ==================
        # snap is the label, and distSnapshot tells me how many time units they are apart (in actual time, not steps)
        def readSimdata(self, snap, verbose0, distSnapshot0=100.0):
                self.verbose = verbose0
                self.distSnapshot = distSnapshot0
                self.strain = self.gammadot*(1.0*snap)*self.distSnapshot
                if (self.verbose):
                        print('L=' + str(self.L))
                        print('Np=' + str(self.N))
                        print('mu=' + str(self.mu))
                        print('xi=' + str(self.xi))
                        print('phi=' + str(self.phi))
                        print('strain=' + str(self.strain))
                filename = foldername + '/PosVel' + str(snap) + '.dat'
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
                        print('error: no position data at snapshot' + str(snap))
                        self.x = 0.0
                        self.y = 0.0
                        self.dx = 0.0
                        self.dy = 0.0
                        self.alpha = 0.0
                        self.N = 0
                else:
                        self.x = data[0, :]
                        self.y = data[1, :]
                        self.alpha = data[2, :]
                        self.dx = data[3, :]
                        self.dy = data[4, :]
                        self.dalpha = data[5, :]
                        del data

                filename = foldername + '/Contact' + str(snap) + '.dat'
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
                        print('error: no position data at snapshot' + str(snap))
                        self.noConf = True
                        self.I = -1
                        self.J = -1
                        self.fullmobi = 0
                        self.Ifull = -1
                        self.Jfull = -1
                else:
                        self.noConf = False
                        self.I = list(data[0, :].astype(int))
                        self.J = list(data[1, :].astype(int))
                        self.fullmobi = data[4, :].astype(int)
                        self.nx = data[2, :]
                        self.ny = data[3, :]
                        self.fnor = data[5, :]
                        self.ftan = data[6, :]+data[7, :]
                        del data
                        self.x -= self.L*np.round(self.x/self.L)
                        self.y -= self.L*np.round(self.y/self.L)
                        self.ncon = len(self.I)

        #========== Experimental data read-in ==================
        def ReadExpdata(self, verbose):
                prefix = self.folder +'/particle_positions.txt'
                try:
                    pos_data = pandas.read_csv(prefix, names=['id', 'x', 'y', 'r', 'n'])
                    pos_data = pos_data[pos_data['n'] == self.step]
                    self.x = pos_data.loc[:,'x'].to_numpy()
                    self.y = pos_data.loc[:,'y'].to_numpy()
                    self.rad = pos_data.loc[:,'r'].to_numpy()
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
                    con_data = pandas.read_csv(prefix, names=['n', 'id1', 'id2', 'ft', 'fn'])
                    con_data= con_data[con_data['n'] == self.step]
                except:
                    print("Error: there is no contact data here")
                    return 1
                
                #Create a list of entries that contain a duplicate
                tups = list(zip(con_data.id1, con_data.id2))
                dups_bool = [None]*len(tups)                       
                for i in range(len(tups)):
                    for j in tups:
                        if (tups[i] != j and tups[i][0] == j[1] and tups[i][1] == j[0]):
                            dups_bool[i] = True   
                
                #Add dups_bool array to the last column of the data frame
                con_data['conbool'] = dups_bool
                con_frame = con_data #This step is to supress copy warning, since adding column creates a copy. 
                
                #Create lists
                self.I=[]
                self.J=[]
                fn0=[]
                ft0=[]
                fm0=[]
                
                #Filter entries and determine if slipping or not
                for k in range(len(dups_bool)):
                    try:
                        if con_frame['conbool'][k] == 'False':
                                self.I.append(con_frame['id1'][k])
                                self.J.append(con_frame['id2'][k])
                                fn0.append(con_frame['fn'][k])
                                ft0.append(con_frame['ft'][k])
                                if (abs(con_frame['ft'][k])/con_frame['fn'][k]>self.mu):
                                    fm0.append(1)
                                else:
                                    fm0.append(0)
                        else:
                            self.I.append(con_frame['id1'][k])
                            self.J.append(con_frame['id2'][k])
                            
                            # Search other force and take the mean
                            norm2 = con_frame[ (con_frame['id1']==con_frame['id2'][k]) & (con_frame['id2']==con_frame['id1'][k])]['fn'].iloc[0]
                            tan2 = con_frame[ (con_frame['id1']==con_frame['id2'][k]) & (con_frame['id2']==con_frame['id1'][k])]['ft'].iloc[0]
                            norm = (con_frame['fn'][k] + norm2)/2
                            tan = (con_frame['fn'][k] + tan2)/2
                            
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
        def AddBoundaryContacts(self,threshold=20,Brad=30.0):
            self.addBoundary=True

            # Amount of contacts before adding boundaries
            self.ncon3 = self.N
            
            # Fist we need to determine the radius of the outer ring
            # This is done by finding the maximum distance between points
            points = np.vstack((self.x, self.y)).T
            # Compute convex hull
            hull = ConvexHull(points)
            # Boundary of convex hull
            hullpoints = points[hull.vertices,:]
            dmax = cdist(hullpoints, hullpoints, metric='euclidean')
            
            # Get the farthest apart points
            dmax_idx = np.unravel_index(dmax.argmax(), dmax.shape)
            #Get coordinates maximal points
            max_points = [[hullpoints[dmax_idx[0]][0], hullpoints[dmax_idx[0]][1]],[hullpoints[dmax_idx[1]][0], hullpoints[dmax_idx[1]][1]]]

            #Get indidices of maximal points
            idx1 = np.where(self.x == max_points[0][0])[0]
            idx1 = idx1[np.where(self.y[idx1] == max_points[0][1])[0][0]]
            idx2 = np.where(self.x == max_points[1][0])[0]
            idx2 = idx2[np.where(self.y[idx2] == max_points[1][1])[0][0]]
            
            #Use the above indices to find radi of particles
            R1_outer = self.rad[idx1]
            R2_outer = self.rad[idx2]

            # Compute distance between these points and add radi
            D_outer = math.dist(max_points[0], max_points[1]) + R1_outer + R2_outer
            
            #Compute radius for outer boundary particles:
            R_outer = D_outer/2 + Brad/2 #Check this!!

            #Determine midpoint of the circle
            xm = (max_points[0][0]+ max_points[1][0])/2
            ym = (max_points[0][1]+ max_points[1][1])/2
            
            #Boundary particles onto outer boundary
            num_outer = 200 #Amount of boundary particles on outer radius
            theta = np.linspace(0, 2*np.pi, num_outer, endpoint=False)
            outer_b = np.zeros((len(theta), 3))
            for i in range(len(outer_b)):
                outer_b[i,:] = [xm + R_outer*np.cos(theta[i]), ym + R_outer*np.sin(theta[i]), Brad]

            Boundaries = outer_b.astype(int) #Approximate boundary coordinates by ints
            
            # Find the particles in contact with the boundary, and label correctly
            self.bindices=np.zeros(len(theta)).astype(int)
            for i in range(len(self.bindices)):
                self.bindices[i] = self.N+i

            fullmobi_add=[]
            fnor_add=[]
            ftan_add=[]
            nx_add=[]
            ny_add=[]

            for i in range(len(Boundaries)):
                for j in range(self.N):
                    d = np.sqrt((self.x[j]-Boundaries[i][0])**2 + (self.y[j]-Boundaries[i][1])**2) - self.rad[j] - Brad
                    if d <= threshold:
                        self.I = np.append(self.I, self.N + i) ##check this
                        self.J = np.append(self.J, j)
                        
                        neii=np.nonzero(self.I[:self.ncon]==j)[0]
                        neij=np.nonzero(self.J[:self.ncon]==j)[0]

                        # This needs to change!!!
                        nx0 = ny0 = 1

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
            self.N+= len(theta)
            #b = 0
            #print(self.I[b], self.x[self.I[b]], self.y[self.I[b]])
            #print(self.J[b], self.x[self.J[b]], self.y[self.J[b]])
            print ("Added boundaries!")      

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
                if self.addBoundary:
                    ival=self.I[k]
                    if ((ival==self.bindices[0]) or (ival==self.bindices[1])): #top or bottom
                        x0=x1
                    if ((ival==self.bindices[2]) or (ival==self.bindices[3])): #left or right
                        y0=y1
                return x0,x1,y0,y1
        
        # same, but based on existing particle labels (in case those come from elsewhere)
        def getConPos2(self,k1,k2):
            if self.experiment and k1 <= self.ncon3 and k2 <= self.ncon3:
                x0=self.x[k1-1]
                x1=self.x[k2-1]
                y0=self.y[k1-1]
                y1=self.y[k2-1]
            elif self.experiment and k1 > self.ncon3 and k2 <= self.ncon3:
                x0=self.x[k1-1]
                x1=self.x[k2]
                y0=self.y[k1-1]
                y1=self.y[k2]
            elif self.experiment and k1 <= self.ncon3 and k2 > self.ncon3:
                x0=self.x[k1]
                x1=self.x[k2-1]
                y0=self.y[k1]
                y1=self.y[k2-1]
            else:
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
            if self.addBoundary:
                if ((k1==self.bindices[0]) or (k1==self.bindices[1])): #top or bottom
                    x0=x1
                if ((k1==self.bindices[2]) or (k1==self.bindices[3])): #left or right
                    y0=y1
            return x0,x1,y0,y1
        
