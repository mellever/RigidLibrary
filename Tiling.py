#Import necessary classes
from Configuration import *

#Import necessary packages
import pickle as pickle
import copy as cp
import numpy as np
import matplotlib.pyplot as plt


class Tiling:
    
    def __init__(self, conf):
        self.conf = conf
        self.I = conf.I 
        self.J = conf.J
        self.x = conf.x
        self.y = conf.y
        self.fn = conf.fnor
        self.ft = conf.ftan
        self.nx = conf.nx
        self.ny = conf.ny
        
    #Plotting the graph
    def graph(self, verbose0):
        if isinstance(self.I, int):
            print('no data')
        else:
            for i in self.I:
                for k in self.J[self.I==i]:
                    x0, x1, y0, y1 = self.conf.getConPos2(i, k)
                    plt.plot([x0, x1], [y0, y1], color='black', marker='o', markersize=15, markerfacecolor='white')
                    if verbose0: plt.annotate(k, (x1,y1))
                if verbose0: plt.annotate(i, (x0,y0))
            plt.title("Contact Network")
            plt.show()
    
    #Plotting the Maxwell Cremona tiling
    def tile(self, verbose0=True):
        if isinstance(self.I, int):
            print('no data')
        else:
            #print(self.fn)
            #print(self.ft)
            #Generate colors for plotting
            color = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)]) for i in range(np.max(self.I)+1)]
            for i in self.I:
                #Determine indices of contacts and put together in one array
                argI = np.argwhere(self.J==i).flatten()
                argJ = np.argwhere(self.I==i).flatten()
                arg = np.concatenate([argI,argJ])
                #Get particle ids from indices and put together in one array
                I = self.I[argI]
                J = self.J[argJ]
                con = np.concatenate([I,J])
                
                #if we have less than three contacts move to the next one
                if len(con) < 3:
                    continue
                
                #Start position of tile, this needs to be modified
                xor = (self.x[i]+self.x[con[0]]+self.x[con[1]])/3
                yor = (self.y[i]+self.y[con[0]]+self.y[con[1]])/3
                xor = yor = 0
            
                #Force data we want to extract
                data=np.zeros((len(arg),4))
                
                #Loop over all contacts 
                for k in range(len(arg)):
                    #Add particle number
                    data[k,0] = con[k]
                    
                    #Compute angle between particles
                    theta = np.arccos(self.nx[arg[k]])
                    data[k,1] = theta
                    
                    #Compute components of forces
                    fx=self.fn[arg[k]]*self.nx[arg[k]]+self.ft[arg[k]]*self.ny[arg[k]]
                    fy=self.fn[arg[k]]*self.ny[arg[k]]-self.ft[arg[k]]*self.nx[arg[k]]
                    data[k,2] = fx
                    data[k,3] = fy
        
                #Sort array from smallest to largest angle
                data = data[data[:, 1].argsort()]
                
                #Compute sum of the components of the sum for debugging
                sumx = np.sum(data[:,2])
                sumy = np.sum(data[:,3])
                print('sum', sumx, sumy)
                if np.abs(sumx) < 1e-5:
                    print(i, con)
                
                #Plot the tile
                for k in range(len(data[:,0])):
                    fx = data[k,2]
                    fy = data[k,3]
                    plt.plot([xor, xor+fx], [yor, yor+fy], color=color[i])
                    xor+=fx
                    yor+=fy
                plt.show()


                
                                
    