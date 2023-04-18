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
        
    #Plotting the graph
    def graph(self):
        if isinstance(self.I, int):
            print('no data')
        else:
            for i in self.I:
                for k in self.J[self.I==i]:
                    x0, x1, y0, y1 = self.conf.getConPos2(i, k)
                    plt.plot([x0, x1], [y0, y1], color='orange', marker=',')
            plt.show()
    
    #Plotting the Maxwell Cremona tiling
    def tile(self):
        if isinstance(self.I, int):
            print('no data')
        else:
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
                
                
                #Plot
                sumx = sumy = 0
                for k in range(len(arg)):
                    theta = np.arctan2(self.y[i]-self.y[con[k]], self.x[i]-self.x[con[k]])
                    if self.ft[arg[k]] > 0:
                        fx = self.fn[arg[k]]*np.cos(theta)-self.ft[arg[k]]*np.sin(theta)
                        fy = self.fn[arg[k]]*np.sin(theta)+self.ft[arg[k]]*np.cos(theta)
                    else:
                        fx = self.fn[arg[k]]*np.cos(theta)+np.abs(self.ft)[arg[k]]*np.sin(theta)
                        fy = self.fn[arg[k]]*np.sin(theta)-np.abs(self.ft)[arg[k]]*np.cos(theta)

                    print(i, k, fx, fy)
                                    
                    sumx += fx
                    sumy += fy
                    
                    #plt.arrow(xor, yor, fx, fy, width=.08, facecolor=color[i])
                    plt.plot([xor, xor+fx], [yor, yor+fy], color=color[i], marker='.')
                    xor += fx
                    yor += fy
                print('sum', sumx, sumy)
            plt.show()                
    