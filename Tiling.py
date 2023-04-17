#Import necessary classes
from Configuration import *

#Import necessary packages
import pickle as pickle
import copy as cp
import numpy as np
import matplotlib.pyplot as plt


class Tiling:
    
    def __init__(self, conf, scale):
        self.conf = conf
        self.I = conf.I 
        self.J = conf.J
        self.scale = scale
        
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
    