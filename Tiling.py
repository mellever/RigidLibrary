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
    
    def plotter(self, data, xor, yor, color):
        orr = np.zeros((len(data[:,0]),3))
        #Plot the tile
        for k in range(len(data[:,0])):
            orr[k, 0] = data[k,0]
            orr[k, 1] = xor
            orr[k, 2] = yor
            fx = data[k,2]
            fy = data[k,3]
            plt.arrow(xor, yor, fx, fy, color=color, alpha=0.5)
            #plt.plot([xor, xor+fx], [yor, yor+fy], color=color, marker='o')
            xor+=fx
            yor+=fy
        return orr
    
    def contact(self, i, flip):
        #Determine indices of contacts and put together in one array
        argI = np.argwhere(self.J==i).flatten()
        argJ = np.argwhere(self.I==i).flatten()
        arg = np.concatenate([argI,argJ])
        #Get particle ids from indices and put together in one array
        I = self.I[argI]
        J = self.J[argJ]
        con = np.concatenate([I,J])
    
        #Force data we want to extract
        data=np.zeros((len(arg),4))
        
        #Loop over all contacts 
        for k in range(len(arg)):
            #Add particle number
            data[k,0] = con[k]
            
            #Compute angle between particles
            theta = np.arccos(self.nx[arg[k]])
            
            #Compute components of forces
            fx=self.fn[arg[k]]*self.nx[arg[k]]+self.ft[arg[k]]*self.ny[arg[k]]
            fy=self.fn[arg[k]]*self.ny[arg[k]]-self.ft[arg[k]]*self.nx[arg[k]]
            if con[k] in J:
                fx=-fx
                fy=-fy
                
     
            #Add to array
            data[k,2] = fx
            data[k,3] = fy
            data[k,1] = theta

        #Sort array from smallest to largest angle
        data = data[data[:, 1].argsort()]
        
        #If in con2 flip the data
        #if flip: data = np.flipud(data)
        
        return arg, con, data
    
    
    #Plotting the Maxwell Cremona tiling
    def tile(self, verbose0=True):
        if isinstance(self.I, int):
            print('no data')
        else:
            #Initial values      
            xor = yor = 0
            a = True
            flip = True
            checklist = np.unique(np.union1d(self.I, self.J))
            l = 0
            i = checklist[0]
            
            #Loop over all particles that have contacts
            while a:
                xor = yor = 0
                if len(checklist) == 0: break #if we have checked all contacts then exit
                if l >= 1: #If we are not in the first iteration check over all contacts
                    for n in range(len(contactlist)):
                        i = contactlist[n]
                        if i in checklist:
                            break
                        if n == len(contactlist)-1: i=checklist[0]
                
                #Remove contact from the checklist
                checklist = checklist[checklist != i]
                
                #Counter and value for next while loop
                k = 0
                b = True
                
                #Get contact data and force data for contact i
                arg1, con1, data1 = self.contact(i, False)
                print(data1)
                
                #Plot and get force vector data
                orr = self.plotter(data1, xor, yor, color='black')
                print(orr)
                
                #Create empty list for all contacts of contacts
                contactlist = []
                
                #Loop over all contacts
                while b:
                    #Generate colors for each tile
                    colors = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)]) for i in range(len(con1))]
                    
                    #Get contact data
                    arg2, con2, data2 = self.contact(con1[k], flip=True)
                    
                    #Create list of all the contacts of the contacts
                    contactlist.extend(con2.tolist())
                    contactlist = [*set(contactlist)] #remove duplicates
                    
                    #Get origin coordinates
                    xor = orr[orr[:,0]==con1[k]][0][1]
                    yor = orr[orr[:,0]==con1[k]][0][2]
                    print(i, con1, con1[k], con2, xor, yor)
                    
                    self.plotter(data1, 0, 0, color='black')
                    self.plotter(data2, xor, yor, color=colors[k])
                    checklist = checklist[checklist != con1[k]]
                    k+=1
                    if k >= len(con1):
                        b = False
                        l += 1
                    plt.show()
                #plt.show()
    
    


                
                                
    