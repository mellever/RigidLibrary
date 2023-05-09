#Import necessary classes
from Configuration import *

#Import necessary packages
import pickle as pickle
import copy as cp
import numpy as np
import matplotlib.pyplot as plt
from queue import Queue


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
        
    #Plotting the contact network
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
    
    #Function that plots tile and returns plotted points
    def plotter(self, data, xor, yor, color, zorder, ls, arrow):
        #Create array for saving the data
        orr = np.zeros((len(data[:,0]),3))
        #Loop over the data array
        for k in range(len(data[:,0])):
            #Save data
            orr[k, 0] = data[k,0]
            orr[k, 1] = xor
            orr[k, 2] = yor
            
            #Retrieve forces
            fx = data[k,2]
            fy = data[k,3]
            
            #Plotting
            if arrow: plt.arrow(xor, yor, fx, fy, color=color, zorder=zorder, ls=ls)
            else: plt.plot([xor, xor+fx], [yor, yor+fy], color=color, marker='o', zorder=zorder, ls=ls)
            
            #Move to next point on the tile
            xor+=fx
            yor+=fy
        return orr
    
    def contact(self, contact, i):
        #Determine indices of contacts and put together in one array
        argI = np.argwhere(self.J==contact).flatten()
        argJ = np.argwhere(self.I==contact).flatten()
        arg = np.concatenate([argI,argJ])
        #Get particle ids from indices and put together in one array
        I = self.I[argI]
        J = self.J[argJ]
        con = np.concatenate([I,J])
    
        #Force data we want to extract
        data=np.zeros((len(arg),5))
        
        #Loop over all contacts 
        for k in range(len(arg)):
            #Add particle number
            data[k,0] = con[k]
            
            #Add contact
            data[k,4] = contact
            
            #Directions and forces
            nx = self.nx[arg[k]]
            ny = self.ny[arg[k]]
            fn = self.fn[arg[k]]
            ft = self.ft[arg[k]]
            
            #Make sure that contacts have equal but opposite forces
            if con[k] in I:
                nx = -nx
                ny = -ny
            
            #Compute angle between particles
            theta = np.arctan2(ny, nx) #%(2*np.pi)
            if theta < 0: 
                theta+=2*np.pi
            
            #Compute components of forces
            fx=fn*nx+ft*ny
            fy=fn*ny-ft*nx
     
            #Add to array
            data[k,2] = fx
            data[k,3] = fy
            data[k,1] = theta

        #Sort array from smallest to largest angle
        data = data[data[:, 1].argsort()]
        
        #Permute array if necessary
        startidx = np.argwhere(data[:,0]==i).flatten()
        if len(startidx)==1:
            ang = data[startidx[0], 1]
            data[:,1] = (data[:,1] - ang)%(2*np.pi)
            data = data[data[:, 1].argsort()]
        
        return arg, con, data
    
    #Function that returns all the contacts of a given particle
    def adjacency(self, i):
        argI = np.argwhere(self.J==i).flatten()
        argJ = np.argwhere(self.I==i).flatten()
        arg = np.concatenate([argI,argJ])
        #Get particle ids from indices and put together in one array
        I = self.I[argI]
        J = self.J[argJ]
        con = np.concatenate([I,J]).tolist()
        return con
    
    #Function that creates adjacency list dictionary
    def adjacency_list(self, checklist):
        self.adj_list = {}
        for i in checklist:
            self.adj_list[i] = self.adjacency(i)
        
    
    #Function that performs breadth first search, such that all contacts get plotted
    #This has been taken from https://towardsdatascience.com/introduction-to-graph-algorithm-breadth-first-search-algorithm-in-python-8644b6d31880
    #Currently visited and level are not in use, so these can be deleted
    def BFS(self, s):
        visited = {}
        level = {}
        parent = {}
        traversal_output = []
        queue = Queue()
        for node in self.adj_list.keys():
            visited[node] = False
            parent[node] = None
            level[node] = -1
        visited[s] = True
        level[s] = 0
        queue.put(s)
        while not queue.empty():
            u = queue.get()
            traversal_output.append(u)
            for v in self.adj_list[u]:
                if not visited[v]:
                    visited[v] = True
                    parent[v] = u
                    level[v] = level[u] + 1
                    queue.put(v)
        return traversal_output, visited, level, parent
    
    
    #Function for maxwell cremona tiling
    def tile(self, arrow):
        if isinstance(self.I, int):
            print('no data')
        else:
            #Stating values
            xor1 = yor1 = l = 0 #Staring position x, starting position y, counter for amount of tiles
            checklist = np.unique(np.union1d(self.I, self.J)) #Checklist for checking if all contacts are plotted
            s = checklist[0] #Starting vertex
            contact = np.max(checklist)+1  #Start with an element that is for sure not in the checklist
            
            #Generate colors for each tile
            colors = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)]) for i in range(len(checklist))]
            
            #Create lists to store data from which origin and angle can be recovered
            orr = []

            #Create adjacency dictionary
            self.adjacency_list(checklist)
            
            #Perform breadth first search
            traversal_output, visited, level, parent = self.BFS(s)
            
            for i in traversal_output:
                #If all contacts have been plotted stop the code
                if len(checklist)==0: break
                
                #Only execute when a contact has not been plotted
                if i not in checklist: continue
                   
                #Remove contact from the checklist
                checklist = checklist[checklist != i]    

                #If we are not in the first iteration
                if l>=1:
                    #Get the contact via the parent node
                    contact = parent[i]
                    
                    #Extract the origin coordinates
                    for n in range(len(orr)):
                        arr = orr[n]
                        if arr[0]==contact:
                            arr1 = arr[1]
                            x = np.argwhere(arr1[:,0]==i).flatten()[0]
                            n = (x+1)%len(arr1)
                            xor1 = arr1[n,1]
                            yor1 = arr1[n,2]
                            break
               
                #Get contact data and force data for contact i
                arg1, con1, data1 = self.contact(i, i=contact)
                
                #Plot and get force vector data
                orr1 = self.plotter(data1, xor1, yor1, color=colors[l], zorder=0, ls='-', arrow=arrow)
                
                #Counter for amount of tiles
                l+=1
                
                #Loop over all contacts
                for k in range(len(con1)):
                    #Skip is necesarry
                    if con1[k] not in checklist: continue
                    
                    #Get contact data
                    arg2, con2, data2 = self.contact(con1[k], i=i)
                    
                    #Get origin coordinates
                    x = np.argwhere(orr1[:,0]==con1[k]).flatten()[0]
                    n = (x+1)%len(con1)
                    xor2 = orr1[n,1]
                    yor2 = orr1[n,2]
                    
                    #Plot the result and get origin coordinates
                    #orr1 = self.plotter(data1, xor1, yor1, color='black', zorder=1, ls=':', arrow=arrow)
                    orr2 = self.plotter(data2, xor2, yor2, color=colors[l], zorder=0, ls='-', arrow=arrow)
                    
                    #Save data
                    orr.append([con1[k],orr2])
                    
                    #Remove contact from the checklist
                    checklist = checklist[checklist != con1[k]]
                    
                    #Counter for amount of tiles
                    l+=1
              
            plt.title("Maxwell-Cremona Tiling")
            plt.show()
                
                
            
    
    


                
                                
    