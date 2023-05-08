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
    def plotter(self, data, xor, yor, color, zorder, ls, alpha, arrow, plot):
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
            if plot:
                if arrow: plt.arrow(xor, yor, fx, fy, color=color, zorder=zorder, ls=ls, alpha=alpha)
                else: plt.plot([xor, xor+fx], [yor, yor+fy], color=color, marker='o', zorder=zorder, ls=ls, alpha=alpha)
            
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
    
    
    #Plotting the Maxwell Cremona tiling
    def tile(self, arrow):
        if isinstance(self.I, int):
            print('no data')
        else:
            #Initial values
            xor1 = yor1 = l = 0 #Staring position x, starting position y, starting angle, counter
            checklist0 = np.unique(np.union1d(self.I, self.J)) #Checklist that does not change
            checklist = np.unique(np.union1d(self.I, self.J)) #Checklist for checking if all contacts are plotted
            contact = np.max(checklist)+1  #Start with an element that is for sure not in the checklist
            a = True #So we enter the first while loop
            i = checklist[0] #Start with first entry of the array, choosing different starting point should not change result
            priority = False #value used for searching connected component
            plot = True #Plot the tiles if true
            
            #Loop over all particles that have contacts
            while a:
                if len(checklist) == 0: break #if we have checked all contacts then exit
                if l >= 1: #If we are not in the first iteration check over all contacts
                    print(i, checklist)
                    #If we dont have any new contacts to plot
                    if len(contactlist)==0: 
                        #Set a goal
                        inew=checklist[0]
                        priority = True
                        
                        #Start over to ensure correct positions
                        i = checklist0[0]
                        xor1 = yor1 = 0 
                    else:
                        #Search for new contact
                        for n in range(len(contactlist)):
                            i = contactlist[n]
                            if i in checklist: break #If we found a new contact in the contactslist break out of the loop
                            
                            #If a new contact can not be found in the contactlist, get a new entry from the checklist 
                            if n == len(contactlist)-1: 
                                i=checklist[0]
                        
                        #Search for the priority 
                        if priority:
                            for n in range(len(contactlist)):
                                if contactlist[n]==inew: 
                                    i = contactlist[n]
                                    priority = False
                                    print('priority')
                                    break
                            
                        #Retrieve previous contact to get correct origin position
                        for n in range(len(con)):
                            arr = con[n]
                            if i in arr[1]:
                                contact = arr[0]
                                break
                        
                        #Retrieve new angle from previous tiles
                        for n in range(len(data)):
                            arr = data[n]
                            x = np.argwhere(arr[:,0]==i).flatten()
                            #If we have an entry equal to i and this is in contact with the previous entry
                            if len(x) != 0 and arr[x].flatten()[-1]==contact:
                                phi1 = arr[x].flatten()[1] #Retrieve angle 
                                break
                        
                        #Retrieve starting position from previous tiles                 
                        arr = orr[n]
                        x = np.argwhere(arr==i).flatten()
                        if len(x)!=0: 
                            x=x[0]
                            n = (x+1)%len(arr)
                            xor1 = arr[n,1]
                            yor1 = arr[n,2]
                        else:
                            print('warning, in a new cluster')
                            xor1=yor1=1
            
                #Remove contact from the checklist
                checklist = checklist[checklist != i]
                
                #Counter and value for next while loop
                k = 0
                b = True
                
                #Get contact data and force data for contact i
                arg1, con1, data1 = self.contact(i, i=contact)
                
                #Plot and get force vector data
                orr1 = self.plotter(data1, xor1, yor1, color='black', zorder=1, ls=':', alpha=1, arrow=arrow, plot=plot)
                
                #Create empty list for all contacts of contacts
                contactlist = []
                
                #Create lists to store data from which origin and angle can be recovered
                data = []
                orr = []
                con = []
                
                #Loop over all contacts
                while b:
                    if (con1[k] in checklist) or priority: plot=True
                    #Generate colors for each tile
                    colors = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)]) for i in range(len(con1))]
                    
                    #Get contact data
                    arg2, con2, data2 = self.contact(con1[k], i=i)
                    #print(i, con1, con1[k], data2)
                    
                    #Create list of all the contacts of the contacts
                    contactlist.extend(con2.tolist())
                    contactlist = [*set(contactlist)] #remove duplicates
                    
                    #Get origin coordinates
                    x = np.argwhere(orr1[:,0]==con1[k]).flatten()[0]
                    n = (x+1)%len(con1)
                    xor2 = orr1[n,1]
                    yor2 = orr1[n,2]
                    
                    #Plot the result and get origin coordinates
                    orr1 = self.plotter(data1, xor1, yor1, color='black', zorder=1, ls=':', alpha=1, arrow=arrow, plot=plot)
                    orr2 = self.plotter(data2, xor2, yor2, color=colors[k], zorder=0, ls='-', alpha=0.7, arrow=arrow, plot=plot)
                    
                    #Save data
                    con.append([con1[k],con2])
                    orr.append(orr2)
                    data.append(data2)

                    #Remove checked contacts from the checklist
                    checklist = checklist[checklist != con1[k]]
                    
                    #Move to next contact
                    k+=1
                    
                    #If all contacts have been checked, break of the while loop
                    if k >= len(con1):
                        b = False #This lets us break from the while loop
                        l += 1 #Counter for making sure that we are not in the first iteration
            
            plt.title("Maxwell-Cremona Tiling")
            plt.show()
            
    
    


                
                                
    