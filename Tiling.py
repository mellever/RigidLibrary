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
    
    def contact(self, i, flip, phi):
        #Determine indices of contacts and put together in one array
        argI = np.argwhere(self.J==i).flatten()
        argJ = np.argwhere(self.I==i).flatten()
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
            data[k,4] = i
            
            #Compute angle between particles
            theta = (np.arccos(self.nx[arg[k]])-phi)%(2*np.pi)
            
            #Compute components of forces
            fx=self.fn[arg[k]]*self.nx[arg[k]]+self.ft[arg[k]]*self.ny[arg[k]]
            fy=self.fn[arg[k]]*self.ny[arg[k]]-self.ft[arg[k]]*self.nx[arg[k]]
            
            #Make sure that contacts have equal but opposite forces
            if con[k] in J:
                fx=-fx
                fy=-fy
     
            #Add to array
            data[k,2] = fx
            data[k,3] = fy
            data[k,1] = theta

        #Sort array from smallest to largest angle
        data = data[data[:, 1].argsort()]
        
        #Flipping the data, don't think this is necessary
        if flip: data = np.flipud(data)
        
        return arg, con, data
    
    
    #Plotting the Maxwell Cremona tiling
    def tile(self, arrow):
        if isinstance(self.I, int):
            print('no data')
        else:
            #Initial values      
            xor1 = yor1 = phi1 = l = 0 #Staring position x, starting position y, starting angle, counter
            checklist0 = np.unique(np.union1d(self.I, self.J)) #List all original contacts
            checklist = np.unique(np.union1d(self.I, self.J)) #Checklist for checking if all contacts are plotted
            a = True #So we enter the first while loop
            i = checklist[0] #Start with first entry of the array, choosing different starting point should not change result
            priority = False #value used for searching connected component
            plot = True #Plot the tiles if true
            
            #Loop over all particles that have contacts
            while a:
                if len(checklist) == 0: break #if we have checked all contacts then exit
                if l >= 1: #If we are not in the first iteration check over all contacts
                    #Retrieve new contacts from the contactlist
                    for n in range(len(contactlist)):
                        i = contactlist[n]
                        if priority:
                            if i==inew: 
                                priority = False
                                plot = True
                                break #If we found the priority contact break out of the loop
                        
                        if priority:
                             if i in checklist0[checklist0!=i]: break #If we are in priority checking needs to be more lenient                      
                        
                        if not priority:
                            if i in checklist: break #If we found a new contact in the contactslist break out of the loop
                        
                        #If a new contact can not be found in the contactlist, get a new entry from the checklist 
                        if n == len(contactlist)-1: 
                            i=checklist[0]
                            temp, contacts, temp = self.contact(i, False, 0) #We are only interested in the contacts
                            #If particle is connected to current cluster we need to find correct origin coordinates
                            if contacts.any() in checklist0 and contacts.any() not in checklist: 
                                print('Warning -> adding connected tile of particle', i)
                                inew = i
                                priority = True
                                plot = False
                                l = 0
                                i = checklist0[0] #restart the sequence, but now end up in inew. This ensures correct xor1, yor1 and phi1
                            
                            #Particle is disconnected, so choose new coordinates 
                            else: 
                                print('Warning -> adding disconnected tile of particle', i)
                                xor1 = yor1 = l
                                phi1 = 0
                    
                    #Retrieve previous contact to get correct origin position
                    for n in range(len(con)):
                        arr = con[n]
                        if i in arr[1]:
                            contact = con[0]
                    
                    #Retrieve new angle from previous tiles
                    for n in range(len(data)):
                        arr = data[n]
                        x = np.argwhere(arr==i)
                        #If we have an entry equal to i and this is in contact with the previous entry
                        if len(x) != 0 and arr[x].flatten()[-1]==contact:
                            phi1 = arr[x].flatten()[1] #Retrieve angle
                            break
                    
                    #Retrieve starting position from previous tiles
                    for n in range(len(orr)):
                        arr = orr[n]
                        if i in arr:
                            x = np.argwhere(arr==i).flatten()[0]
                            n = (x+1)%len(arr)
                            xor1 = arr[n,1]
                            yor1 = arr[n,2]
                            break

                #Remove contact from the checklist
                checklist = checklist[checklist != i]
                
                #Counter and value for next while loop
                k = 0
                b = True
                
                #Get contact data and force data for contact i
                arg1, con1, data1 = self.contact(i, flip=False, phi=phi1)
                
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
                    #Generate colors for each tile
                    colors = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)]) for i in range(len(con1))]
                    
                    #Extract angle to rotate plotting
                    phi2 = data1[data1[:,0]==con1[k]].flatten()[1] + phi1 #add phi to make sure we plot in the correct order
                    
                    #Get contact data
                    arg2, con2, data2 = self.contact(con1[k], flip=False, phi=phi2)
                    
                    #Create list of all the contacts of the contacts
                    contactlist.extend(con2.tolist())
                    contactlist = [*set(contactlist)] #remove duplicates
                    
                    #Get origin coordinates
                    x = np.argwhere(orr1==con1[k]).flatten()[0]
                    n = (x+1)%len(con1)
                    xor2 = orr1[n,1]
                    yor2 = orr1[n,2]
                    #print(i, con1, con1[k], con2, xor2, yor2)
                    
                    #Plot the result and get origin coordinates
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
                plt.show()
            
    
    


                
                                
    