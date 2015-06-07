import numpy as np
import math
from scipy.integrate import odeint
import matplotlib.pyplot as plt

import MBstructs as MB
import MultiBodyFuncts as MBF

def openR(n,i,bodies,joints,BC1,BC2): 
    """
    This function uses the DCA to form and solve the equations of motion.
    """
    
    #This if statements is the recursive case of this function.
    #Inside this if statement the current list of bodies is
    #assembled into a smaller list and passed back to this function.
    if len(bodies) !=1: #if there are still more than one bodies in the system
        j=0
        odd = 0
        
        #check if there are an odd number of bodies and
        #create a new list called "newbds" containing the correct number of
        #new bodies
        newbds=[]
        if (len(bodies)%2)==0:
            for k in range (0,math.trunc(len(bodies)/2)):
                newbds.append(MB.Body())
        else:
            odd=1
            for k in range (0,math.trunc(len(bodies)/2)+1):
                newbds.append(MB.Body())
            
        #Loop through all of the new bodies, assembling the correct two
        #"old" bodies to create the new.
        while j < len(newbds):
            
            #if there are an odd # of bodies and we're at the last 
            #body in newbds, add the last body to newbds
            if j == len(newbds)-1 and odd == 1:
                newbds[j].z11=np.zeros((6,6))
                newbds[j].z12=np.zeros((6,6))
                newbds[j].z13=np.zeros((6))
                newbds[j].z21=np.zeros((6,6))
                newbds[j].z22=np.zeros((6,6))
                newbds[j].z23=np.zeros((6))
                newbds[j].z11[:,:]=bodies[2*j].z11[:,:]
                newbds[j].z12[:,:]=bodies[2*j].z12[:,:]
                newbds[j].z21[:,:]=bodies[2*j].z21[:,:]
                newbds[j].z22[:,:]=bodies[2*j].z22[:,:]
                newbds[j].z13[:]=bodies[2*j].z13[:]
                newbds[j].z23[:]=bodies[2*j].z23[:]
                newbds[j].m=bodies[2*j].m
                
            #Otherwise, calculate the new zetas for the newly assembled bodies
            #according to the formulae
            else:
                newbds[j].m=bodies[2*j].m+bodies[2*j+1].m         
                #X, Q, and Y  are  intermediate quantities
                newbds[j].X=np.dot(np.transpose(joints[2*j+1].D),np.dot(bodies[2*j+1].z11+bodies[2*j].z22,joints[2*j+1].D))
                newbds[j].Xinv=np.linalg.inv(newbds[j].X)
                newbds[j].W=np.dot(joints[2*j+1].D,np.dot(newbds[j].Xinv,np.transpose(joints[2*j+1].D)))
                newbds[j].Y=np.dot(newbds[j].W,bodies[2*j].z23-bodies[2*j+1].z13)#ommitted pdot*u
                #assemble the bodies based on the formulas
                newbds[j].z11=bodies[2*j].z11-np.dot(bodies[2*j].z12,np.dot(newbds[j].W,bodies[2*j].z21))
                newbds[j].z12=np.dot(bodies[2*j].z12,np.dot(newbds[j].W,bodies[2*j+1].z12))
                newbds[j].z21=np.dot(bodies[2*j+1].z21,np.dot(newbds[j].W,bodies[2*j].z21))
                newbds[j].z22=bodies[2*j+1].z22-np.dot(bodies[2*j+1].z21,np.dot(newbds[j].W,bodies[2*j+1].z12))
                newbds[j].z13=bodies[2*j].z13-np.dot(bodies[2*j].z12,newbds[j].Y)
                newbds[j].z23=bodies[2*j+1].z23+np.dot(bodies[2*j+1].z21,newbds[j].Y)
            j=j+1
            
        
        #Find the new joints corresponding to the new bodies
        newjs=[]
        for k in range (0,len(newbds)):
            # newjs.append(MB.Joint())
            if j==len(newjs)-1 and odd ==1:
                newjs.append(MB.Joint(joints[len(joints)-1].P,np.max(joints[len(joints)-1].P)))
        
            else:
                newjs.append(MB.Joint(joints[2*k].P,np.max(joints[2*k].P)))
        
        # #loop that brings along the P and D matrices from the previous joints
        # for j in range(0,len(newjs)):
        #     if j==len(newjs)-1 and odd ==1:
        #         newjs[j].D=joints[len(joints)-1].D
        #         newjs[j].P=joints[len(joints)-1].P
        #        
        #     else:
        #         newjs[j].D=joints[2*j].D
        #         newjs[j].P=joints[2*j].P
               
                
                
        #This is the recursive call.  This will return a list of the form:
        #[A11,F1c1,A12,F1c2,A21,F2c1,...,An2,Fnc2] where the Axy's and Fxcy's 
        #correspond to the accellerations and constraint forces in the 
        #xth body on the yth handle. The function is given the original 
        #number of bodies, the list of new bodies, the list of new joints,
        #the boundary conditions,and the state variables at that timestep.
        #At this point this function will repeat itself in a loop until there is
        #one body left

        sol=openR(n,i+1,newbds,newjs,BC1,BC2)
       
        
        #Forces and Accelerations at the new joints are found,
        #so these values can now be found for the old joints.
        
        #newsol will contain the new solution
        newsol=[];
        newsol.append(sol[0])
        newsol.append(sol[1])
        
       
        
        #In this loop I start with the first body and find the force and 
        #acceleration on its other handle. These values are added to the new
        #solution list, along with the next values in sol.
        flag=0
        for j in range(0,len(newbds)):
            
            #Don't enter this if there are an odd number of bodies and it 
            #is the last time through the loop, otherwise enter.
            if not(odd==1 and j==len(newbds)-1):
                
                #index counter
                k=len(newsol)-2
                
                #The force in the joint between two assembled bodies
                #F=np.dot(np.linalg.inv(bodies[2*j].z12),newsol[k]-bodies[2*j].z13-np.dot(bodies[2*j].z11,newsol[k+1]))
                F=-1*(np.dot(joints[j+1].D,np.dot(newbds[j].Xinv,np.dot(np.transpose(joints[j+1].D),(np.dot(bodies[2*j].z21,newsol[k+1])-np.dot(bodies[2*j+1].z12,sol[4*j+3])+bodies[2*j].z23-bodies[2*j+1].z13)))))
                
                #A2 is the acceleration of handle 2 of body k
                A2=np.dot(bodies[2*j].z21,newsol[k+1])+np.dot(bodies[2*j].z22,F)+bodies[2*j].z23
               
                #A1 is the acceleration of handle 1 of body k+1
                A1=np.dot(bodies[2*j+1].z11,-1*F)+np.dot(bodies[2*j+1].z12,sol[4*j+3])+bodies[2*j+1].z13
                 
                #Add the newly found values to new solution list,
                #maintiaining the format described above.
                newsol.append(A2)
                newsol.append(F)
                newsol.append(A1)
                newsol.append(-1*F)
                
                #Add the next two values in sol to the new solution list
                newsol.append(sol[4*j+2])
                newsol.append(sol[4*j+3])
            
            #if there are an odd number of bodies and this is the last
            #time through the loop, append the last force and accelleration 
            #in the old solution to the new one
            else:
                newsol.append(sol[4*j+2])
                newsol.append(sol[4*j+3])
                flag=1
            
            #If this is not the last time through the loop, append the 
            #next two values in the old solution to the new solution
            if j!=len(newbds)-1  and flag !=1:   
                newsol.append(sol[4*j+4])
                newsol.append(sol[4*j+5])
        
        #If this is the 0th level of recursion, delete all the forces 
        #so only the accellerations are returned
        if i ==0:
            for j in range (1,2*n+1):
                del newsol[j]
                  
            return newsol
        
        #If this is not the 0th level of recursion, return the solution
        #with all of the forces included
        else:
            return newsol
    
    #This is the base case of this recursive function.  This is where
    #the recursion ends and reassembly begins
    elif BC1==2 and BC2==1:
        
        #Fc2 is zero because that end is free
        Fc2=np.zeros((6))
        
        #Because the first joint is a pin, it will
        #have no translational motion, and support no 
        #moments, so mat is used to solve for the two translational
        #forces of Fc1
        mat=np.zeros((3,4))
        mat[:,:3]=bodies[0].z11[3:,3:]
        mat[:,3]=-1*bodies[0].z13[3:]
        
        #Loop to put a matrix in reduced row echelon form
        for s in range(0,3):
            mat[s,:]=mat[s,:]/mat[s,s]
            for j in range(0,3):
                if s !=j:
                    mat[j,:]=mat[j,:]-(mat[s,:]*mat[j,s])
    
        Fc1=np.zeros_like(Fc2)
        Fc1[3:]=mat[:,3]

        #solve for the A's given Fc2=0
        A1=np.dot(bodies[0].z11,Fc1)+bodies[0].z13
        A2=np.dot(bodies[0].z21,Fc1)+bodies[0].z23
        sol=[A1,Fc1,A2,Fc2]
        
        #return the forces and accellerations at the end joints
        #and begin disassembly
        return sol  
