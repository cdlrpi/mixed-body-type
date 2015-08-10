# Jeremy Laflin
# Modified version of MultiBodyFuncts.py by Mike Hans
# In this version
#   - Rotation matricies are reduced to 2x2 simple rotation about 3rd axis
import numpy as np
#This script will hold some functions that will be used in other scripts
def DCM(theta):
    C = np.array([[np.cos(theta), -np.sin(theta)],\
                  [np.sin(theta),  np.cos(theta)]])
    return C

# Once these work they should be moved to class functions (maybe)
# Plotting Functions
def get_topology_2DRigid(q,bodies):
    """
    This fuction determined the location of the endpoints of a rigid pendlum.  
    """

    theta = np.cumsum(q,1)
    C = [[np.array([[np.cos(angle),-np.sin(angle)],\
                    [np.sin(angle), np.cos(angle)]]) 
      for angle in tstep] for tstep in theta]
    R = [[np.dot(CBN,np.array([body.l,0])) for CBN in timestep] \
                                      for (timestep,body) in zip(C,bodies)]
    R = [np.cumsum(r,0) for r in R]

    x = [np.hstack((np.array([0]),r[:,0])) for r in R]
    y = [np.hstack((np.array([0]),r[:,1])) for r in R]

    return x,y
    
def get_topology_2DGEBF(q,r0,nelements,npoints):
    """
    This function determines the positions of all the nodes of the gebf elements of the system. 
    """
    
    # delete rotational degrees of freedom
    # remaining dof's are theta1,x1,y1,theta2,x2,y2 --> x1,y1,x2,y2 
    ntsteps, len_q = np.shape(q)
    
    displacements = np.delete(q,range(0,len_q,3),1)
    
    # compute position of each node
    position = np.vstack((r0,displacements))
    position = np.cumsum(position,0)
    
    position_element = np.array_split(position,nelements,axis=1)
    
    
    # interpolate 
    x = np.linspace(-1,1,npoints)
    h1 = (1/2)*(1 - x)
    h2 = (1/2)*(1 + x)

    # Compute shape function matrix
    H = np.array([np.hstack((h1*np.eye(2), h2*np.eye(2))) 
                  for h1,h2 in zip(h1,h2)]).reshape(2*npoints,4).T
    
    # Interpolate the X and Y cooridnates
    xy = np.hstack(np.array([np.dot(position_element,H) 
                             for position_element in position_element]))
    
    x = np.array_split(xy[:,0::2],ntsteps)
    y = np.array_split(xy[:,1::2],ntsteps)
    return x,y
    

# Energy Functions
def get_energy_Rigid(bodies,state):
    """
    This function determines the kinetic energy of the rigid-bodies of a pendulum given the generalized coordinates of the system.
    """
    
    ntsteps,nbodies = np.shape(state)
    q = state[:,:nbodies/2]
    u = state[:,nbodies/2:]

    theta = np.cumsum(q,1)
    omega = np.cumsum(u,1)

    C = [[np.array([[np.cos(angle),-np.sin(angle)],\
                    [np.sin(angle), np.cos(angle)]]) for angle in tstep] 
                                                     for tstep in theta]

    ## Potential Energy Calculation ##
    # vector location of joints for all bodies
    Rh2Body = [[np.dot(CBN,np.array([body.l,0])) 
                                      for (CBN,body) in zip(timestep,bodies)] 
                                      for timestep in C]
    Rh2 = [np.cumsum(r,0) for r in Rh2Body]

    # vector location of center of mass for all bodies
    # from handle two
    rcmBody = [[np.dot(CBN,np.array([-body.l/2,0])) 
                                        for (CBN,body) in zip(timestep,bodies)]
                                        for timestep in C]

    # locate the centers of mass w.r.t. origin
    rcm = [np.array(Rbody + rbody) for Rbody,rbody in zip(Rh2,rcmBody)]

    # slice out c.m. y position
    ycm = [y[:,1] for y in rcm]

    # compute heights of centers of mass
    hcmRel = [body.l/2 for body in bodies]
    lcm = np.cumsum([body.l for body in bodies])
    lcm = lcm.tolist()
    lcm.pop(len(bodies)-1)
    lcm.insert(0,0)
    lcm = np.array(lcm)
    hcm = [lcm+hcmRel+y for y in ycm]

    pe = np.sum([[9.81*body.m*h 
                for body, h in zip(bodies,timestep)]
                for timestep in hcm],1)

    ## Kinetic Energy Calculation ##
    vh2Body = np.cumsum([[qdot*np.array([-R[1],R[0]]) 
                for qdot,R in zip(timestep,R)] 
                for timestep,R in zip(omega,Rh2Body)],1)

    vcmBody = [[qdot*np.array([-r[1],r[0]]) 
                for qdot,r in zip(tstep,r_tstep)] 
                for tstep,r_tstep in zip(omega,rcmBody)]

    vcm = vcmBody+vh2Body

    keT = np.sum([[1/2*body.m*np.dot(v,v) 
                   for body,v in zip(bodies,timestep)] 
                   for timestep in vcm],1)

    keR = np.sum([[1/2*body.I*qdot**2 
                for body,qdot in zip(bodies,timestep)]
                for timestep in omega],1)
    ke = keT + keR
    te = ke + pe
    return ke,pe,te

# Functions to determine generalized accelerations from the spatial accelerations
# DCA returns Spatial Accelerations
# Therefore, need to extract generalized accelerations
def get_gen_accel_Rigid(nbodies, joints, accel):
    """
    This function determines the generalized acclerations between rigid-bodies 
    """
    udot = np.zeros((nbodies))
    for j in range(nbodies):
        if j == 0:
            A1 = accel.pop(0)
            udot[j] = np.dot(np.transpose(joints[j].P),A1)
        else:
            A2 = accel.pop(0)
            A1 = accel.pop(0)
            udot[j] = np.dot(np.transpose(joints[j].P),(A1-A2))
    
    #add the velocities to d_dt and return to the integrator
    return udot 

def get_gen_accel_GEBF(nGEBF, accel):
    """
    This function determines the generalized accelerations between GEBF elements
    """
    ndofsGEBF = 6
    state_dot = np.array(accel,dtype=np.double).reshape(nGEBF*ndofsGEBF)
    return state_dot 

def kinematics_Rigid2D(bodies,q,u):
    """
    This function determines the kinematics of a 2D rigid pendulum
    """

    for body,theta,omega in zip(bodies,np.cumsum(q),np.cumsum(u)):
        body.omega = omega
        body.CBN = DCM(theta)
            
def kinematics_GEBF2D(bodies,q,u):
    """
    This function determines the kinematics of a 2D gebf pendulum
    """
    # slice state into 'qs' for each element 
    thetae = q[::3]
    theta = np.cumsum(thetae)

    # slice into the two rotational coordinates for each body
    theta1 = theta[0::2]
    theta2 = theta[1::2]
    for body,q0,q3 in zip(bodies, theta1, theta2):
        body.theta1 = q0
        body.theta2 = q3

