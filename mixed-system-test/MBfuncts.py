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
def get_topology_2DRigid(q,l):
    """
    This fuction determined the location of the endpoints of a rigid pendlum.  
    """

    theta = np.cumsum(q,1)
    C = [[np.array([[np.cos(angle),-np.sin(angle)],\
                    [np.sin(angle), np.cos(angle)]]) 
      for angle in tstep] for tstep in theta]
    R = [[np.dot(CBN,np.array([l,0])) for CBN in timestep] \
                                      for timestep in C]
    R = [np.cumsum(r,0) for r in R]

    x = [np.hstack((np.array([0]),r[:,0])) for r in R]
    y = [np.hstack((np.array([0]),r[:,1])) for r in R]

    return x,y
    
def get_topology_2DGEBF(q,r0,l,nelements,npoints):
    """
    This function determines the positions of all the nodes of the gebf elements of the system. 
    """
    
    # delete rotational degrees of freedom
    # remaining dof's are theta1,x1,y1,theta2,x2,y2 --> x1,y1,x2,y2 
    ntsteps, len_q = np.shape(q)
    
    displacements = np.delete(q,range(0,len_q,3),1)
    
    # compute position of each node
    position = r0 + displacements
    position_element = np.array_split(position,nelements,axis=1)
    
    # interpolate 
    x = np.linspace(0,l,npoints)
    h1 = (1 - x/l) 
    h2 = (x/l)

    # Compute shape function matrix
    H = np.array([np.hstack((h1*np.eye(2), h2*np.eye(2))) 
                  for h1,h2 in zip(h1,h2)]).reshape(2*npoints,4).T
    
    # Interpolate the X and Y cooridnates
    xy = np.hstack(np.array([np.dot(position_element,H) 
                             for position_element in position_element]))
    
    x = np.array_split(xy[:,0::2],ntsteps+1)
    y = np.array_split(xy[:,1::2],ntsteps+1)

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
    hcm = ycm
    # hcm = [np.array(lcm)+np.array(hcmRel)]

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
            udot[j] = np.dot(joints[j].P.T,A1)
        else:
            A2 = accel.pop(0)
            A1 = accel.pop(0)
            udot[j] = np.dot(joints[j].P.T,(A1-A2))
    
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
    Needs fixing !!!!!!!!!!!!!!!!!!!!!!
    """
    #  # slice state into 'qs' for each element 
    #  thetae = q[::3]
    #  theta = np.cumsum(thetae)

    #  # slice into the two rotational coordinates for each body
    #  theta1 = theta[0::2]
    #  theta2 = theta[1::2]
    #  for body,q0,q3 in zip(bodies, theta1, theta2):
    #      body.theta1 = q0
    #      body.theta2 = q3
    raise Exception('Fix me!!!!')

def myRK4(func, state0, tspan, *args):
    """
    Inputs:
        func - function that returns the derivatives of the state variables
          x0 - starting conditions
           t - array of time values 
   Outputs:
       w - approximate solution at the mesh-points (time-steps)
    """
    dt = np.max(np.diff(tspan))
    nsteps = len(tspan)
    state = np.zeros((nsteps-1,len(state0)))
    state = np.vstack((state0,state))
        
    for i in range(len(tspan)-1):
        k1 = dt*func(state[i,:],        tspan,*args)
        k2 = dt*func(state[i,:] + k1/2, tspan,*args)
        k3 = dt*func(state[i,:] + k2/2, tspan,*args)
        k4 = dt*func(state[i,:] + k3,   tspan,*args)
        
        state[i+1,:] = state[i,:] + (1/6*((k1)+(2*k2)+(2*k3)+(k4))).reshape(1,len(state0))
    return state
