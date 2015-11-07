
# coding: utf-8

# ## This notebook simulates a pendulum created of some rigid-bodies and some highly flexible GEBF bodies

# In[1]:

import pdb as pdb 
import math
import time
import numpy as np
import scipy as sc
import cProfile
import pdb as pdb 

from scipy.integrate import odeint

import matplotlib.pyplot as plt
import matplotlib.animation as animation

import MBstructs as MBS
import MBfuncts as MBF
import DCA

from IPython.display import display
from __future__ import division
from sympy.interactive import printing
printing.init_printing(use_latex='mathjax')
np.set_printoptions(precision=4,suppress=True)


# In[2]:

def simulate(state,tspan,nbodies,bodiesGEBF,bodiesRigid,joints,BC1,BC2):
    '''
    This function extracts the generalized coordinates from the solution 
    of the equations of motion after calling the DCA to solve the equations of motion 
    '''
    # only pin joints connecting rigid bodies
    q = state[:nRIGID*ndofsRigid + nGEBF*ndofsGEBF]
    u = state[nRIGID*ndofsRigid + nGEBF*ndofsGEBF:]
    
    # slice out rigid generalized coordinates and speeds from GEBF
    qRigid = q[:nRIGID*ndofsRigid]
    uRigid = u[:nRIGID*ndofsRigid]
    qGEBF  = q[nRIGID*ndofsRigid:]
    uGEBF  = u[nRIGID*ndofsRigid:]
    
    # Update the Kinematics 2D n-link Pendulum
    # Now that these functions "woe
    # This should be a class function for each body type not an absolute function
    # e.g., [body.updateKin for body in bodies]
    MBF.kinematics_Rigid2D(bodiesRigid,qRigid,uRigid)  
#     MBF.kinematics_GEBF2D(bodiesGEBF,qGEBF,uGEBF)
    
#     # slice out generalized cooridnates
#     # and make a sublist for each body
#     # these u's are not generalized speeds these are the displacements of the nodes
#     u = qGEBF.tolist()
#     u = [u[i*ndofsGEBF:(i*ndofsGEBF)+ndofsGEBF] for i in range((int(len(u)/ndofsGEBF)))]

#     # 'pop' off the rotational coordinates 
#     # those are dealt with through the 'kinematic' sweep
#     for ue in u:
#         ue.pop(3)
#         ue.pop(0)
        
    # compute the inverse inertial properties of the body 
    # with the updated generalized speeds
    
    for body,q,u in zip(bodiesGEBF, np.array_split(qGEBF,nGEBF), np.array_split(uGEBF,nGEBF)):
        body.intProps('gebf', q, u)

    # Mass-matricies constant (don't need to pass generalized coords)
    for body in bodiesRigid:
        body.intProps('rigid')

    # join the lists of bodies for a total list 
    # order matters
    bodies = bodiesRigid + bodiesGEBF
    
    # Call the Recursive DCA Algorithm
    # This returns a list of the form:
    # [A11,A12,A21,A22,...,An1,An2]
    # where Axy corresponds to the acceleration
    # of the yth handle of the xth body
    accel = DCA.solve(nbodies,0,bodies,joints,BC1,BC2)
    
    accelRigid = accel[:2*len(bodiesRigid)]
    accelGEBF  = accel[2*len(bodiesRigid):]

    # compute the generalized accelerations for kinematic joints
    udot_Rigid = MBF.get_gen_accel_Rigid(len(bodiesRigid), joints, accelRigid)
    udot_GEBF  = MBF.get_gen_accel_GEBF(len(bodiesGEBF), accelGEBF)
    udot = np.hstack((udot_Rigid,udot_GEBF))

    state_dot = np.hstack((state[nRIGID*ndofsRigid + nGEBF*ndofsGEBF:],udot))
    return state_dot


# #### System Initilization

# In[3]:

# number of bodies and GEBF elements and kinematic joints
# *** NOTE: For testing pupposes DCA does NOT effectively work for one body
#          - For one body DCA returns [A1 F1c A2 F2c]
#          - For two or more bodies DCA returns [A1 A2] ***
nRIGID = 2
nGEBF  = 6
nbodies = nGEBF + nRIGID


nJpin = 1
nJfixed = (nbodies - nJpin) 

ndofsGEBF = 6   # (x,y,theta per joint)
ndofsRigid  = 1 # (one pin joint)

# GEBF-Properties
# Physical Properties
A   =  0.0018
I   =  1.215e-8
L   =  1.2
r   =  math.sqrt(A/math.pi)
l = 0.12
g = 9.81

# Rigid-body Properties
# Physical Properties
mR = 1.0
IR = 1/12*mR*l**2

# Material Properties
E   =  0.7e6
rho =  5540

# Length of time of the simulation and time-step
t_final = 0.002
dt = 0.0001
tspan = np.arange(0,t_final,dt)


# #### Simulation Specifications

# In[4]:

# Initial conditions

# sublists handy for initilization of bodies 
q0GEBF = [[0 for i in range(ndofsGEBF)] for i in range(nGEBF)]
# q0GEBF[0][0] = -np.pi/2
u0GEBF = [[0 for i in range(ndofsGEBF)] for i in range(nGEBF)]
# u0GEBF[0][0] = -1

# All rigid bodies are connected by pin joints
q0Rigid = [[0 for i in range(ndofsRigid)] for i in range(nRIGID)]
u0Rigid = [[0 for i in range(ndofsRigid)] for i in range(nRIGID)]
# u0Rigid[0][0] = -1 

q0 = q0Rigid + q0GEBF
u0 = u0Rigid + u0GEBF

# flatten sub-lists into a single state vector
q0 = np.array([q for q_body in q0 for q in q_body],dtype=np.double).squeeze()
u0 = np.array([u for u_body in u0 for u in u_body],dtype=np.double).squeeze()

state0 = np.hstack((q0,u0))


# In[5]:

# Initialize a list of Joints
jointsPin = [MBS.Joint('revolute2D') for i in range(nJpin)]
jointsFixed = [MBS.Joint('fixed') for i in range (nJfixed)]
joints = jointsPin + jointsFixed

# Initialize the list of bodies 
bodiesGEBF = [MBS.Body() for i in range(nGEBF)]
bodiesRigid = [MBS.Body() for i in range(nRIGID)]
for body in bodiesGEBF:
    body.initialize('gebf', E, A, I, r, rho, l, g, q0GEBF)
for body in bodiesRigid:
    body.initialize('rigid', mR,l, IR, g)


# In[6]:

# odeint is the numerical integrator used
# state = MBF.myRK4(simulate,state0,tspan,nbodies,bodiesGEBF,bodiesRigid,joints,2,1)
# %%prun -s cumulative 
state = odeint(simulate,state0,tspan,(nbodies,bodiesGEBF,bodiesRigid,joints,2,1))


# In[7]:

q = state[:,:nRIGID*ndofsRigid + nGEBF*ndofsGEBF]
u = state[:,nRIGID*ndofsRigid + nGEBF*ndofsGEBF:]

# slice out rigid generalized coordinates and speeds from GEBF
qRigid = q[:,:nRIGID*ndofsRigid]
uRigid = u[:,:nRIGID*ndofsRigid]
qGEBF  = q[:,nRIGID*ndofsRigid:]
uGEBF  = u[:,nRIGID*ndofsRigid:]

r0 = np.array([[np.cos(q0GEBF[0][0])*(l*i+nRIGID*l),   np.sin(q0GEBF[0][0])*(l*i+nRIGID*l), 
                np.cos(q0GEBF[0][0])*(l*i+l+nRIGID*l), np.sin(q0GEBF[0][0])*(l*i+l+nRIGID*l)] 
               for i in range(0,nGEBF)],dtype=np.double).reshape(1,ndofsGEBF*nGEBF - 2*nGEBF)

xRigid,yRigid = MBF.get_topology_2DRigid(qRigid,0.12)
# Interpolate between nodes with npoints, and compute absolute positions of nodes for plotting
npoints = 10
xGEBF,yGEBF = MBF.get_topology_2DGEBF(qGEBF, r0, l, nGEBF, npoints)

x  = [np.hstack((xR.reshape(1,len(xR)),xG)) for xR,xG in zip(xRigid,xGEBF)]
y = [np.hstack((yR.reshape(1,len(yR)),yG)) for yR,yG in zip(yRigid,yGEBF)]

# x = xGEBF
# y = yGEBF


# In[8]:

# ntsteps, len_q = np.shape(q)
    
# displacements = np.delete(q,range(0,len_q,3),1)

# # compute position of each node
# position = r0 + displacements
# position_element = np.array_split(position,nelements,axis=1)

# # interpolate 
# x = np.linspace(-1,1,npoints)
# h1 = (1 - x/l)
# h2 = (x/l)

# # Compute shape function matrix
# H = np.array([np.hstack((h1*np.eye(2), h2*np.eye(2))) 
#               for h1,h2 in zip(h1,h2)]).reshape(2*npoints,4).T

# # Interpolate the X and Y cooridnates
# xy = np.hstack(np.array([np.dot(position_element,H) 
#                          for position_element in position_element]))

# x = np.array_split(xy[:,0::2],ntsteps+1)
# y = np.array_split(xy[:,1::2],ntsteps+1)

# return x,y


# In[9]:

fig = plt.figure()
ax = fig.add_subplot(111, autoscale_on=False, xlim=(-(nbodies+5)*l,(nbodies+5)*l), 
                                              ylim=(-(nbodies+5)*l,(nbodies+5)*l))
# ax = fig.add_subplot(111, autoscale_on=False, xlim=(-10,10), 
#                                               ylim=(-10,10))

ax.grid()

line, = ax.plot([], [], 'o-', lw=2)
time_template = 'time = %0.4fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

def init():
    line.set_data([], [])
    time_text.set_text('')
    return line, time_text

def animate(i):
#     thisx = [0, x[i][1], x[i][2]]
#     thisy = [0, y[i][1], y[i][2]]
    thisx = x[i]
    thisy = y[i]

    line.set_data(thisx, thisy)
    time_text.set_text(time_template%(i*dt))
    return line, time_text

ani = animation.FuncAnimation(fig, animate, np.arange(1, len(y)),
    interval = 1, blit=False, init_func=init)

# ani.save('rigid_rope.mp4', fps=60)
plt.show()


# In[10]:

# ke,pe,te = MBF.get_energy_Rigid(bodiesRigid,state)


# In[11]:

# # The energy of the system is calculated and plotted
# %matplotlib inline 
# plt.plot(tspan,te-te[0])
# plt.xlabel("Time [s]")
# plt.ylabel("energy")
# plt.title("System Energy")
# plt.show()

# plt.plot(tspan,pe,tspan,ke)
# plt.xlabel("Time[s]")
# plt.ylabel("energy")
# plt.title("Kinetic and Potential Energy")
# plt.show


# In[12]:

# %qtconsole


# In[ ]:



