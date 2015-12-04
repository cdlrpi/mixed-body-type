
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
       

    # Mass-matricies constant (don't need to pass generalized coords - kinematics updated
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


nJpin = 3
nJfixed = (nbodies - nJpin) 

ndofsGEBF = 6   # (x,y,theta per joint)
ndofsRigid  = 1 # (one pin joint)

# GEBF-Properties
# Physical Properties
r   =  0.02393653682408596   #math.sqrt(A/math.pi)
A   =  2*r #0.0018
I   =  2/3*r**3 #1.215e-8
L   =  1.2

l  = 0.12
lr = 4*l
g  = 9.81

# Rigid-body Properties
# Physical Properties
mR = 1.0
IR = 1/12*mR*l**2

# Material Properties
E   =  0.7e5
rho =  5540

# Length of time of the simulation and time-step
t_final = 1.0
dt = 0.0001
tspan = np.arange(0,t_final,dt)


# In[4]:

math.sqrt(A/math.pi)


# #### Simulation Specifications

# In[5]:

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


# In[6]:

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
    body.initialize('rigid', mR,lr, IR, g)


# In[7]:

# odeint is the numerical integrator used
# state = MBF.myRK4(simulate,state0,tspan,nbodies,bodiesGEBF,bodiesRigid,joints,2,1)
# %%prun -s cumulative 
state = odeint(simulate,state0,tspan,(nbodies,bodiesGEBF,bodiesRigid,joints,2,1))


# In[8]:

q = state[:,:nRIGID*ndofsRigid + nGEBF*ndofsGEBF]
u = state[:,nRIGID*ndofsRigid + nGEBF*ndofsGEBF:]

# slice out rigid generalized coordinates and speeds from GEBF
qRigid = q[:,:nRIGID*ndofsRigid]
uRigid = u[:,:nRIGID*ndofsRigid]
qGEBF  = q[:,nRIGID*ndofsRigid:]
uGEBF  = u[:,nRIGID*ndofsRigid:]

r0 = np.array([[np.cos(q0GEBF[0][0])*(l*i+nRIGID*lr),   np.sin(q0GEBF[0][0])*(l*i+nRIGID*lr), 
                np.cos(q0GEBF[0][0])*(l*i+l+nRIGID*lr), np.sin(q0GEBF[0][0])*(l*i+l+nRIGID*lr)] 
               for i in range(0,nGEBF)],dtype=np.double).reshape(1,ndofsGEBF*nGEBF - 2*nGEBF)

xRigid,yRigid = MBF.get_topology_2DRigid(qRigid,lr)

# Interpolate between nodes with npoints, and compute absolute positions of nodes for plotting
npoints = 10
xGEBF,yGEBF = MBF.get_topology_2DGEBF(qGEBF, r0, l, nGEBF, npoints)

x  = [np.hstack((xR.reshape(1,len(xR)),xG)) for xR,xG in zip(xRigid,xGEBF)]
y = [np.hstack((yR.reshape(1,len(yR)),yG)) for yR,yG in zip(yRigid,yGEBF)]

# x = xGEBF
# y = yGEBF


# In[9]:

fig = plt.figure()
ax = fig.add_subplot(111, autoscale_on=False, xlim=(-0.5,2.5), 
                                              ylim=(-2.5,0.5))
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

# font = {'family' : 'serif',
#         'weight' : 'bold',
#         'size'   : 16}

# plt.rc('font', **font)


# In[11]:

ke_RIGID,pe_RIGID,te_RIGID = MBF.get_energy_Rigid(bodiesRigid,state)

# comptue energy from GEBF elements
se_GEBF = []
pe_GEBF = []
ke_GEBF = []
for q,u in zip(qGEBF,uGEBF):
    SE_tot_tstep = 0
    PE_tot_tstep = 0
    KE_tot_tstep = 0
    for body,qbody,ubody in zip(bodiesGEBF, np.array_split(q,nGEBF), np.array_split(u,nGEBF)):
        [SE, PE, KE] = body.energy('gebf', qbody, ubody)
        SE_tot_tstep = SE_tot_tstep + SE
        PE_tot_tstep = PE_tot_tstep + PE
        KE_tot_tstep = KE_tot_tstep + KE
    se_GEBF.append(SE_tot_tstep)
    pe_GEBF.append(PE_tot_tstep)
    ke_GEBF.append(KE_tot_tstep)
    
ke = ke_RIGID + np.array(ke_GEBF)
pe = pe_RIGID + np.array(pe_GEBF)
se = np.array(se_GEBF)

# ke = np.array(ke_GEBF)
# pe = np.array(pe_GEBF)
# se = np.array(se_GEBF)


# In[12]:

plt.rc('text', usetex=True)
plt.plot(np.array(tspan), ke, 'r', 
         np.array(tspan), pe, 'b',
         np.array(tspan), se, 'g',
         np.array(tspan), ke+pe+se, 'k')
plt.ylabel('Energy')
plt.xlabel('Time')
plt.show()


# In[13]:

# %qtconsole


# In[14]:

nsteps = tspan.size
energy = np.hstack((tspan.reshape(nsteps,1),
                  ke.reshape(nsteps,1),
                  pe.reshape(nsteps,1),
                  se.reshape(nsteps,1)))
np.savetxt('gebf-mixed-energy.mat', energy, fmt='%3.16f') 


# In[21]:

nsteps = tspan.size

for m in np.linspace(0,nsteps,5)
data = np.hstack((np.squeeze(x[m]).reshape(63,1), np.squeeze(y[m]).reshape(63,1)))
np.savetxt('gebf-mixed-dyn-%d.mat' %d, data, fmt='%3.16f') 
# data = np.hstack((np.squeeze(x[1000]).reshape(63,1), np.squeeze(y[1000]).reshape(63,1)))
# np.savetxt('gebf-mixed-dyn-1000.mat', data, fmt='%3.16f') 
# data = np.hstack((np.squeeze(x[1500]).reshape(63,1), np.squeeze(y[1500]).reshape(63,1)))
# np.savetxt('gebf-mixed-dyn-1500.mat', data, fmt='%3.16f') 
# data = np.hstack((np.squeeze(x[2000]).reshape(63,1), np.squeeze(y[2000]).reshape(63,1)))
# np.savetxt('gebf-mixed-dyn-2000.mat', data, fmt='%3.16f') 
# data = np.hstack((np.squeeze(x[2500]).reshape(63,1), np.squeeze(y[2500]).reshape(63,1)))
# np.savetxt('gebf-mixed-dyn-2500.mat', data, fmt='%3.16f') 
# data = np.hstack((np.squeeze(x[3000]).reshape(63,1), np.squeeze(y[3000]).reshape(63,1)))
# np.savetxt('gebf-mixed-dyn-3000.mat', data, fmt='%3.16f') 
# data = np.hstack((np.squeeze(x[3500]).reshape(63,1), np.squeeze(y[3500]).reshape(63,1)))
# np.savetxt('gebf-mixed-dyn-3500.mat', data, fmt='%3.16f') 
# data = np.hstack((np.squeeze(x[4000]).reshape(63,1), np.squeeze(y[4000]).reshape(63,1)))
# np.savetxt('gebf-mixed-dyn-3000.mat', data, fmt='%3.16f') 
# data = np.hstack((np.squeeze(x[4500]).reshape(63,1), np.squeeze(y[4500]).reshape(63,1)))
# np.savetxt('gebf-mixed-dyn-3000.mat', data, fmt='%3.16f') 
# data = np.hstack((np.squeeze(x[4999]).reshape(63,1), np.squeeze(y[4999]).reshape(63,1)))
# np.savetxt('gebf-mixed-dyn-3000.mat', data, fmt='%3.16f') 


# In[16]:

# plt.rc('text', usetex=True)
# fig = plt.figure()
# ax = fig.add_subplot(111, autoscale_on=False, xlim=(-3,3), 
#                                               ylim=(-3,3))
# plt.plot(np.squeeze(x[4999]), np.squeeze(y[3500]), 'o-', lw=2)
# # plt.ylabel('Energy')
# # plt.xlabel('Time')
# plt.show()


# In[17]:

#  data


# In[3]:

np.linspace(0,5000,5)


# In[ ]:



