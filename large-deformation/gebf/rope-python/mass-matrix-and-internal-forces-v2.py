
# coding: utf-8

# #Beam Element - This notebook computes the mass matrix and internal forces and outputs them for later use

# In[1]:

import numpy as np
import scipy as sp
import sympy as sym
import pickle

from sympy import cos, sin

from IPython.display import display
from __future__ import division
from sympy.interactive import printing
printing.init_printing(use_latex='mathjax')


# In[2]:

sym.mpmath.mp.pretty=True


# In[3]:

def Skew_sym(v):
    """
    This function returns the skew symetric matrix 
    of the vector 'v' to affect the cross product of 'v'x'u'
    """
    v_skew = sym.Matrix([[  0 , -v[2],  v[1]],
                         [v[2],     0, -v[0]],
                         [-v[1],  v[0],    0]])
    return v_skew


# In[4]:

def Rotate_sym(theta):
    """
    This function returns the symbolic rotation matrix 
    for the simple 2D-rotation about the third axis
    """
    R = sym.Matrix([[cos(theta), sin(theta), 0],
                   [-sin(theta), cos(theta), 0],[0,0,1]])
    return R


# #### Define symbolic quantites

# In[5]:

# nodal positons, velocities, and accelerations
"""
"""

# symbolic system parameters 
E, I, A, rho, x, l, tau_a = sym.symbols('E I A rho x l, tau_a')

# Beam coordinates at cross-section 1 and 2
r1 = sym.Matrix(['r_1x','r_1y'])
r2 = sym.Matrix(['r_2x','r_2y'])

# Deformation coordinates
u = sym.Matrix(['u_1x','u_1y','u_1z','u_2x','u_2y','u_2z'])
# u1 = sym.Matrix(['u_1x','u_1y','u_1z'])
# u2 = sym.Matrix(['u_2x','u_2y','u_2z'])
udot = sym.Matrix(['udot_1x','udot_1y','udot_1z','udot_2x','udot_2y','udot_2z'])
# u1dot = sym.Matrix(['udot_1x','udot_1y','udot_1z'])
# u2dot = sym.Matrix(['udot_2x','udot_2y','udot_2z'])
uddot = sym.Matrix(['uddot_1x','uddot_1y','uddot_1z','uddot_2x','uddot_2y','uddot_2z'])
# u1ddot = sym.Matrix(['uddot_1x','uddot_1y','uddot_1z'])
# u2ddot = sym.Matrix(['uddot_2x','uddot_2y','uddot_2z'])

theta1 = sym.Matrix(sym.symarray('theta',4))
theta2 = sym.Matrix(sym.symarray('phi',4))
theta1_dot = sym.Matrix(sym.symarray('thetadot',4))
theta2_dot = sym.Matrix(sym.symarray('phidot',4))
theta1_ddot = sym.Matrix(sym.symarray('thetaddot',4))
theta2_ddot = sym.Matrix(sym.symarray('phiddot',4))

# locating the point in the cross section
s1 = sym.Matrix([0, r1])
s2 = sym.Matrix([0, r2])
s = sym.Matrix.vstack(s1,s2)

# define state variables
# e = sym.Matrix([u1, [theta1[3]], u2, [theta2[3]]])
e = sym.Matrix([u[0:3,0], [theta1[3]], u[3:6,0], [theta2[3]]])
# define state variables
edot = sym.Matrix([udot[0:3,0], [theta1_dot[3]], udot[3:6,0], [theta2_dot[3]]])
# edot = sym.Matrix([u1dot, [theta1_dot[3]], u2dot, [theta2_dot[3]]])
eddot = sym.Matrix([uddot[0:3,0], [theta1_ddot[3]], uddot[3:6,0], [theta2_ddot[3]]])
# eddot = sym.Matrix([u1ddot, [theta1_ddot[3]], u2ddot, [theta2_ddot[3]]])
display(e)


# ### Needed Matrix Quantities

# In[6]:

# angular velocity for 2D case is trivial
omega1_skew = Skew_sym(theta1_dot[1:4])
omega2_skew = Skew_sym(theta2_dot[1:4])
alpha1_skew = Skew_sym(theta1_ddot[1:4])
alpha2_skew = Skew_sym(theta2_ddot[1:4])

R1 = Rotate_sym(theta1[3])
R2 = Rotate_sym(theta2[3])

R = sym.Matrix.vstack(sym.Matrix.hstack(R1,sym.zeros(3)),                       sym.Matrix.hstack(sym.zeros(3),R2))

Omega_skew = sym.Matrix.vstack(sym.Matrix.hstack(omega1_skew,sym.zeros(3)),                                sym.Matrix.hstack(sym.zeros(3),omega2_skew))

Alpha_skew = sym.Matrix.vstack(sym.Matrix.hstack(alpha1_skew,sym.zeros(3)),                                sym.Matrix.hstack(sym.zeros(3),alpha2_skew))
# display(Omega)


# ### Define Kinematics

# In[7]:

# Define velocity of element endpoints (nodes)
v = udot + R*Omega_skew*s
# v1 = u1dot + R1*omega1_skew*s1
# v2 = u2dot + R2*omega2_skew*s2
# v = sym.Matrix([v1,v2])
print('v = ')
display(v)

# Define acceleration of element endpoints (nodes)
a = uddot + R*Omega_skew*Omega_skew*s + R*Alpha_skew*s
# a1 = u1ddot + R1*omega1_skew*omega1_skew*s1 + R1*alpha1_skew*s1
# a2 = u2ddot + R2*omega2_skew*omega2_skew*s2 + R2*alpha2_skew*s2
# a = sym.Matrix([a1,a2])
print('\na = ')
display(a)


# ### Compute the Mass Matrix

# In[21]:

# Define shape function for element with one node at each end
h = sym.symarray('h', 2)

h[0] = 1/2*(1 - x)
h[1] = 1/2*(1 + x)

# Compute shape function matrix
H = sym.expand(sym.Matrix([h[0]*sym.eye(3), h[1]*sym.eye(3)])).T
# print('\nH = ')
# display(H)

# Define velocity of any point 
Vp = H*v
# print('\nV = ')
# display(V)

# Define velocity of any point 
Ap = H*a
# print('\nA = ')
# display(Accel)

# Compute partial velocities of the nodes
Vr = sym.Matrix([[sym.diff(Vp[j],edot[i]) 
                         for j in range(len(V))] for i in range(len(edot))]).T

# v_r = H
print('\nVr = ')
display(Vr)

# Compute mass matrix
M = sym.Matrix([[sym.expand(sym.integrate(Vr[:,i].dot(Ap)*rho*A,(x,0,l))).coeff(eddot[j]) for i in range(len(eddot))]
                for j in range(len(eddot))])
print('\nM = \n')
display(M)
pickle.dump( M, open( "mass-matrix.dump", "wb" ) )


# ### Compute Internal forces 

# #### Compute Longitudinal Strain Energy

# In[9]:




# #### Compute Transverse Strain Energy

# In[9]:




# #### Compute Internal Forces $Q_e = \frac{\partial U}{\partial e}$

# In[9]:




# ### Applied and body force vector

# In[9]:



