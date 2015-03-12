
# coding: utf-8

# #Beam Element - This notebook computes the mass matrix and internal forces and outputs them for later use

# In[ ]:

import numpy as np
import scipy as sp
import sympy as sym
import pickle

from sympy import cos, sin

from IPython.display import display
from __future__ import division
from sympy.interactive import printing
printing.init_printing(use_latex='mathjax')


# In[ ]:

sym.mpmath.mp.pretty=True


# In[ ]:

def Skew_sym(v):
    """
    This function returns the skew symetric matrix 
    of the vector 'v' to affect the cross product of 'v'x'u'
    """
    v_skew = sym.Matrix([[  0 , -v[2],  v[1]],
                         [v[2],     0, -v[0]],
                         [-v[1],  v[0],    0]])
    return v_skew


# In[ ]:

def Rotate_sym(theta):
    """
    This function returns the symbolic rotation matrix 
    for the simple 2D-rotation about the third axis
    """
    R = sym.Matrix([[cos(theta), sin(theta), 0],
                   [-sin(theta), cos(theta), 0],[0,0,1]])
    return R


# #### Define symbolic quantites

# In[ ]:

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
udot = sym.Matrix(['udot_1x','udot_1y','udot_1z','udot_2x','udot_2y','udot_2z'])
uddot = sym.Matrix(['uddot_1x','uddot_1y','uddot_1z','uddot_2x','uddot_2y','uddot_2z'])

theta = sym.Matrix(['theta_1x','theta_1y','theta_1z','theta_2x','theta_2y','theta_2z'])
omega = sym.Matrix(['omega_1x','omega_1y','omega_1z','omega_2x','omega_2y','omega_2z'])
alpha = sym.Matrix(['alpha_1x','alpha_1y','alpha_1z','alpha_2x','alpha_2y','alpha_2z'])

# locating the point in the cross section
s1 = sym.Matrix([0, r1])
s2 = sym.Matrix([0, r2])
s = sym.Matrix.vstack(s1,s2)

# define state variables
e = sym.Matrix([u[3:6,0], theta[0:3,0], u[3:6,0], theta[3:6,0]])
edot = sym.Matrix([udot[0:3,0], omega[0:3,0], udot[3:6,0], omega[3:6,0]])
eddot = sym.Matrix([uddot[0:3,0], alpha[0:3,0], uddot[3:6,0], alpha[3:6,0]])
# e = u
# edot = udot
# eddot = uddot
display(s)
display([e,edot,eddot])


# ### Needed Matrix Quantities

# In[ ]:

# angular velocity for 2D case is trivial
omega1_skew = Skew_sym(omega[0:3,0])
omega2_skew = Skew_sym(omega[3:6,0])
alpha1_skew = Skew_sym(alpha[0:3,0])
alpha2_skew = Skew_sym(alpha[3:6,0])

# Rotation for a 2D problem
R1 = Rotate_sym(theta[2])
R2 = Rotate_sym(theta[5])
# display(R1,R2)

# "spatial" rotation matrix
R = sym.Matrix.vstack(sym.Matrix.hstack(R1,sym.zeros(3)),                       sym.Matrix.hstack(sym.zeros(3),R2))

# "spatial" angular velocity matrix
Omega_skew = sym.Matrix.vstack(sym.Matrix.hstack(omega1_skew,sym.zeros(3)),                                sym.Matrix.hstack(sym.zeros(3),omega2_skew))

# "spatial" angular acceleration matrix
Alpha_skew = sym.Matrix.vstack(sym.Matrix.hstack(alpha1_skew,sym.zeros(3)),                                sym.Matrix.hstack(sym.zeros(3),alpha2_skew))
# display(Omega)


# ### Define Kinematics

# In[ ]:

# Define velocity of element endpoints (nodes)
v = sym.simplify(udot + R*Omega_skew*s)
print('v = ')
display(v)

# Define acceleration of element endpoints (nodes)
a = sym.simplify(uddot + R*Omega_skew*Omega_skew*s + R*Alpha_skew*s)
print('\na = ')
display(a)


# ### Compute the Mass Matrix

# In[ ]:

# Define shape function for element with one node at each end
h = sym.symarray('h', 2)

h[0] = sym.Rational(1,2)*(1 - x)
h[1] = sym.Rational(1,2)*(1 + x)

# Compute shape function matrix
H = sym.expand(sym.Matrix([h[0]*sym.eye(3), h[1]*sym.eye(3)])).T
print('\nH = ')
display(H.T)

# Define velocity of any point 
Vp = H*v
print('\nV = ')
display(Vp)

# Define velocity of any point 
Ap = H*a
# print('\nA = ')
# display(Ap)

# Compute partial velocities of the nodes
Vr = sym.Matrix([[sym.diff(Vp[j],edot[i]) 
                         for j in range(len(Vp))] for i in range(len(edot))]).T
# v_r = H
print('\nVr = ')
display(Vr.T)
print(Vr.shape)

# Compute mass matrix
M = sym.Matrix([[sym.expand(sym.integrate(Vr[:,i].dot(Ap)*rho*A,(x,0,l))).coeff(eddot[j]) for i in range(len(eddot))]
                for j in range(len(eddot))])
print('\nM = \n')
display(M)
pickle.dump( M, open( "mass-matrix.dump", "wb" ) )


# In[ ]:

# sym.factor(M,A*rho*l)


# ### Compute Internal forces 

# #### Compute Longitudinal Strain Energy

# In[ ]:




# #### Compute Transverse Strain Energy

# In[ ]:




# #### Compute Internal Forces $Q_e = \frac{\partial U}{\partial e}$

# In[ ]:




# ### Applied and body force vector

# In[ ]:



