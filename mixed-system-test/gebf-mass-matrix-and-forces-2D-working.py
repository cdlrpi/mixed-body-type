
# coding: utf-8

# # Beam Element - This notebook computes the mass-matrix and body-forces for later use

# In[1]:
from __future__ import division

import numpy as np
import scipy as sp
import sympy as sym
import pickle

from scipy import linalg
from sympy import mpmath
from sympy import cos, sin
from sympy import lambdify


# #### Define Needed Functions

# In[2]:

def Skew_sym(v):
    """
    This function returns the skew symetric matrix 
    of the angular velocity vector for a 2D problem
    """
    v_matrix = sym.Matrix([[  0 , -v[2]],
                          [v[2],     0]])
    return v_matrix


# In[3]:

def Axial_sym(A):
    '''
    This funtcion returns the vector of the skew-symmectric matrix in 2D
    '''
    a_vec = 1/2*sym.Matrix([A[1,0] - A[0,1]])
    return a_vec


# In[4]:

def Rotate_sym(theta):
    """
    This function returns the symbolic rotation matrix 
    for the simple 2D-rotation about the third axis
    """
    R = sym.Matrix([[sym.cos(theta),-sym.sin(theta)],
                    [sym.sin(theta), sym.cos(theta)]])
    return R


# #### Define symbolic quantites

# In[5]:

# symbolic system parameters 
E, G, I, A, rho, x, l, r, g, m, a, b, c, pi, phi  = sym.symbols('E G I A rho x l r g m a b c pi phi')

# theta = sym.Matrix(['theta_1','theta_2'])
theta_dot = sym.Matrix(['thetadot_1','thetadot_2'])
theta_ddot = sym.Matrix(['thetaddot_1','thetaddot_2'])
omega = sym.Matrix(['omega_1','omega_2'])
alpha = sym.Matrix(['alpha_1','alpha_2'])

# Symbolic constraint forces and moments at the handles
f1c = sym.Matrix(['f_1_1','f_1_2'])
f2c = sym.Matrix(['f_2_1','f_2_2'])
tau1c = sym.Matrix(['tau_1'])
tau2c = sym.Matrix(['tau_2'])

# coordinates of the point in the 2D cross-section
# of nodes one and two 
s = sym.Matrix(['r_1','r_2']) 
s[0] = 0

# generalized coordinates
# one rotation and two displacements per-node (two nodes per element)
# in this version generalzied speeds are NOT ALWAYS qdots
q = sym.Matrix(sym.symarray('q',6))
qdot = sym.Matrix(sym.symarray('qdot',len(q)))
qddot = sym.Matrix(sym.symarray('qddot',len(q)))

# Deformations of Nodes (u's are not generalized speeds) 
delta = sym.Matrix([q[1:3,0], q[4:6,0]])
deltadot = sym.Matrix([qdot[1:3,0], qdot[4:8,0]])
deltaddot = sym.Matrix([qddot[1:3,0], qddot[4:6,0]])

# Define shape function for element with one node at each end
h = sym.symarray('h', 2)
h[0] = (1 - x/l)
h[1] = (x/l)


# Compute shape function matrix
H = sym.Matrix([h[0]*sym.eye(2), h[1]*sym.eye(2)]).T
dHdx = H.diff(x)

# In[7]:

# Rotation Matricies for each node
theta = sym.Matrix([q[0], q[0] + q[3]])

R1 = Rotate_sym(theta[0])
R2 = Rotate_sym(theta[1])

# Angular Velocities and Accelerations are trivial for the 2D case
# For each node
omega1_skew = Skew_sym([0,0,omega[0]])
omega2_skew = Skew_sym([0,0,omega[1]])
alpha1_skew = Skew_sym([0,0,alpha[0]])
alpha2_skew = Skew_sym([0,0,alpha[1]])

R_interp = H*sym.Matrix.vstack(R1,R2)

# ### Define Kinematics

# In[8]:

# Define velocity of generic point on element and handles
vP = H*deltadot + H*sym.Matrix.vstack(omega1_skew,omega2_skew)*R_interp*s
v1 = H.subs(x,0)*deltadot
v2 = H.subs(x,l)*deltadot

# Angular velocites of the cross-sections of the nodes
omega1 = sym.Matrix([omega[0]])
omega2 = sym.Matrix([omega[1]])

# Create 'Spatial' velocities for the nodes (useful for constraint force coefficient matrix)
V1 = sym.Matrix.vstack(omega1,v1)
V2 = sym.Matrix.vstack(omega2,v2)

# Define acceleration of generic point on element
aP = H*deltaddot + H*sym.Matrix.vstack(alpha1_skew,alpha2_skew)*R_interp*s +   H*sym.Matrix.vstack(omega1_skew*omega1_skew,omega2_skew*omega2_skew)*R_interp*s


# ### Compute Partial Velocities

# In[9]:

# generalized speeds need to be omega_i not qdot_i 
gen_speed = sym.Matrix.vstack(sym.Matrix.vstack(sym.Matrix([omega[0]]), sym.Matrix(qdot[1:3])), 
    sym.Matrix.vstack(sym.Matrix([omega[1]]), sym.Matrix(qdot[4:6])))

# Compute partial velocities of point 'P' and the nodes
VrP = [sym.Matrix([sym.diff(v,u) for v in vP]) for u in gen_speed]
Vr1 = [sym.Matrix([sym.diff(v,u) for v in V1]) for u in gen_speed]
Vr2 = [sym.Matrix([sym.diff(v,u) for v in V2]) for u in gen_speed]


# ### Compute the Mass Matrix

# In[10]:

gen_accel = sym.Matrix.vstack(sym.Matrix.vstack(sym.Matrix([alpha[0]]), sym.Matrix(qddot[1:3])), 
    sym.Matrix.vstack(sym.Matrix([alpha[1]]), sym.Matrix(qddot[4:6])))

RHS = sym.Matrix([sym.integrate(rho*vr.dot(aP), ('r_2',-r,r), (x,0,l)) for vr in VrP])
# Compute mass matrix
M = sym.Matrix([[sym.simplify(sym.expand(rhs).coeff(uddot)) for uddot in gen_accel] for rhs in RHS])
RHS = sym.simplify(RHS - M*gen_accel)


# ### Compute Internal forces 

# #### 1. Transverse (Bending) Strain

# In[11]:

# Orthogonal Matricies Not Extracted to Simplify Algebra
dT = sym.simplify(H.diff(x)*sym.Matrix([R1,R2]))
kappa = sym.simplify(sym.Matrix([Axial_sym(dT*R_interp.T),'0','0','0']))


# #### 2. Longitudinal (Axial) Strian

# In[12]:

# Define Locations of Centroid as a function of beam axis coordinate
x0_B = sym.Matrix(['x','0'])
# Convert to the newtonian basis
x0 = R_interp*x0_B

# Define Newtonian Unit Vector x-dir
n1 = sym.Matrix(['1','0'])

# strain = du/dx
# Derivatives w.r.t longitudinal beam coordinate
d_delta = dHdx*delta
dx0 = x0.diff(x)

# Compute axial strain in 'N'
epsilon = dx0 + d_delta - R_interp*n1


# In[13]:

# epsilon


# #### 3. Compute Internal Forces $Q_e = -\frac{\partial U}{\partial e}$


"""
Note: Sympy bug! Integrating a matrix returns a vector!!!
"""
# Transverse strain energy
kappa_squared = (kappa.T*dHdx.T).dot(dHdx*kappa)
Ut = 1/2*sym.integrate(E*I*kappa_squared, (x,0,l))


# In[ ]:

# stiffness matrix in 'N'  
C = E*A*(R_interp*H*H.T*R_interp.T)
Ul = 1/2*sym.simplify(sym.integrate(epsilon.T.dot(C*epsilon), (x,0,l)))


# In[ ]:

# Compute Total Energy
U = Ul + Ut

# Compute Internal Force Vector
gen_coords = sym.Matrix.vstack(sym.Matrix.vstack(sym.Matrix([theta[0]]), sym.Matrix(q[1:3])), 
    sym.Matrix.vstack(sym.Matrix([theta[1]]), sym.Matrix(q[4:6])))
Qe = sym.Matrix([-U.diff(qi) for qi in q])


# ### Applied and body force vector

# In[ ]:

# Applied forces
# Gravity body force
fg = g*rho*sym.Matrix([0,-1])

# Compute beta
beta = sym.Matrix([sym.simplify(sym.integrate(vr.dot(fg),('r_2',-r,r),(x,0,l))) for vr in VrP]) + Qe + RHS


# ### Compute Constraint Force Coefficient Matricies

# In[ ]:

# Compute partial velocities of handles
r12 = R1*sym.Matrix([ l, 0])
r21 = R2*sym.Matrix([-l, 0])

S12 = sym.Matrix([[1, r12[1], r12[0]],[0, 1, 0],[0, 0, 1]])
S21 = sym.Matrix([[1, r21[1], r21[0]],[0, 1, 0],[0, 0, 1]])

F1c = sym.Matrix.vstack(tau1c,f1c)
F2c = sym.Matrix.vstack(tau2c,f2c)

Vr12 = [S12.T*vr for vr in Vr1]
Vr21 = [S21.T*vr for vr in Vr2]

Gamma1 = sym.Matrix([[sym.expand(vr.dot(F1c)).coeff(f) for f in F1c] for vr in Vr12])
Gamma2 = sym.Matrix([[sym.expand(vr.dot(F2c)).coeff(f) for f in F2c] for vr in Vr21])


# Compute Potential Energy
height = h[0]*delta[1] + h[1]*delta[2]
PE = sym.integrate(rho*g*height,('r_2',-r,r),(x,0,l))
KE = sym.integrate(rho*vP.dot(vP),('r_2',-r,r),(x,0,l))


# ### Output

# In[ ]:

pickle.dump( M,      open( "gebf-mass-matrix.dump",   "wb" ) )
pickle.dump( beta,   open( "gebf-force-vector.dump",  "wb" ) )
pickle.dump( Gamma1, open( "gebf-1c-matrix.dump",     "wb" ) )
pickle.dump( Gamma2, open( "gebf-2c-matrix.dump",     "wb" ) )
pickle.dump( U,      open( "gebf-strain-energy.dump", "wb" ) )
pickle.dump( PE,     open( "gebf-PE.dump",            "wb" ) ) 
pickle.dump( KE,     open( "gebf-KE.dump",            "wb" ) )


