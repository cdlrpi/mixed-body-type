
# coding: utf-8

# #Beam Element - This notebook computes the mass matrix and internal forces and outputs them for later use

# 
# 1. 2D beam generalized coordinates (1 rotation 2 displacements) (Reduced dimensions so all matricies are subsequently invertible, Mii = 0 if some dof's are accounted for but  unused)
# 2. Interpolation of roataion matricies WITHOUT orthogonalization (sympy struggling with too much algebra)
# 3. circular cross-section with radius 'r'

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

from IPython.display import display
from sympy.interactive import printing
printing.init_printing(use_latex='mathjax')
np.set_printoptions(precision=4,suppress=True)


# #### Define Needed Functions

# In[2]:

def Skew_sym(v):
    """
    This function returns the skew symetric matrix 
    of the vector 'v' to affect the cross product of 'v'x'u'
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
    R = sym.Matrix([[cos(theta),-sin(theta)],
                    [sin(theta), cos(theta)]])
    return R


# #### Define symbolic quantites
# symbolic system parameters 
E, G, I, A, rho, x, l, r, g  = sym.symbols('E G I A rho x l r g')

# Kinematic values of previos nodes (generic)
# e.g., omega_node  = omega + qdot
theta = sym.Matrix(['theta_1','theta_2'])
omega = sym.Matrix(['omega_1','omega_2'])
alpha = sym.Matrix(['alpha_1','alpha_2'])

# coordinates of the point in the 2D cross-section
# of nodes one and two 
s1 = sym.Matrix(['r_2','r_3'])
s2 = sym.Matrix(['r_2','r_3'])
s = sym.Matrix.vstack(s1,s2)

# generalized coordinates
# one rotation and two displacements per-node (two nodes per element)
# in this version generalzied speeds are qdots
q = sym.Matrix(sym.symarray('q',6))
qdot = sym.Matrix(sym.symarray('qdot',len(q)))
qddot = sym.Matrix(sym.symarray('qddot',len(q)))

# Deformations of Nodes (u's are not generalized speeds) 
u = sym.Matrix([q[1:3,0], q[4:6,0]])
udot = sym.Matrix([qdot[1:3,0], qdot[4:8,0]])
uddot = sym.Matrix([qddot[1:3,0], qddot[4:6,0]])

# ### Needed Matrix Quantities
""" 
Some cheating here: 
q0,q3 and q0dot,q3dot are really theta_1j, theta_2j and omega_1j, omega_2j
"""
# angular position and velocity for 2D using relative coordinates
# the sum of the respective quantites of bodies 1 - bodyj-1 
# Kinematics are trivial for the planar case:
# theta_j, omega_j =  sum(q1 ... qj), sum(q1dot ... qjdot), j = 1 ... nbodies

R1 = Rotate_sym(q[0])
R2 = Rotate_sym(q[3])

# Only true for planar case
omega1_skew = Skew_sym([0,0,qdot[0]])
omega2_skew = Skew_sym([0,0,qdot[3]])

# Only true for planar case
alpha1_skew = Skew_sym([0,0,qddot[0]])
alpha2_skew = Skew_sym([0,0,qddot[3]])


# "spatial" rotation matrix
R = sym.Matrix.vstack(sym.Matrix.hstack(R1,sym.zeros(2)),                       sym.Matrix.hstack(sym.zeros(2),R2))

# "spatial" angular velocity matrix
Omega_skew = sym.Matrix.vstack(sym.Matrix.hstack(omega1_skew,sym.zeros(2)),                                sym.Matrix.hstack(sym.zeros(2),omega2_skew))

# "spatial" angular acceleration matrix
Alpha_skew = sym.Matrix.vstack(sym.Matrix.hstack(alpha1_skew,sym.zeros(2)),                                sym.Matrix.hstack(sym.zeros(2),alpha2_skew))


# ### Define Kinematics
# Define Locations of Centroid of nodes
X0 = sym.Matrix(['0','0','l','0'])
rp = sym.simplify(X0 + u + R*s)

# Define velocity of element endpoints (nodes)
v = sym.simplify(udot + R*Omega_skew*s)

# Define acceleration of element endpoints (nodes)
a = sym.simplify(uddot + R*Omega_skew*Omega_skew*s + R*Alpha_skew*s)

# ### Compute the Mass Matrix
# Define shape function for element with one node at each end
h = sym.symarray('h', 2)

h[0] = sym.Rational(1,2)*(1 - x)
h[1] = sym.Rational(1,2)*(1 + x)

# Compute shape function matrix
H = sym.Matrix([h[0]*sym.eye(2), h[1]*sym.eye(2)]).T
dHdx = H.diff(x)

# Define velocity of any point 
Rp = H*rp

# Define velocity of any point 
Vp = H*v

# Define velocity of any point 
Ap = H*a

# Compute partial velocities of the nodes
Vr = sym.simplify(sym.Matrix([[sym.diff(Vp,qdot) for Vp in Vp] for qdot in qdot]).T)

# Compute mass matrix
M = sym.simplify(sym.factor(sym.Matrix(
            [[sym.expand(sym.integrate(Vr[:,i].dot(Ap)*rho,('r_2',0,r),('r_3',0,r),(x,0,l))).coeff(qddot[j]) 
              for i in range(len(qddot))] for j in range(len(qddot))])))



# ### Compute Internal forces 
# #### 1. Transverse (Bending) Strain

# Orthogonal Matricies Not Extracted to Simplify Algebra
R_interp = sym.simplify(H*sym.Matrix([R1,R2]))
dT = sym.simplify(H.diff(x)*sym.Matrix([R1,R2]))
kappa = sym.simplify(sym.Matrix([Axial_sym(dT*R_interp.T),'0','0','0']))

# #### 2. Longitudinal (Axial) Strian
# Define Locations of Centroid
x0_B = sym.Matrix(['x','0'])
x0 = R_interp*x0_B

# Define Newtonian Unit Vector x-dir
n1 = sym.Matrix(['1','0'])

# Interpolate Displacemnts
u_terp = H*u

# Derivatives w.r.t longitudinal beam coordinate
du = u_terp.diff(x)
dx0 = x0.diff(x)

# Compute axial strain
u_ax = (dx0 + du - R_interp*n1).simplify()
epsilon = sym.Matrix(['0', u_ax[0], '0', u_ax[1]])
# epsilon = u_ax
# display(epsilon)

# #### 3. Compute Internal Forces $Q_e = -\frac{\partial U}{\partial e}$
"""
Note: Sympy bug! Integrating a matrix returns a vector!!!
"""
# Transverse strain energy
kappa_squared = (kappa.T*dHdx.T).dot(dHdx*kappa)
Ut = 1/2*sym.integrate(E*I*kappa_squared, (x,0,l))

G = E/2.6
C = sym.Matrix([[E*A, 0],[0, 5/6*G*A]])
Ul = 1/2*sym.integrate(epsilon.T*dHdx.T*R_interp*C*R_interp.T*dHdx*epsilon, (x,0,l))[0]

# Compute Total Energy
U = Ul + Ut

# Compute Internal Force Vector
Qe = sym.Matrix([sym.simplify(sym.expand(-sym.diff(U,qi))) for qi in q])


# ####4. Applied and body force vector

# Applied forces
# Gravity body force
fg = g*rho*sym.Matrix([0,-1])

# Applied torque (not considered at this time, no partial angular velocities)
# torque_app = sym.Matrix([0,0,tau_a])

# Compute beta
beta = sym.Matrix([sym.simplify(sym.integrate(Vr[:,j].dot(fg),('r_2',0,r),('r_3',0,r),(x,0,l)))
                   + qe for j,qe in zip(range(len(q)),Qe)])

# Just for debugging purposes
Fg = sym.Matrix([sym.simplify(sym.integrate(Vr[:,j].dot(fg),('r_2',0,r),('r_3',0,r),(x,0,l)))
                   for j in range(len(q))])

pickle.dump( M,    open( "gebf-mass-matrix.dump",   "wb" ) )
pickle.dump( beta, open( "gebf-force-vector.dump",  "wb" ) )
pickle.dump( U,    open( "gebf-strain-energy.dump", "wb" ) )
# pickle.dump(PE, open("potential_enrgy", "wb")) 

M_func    = lambdify((E, A, I, r, rho, l, g, q),    M, "numpy")
beta_func = lambdify((E, A, I, r, rho, l, g, q), beta, "numpy")
U_func    = lambdify((E, A, I, r, rho, l, g, q),    U, "numpy")

# Debugging functions to trace source of error 
Qe_func    = lambdify((E, A, I, r, rho, l, g, q),    Qe, "numpy")
Fg_func    = lambdify((E, A, I, r, rho, l, g, q),    Fg, "numpy")
Ut_func    = lambdify((E, A, I, r, rho, l, g, q),    Ut, "numpy")
Ul_func    = lambdify((E, A, I, r, rho, l, g, q),    Ul, "numpy")

beta_num = beta_func(0.7e6, 0.0018, 1.215e-8, 0.02393, 5540, 0.12, 9.81, np.zeros_like(q))
M_num    = M_func(   0.7e6, 0.0018, 1.215e-8, 0.02393, 5540, 0.12, 9.81, np.zeros_like(q))
U_num    = U_func(   0.7e6, 0.0018, 1.215e-8, 0.02393, 5540, 0.12, 9.81, np.zeros_like(q))


theta = np.linspace(0,2*np.pi,200)
q0GEBF = [np.array([theta,0,0,0,0,0]).reshape(6,1) for theta in theta]
q0GEBF[:][3] = q0GEBF[:][0]
Ul_theta1 = np.array([Ul_func(0.7e6, 0.0018, 1.215e-8, 0.02393, 5540, 0.12, 9.81, q0) for q0 in q0GEBF])
Ut_theta1 = np.array([Ut_func(0.7e6, 0.0018, 1.215e-8, 0.02393, 5540, 0.12, 9.81, q0) for q0 in q0GEBF])
Qe_theta1 = np.array([np.linalg.norm(Qe_func(0.7e6, 0.0018, 1.215e-8, 0.02393, 5540, 0.12, 9.81, q0)) for q0 in q0GEBF])

q0GEBF = [np.array([theta,0,0,0,0,0]).reshape(6,1) for theta in theta]
# q0GEBF[:][3] = q0GEBF[:][0]
Ul_theta2 = np.array([Ul_func(0.7e6, 0.0018, 1.215e-8, 0.02393, 5540, 0.12, 9.81, q0) for q0 in q0GEBF])
Ut_theta2 = np.array([Ut_func(0.7e6, 0.0018, 1.215e-8, 0.02393, 5540, 0.12, 9.81, q0) for q0 in q0GEBF])
Qe_theta2 = np.array([np.linalg.norm(Qe_func(0.7e6, 0.0018, 1.215e-8, 0.02393, 5540, 0.12, 9.81, q0)) for q0 in q0GEBF])

theta.tofile('theta.data')
Ul_theta1.tofile('strainEnergyAxial1_3.data')
Ul_theta2.tofile('strainEnergyAxial2_3.data')
Ut_theta1.tofile('strainEnergyBending1_3.data')
Ut_theta2.tofile('strainEnergyBending2_3.data')
Qe_theta1.tofile('bodyForces1_3.data')
Qe_theta2.tofile('bodyForces2_3.data')
