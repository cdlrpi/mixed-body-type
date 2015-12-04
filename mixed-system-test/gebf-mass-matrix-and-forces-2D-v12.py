
# coding: utf-8

# # Beam Element - This notebook computes the mass-matrix and body-forces for later use

# In[1]:

import numpy as np
import scipy as sp
import sympy as sym
import pickle

from scipy import linalg
from sympy import mpmath
from sympy import cos, sin
from sympy import lambdify

from IPython.display import display
from __future__ import division
from sympy.interactive import printing
printing.init_printing(use_latex=True)
np.set_printoptions(precision=4,suppress=True)
from sympy.interactive import printing


import matplotlib.pyplot as plt
get_ipython().magic('matplotlib inline')


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

# display([q,qdot,qddot])
# display([u,udot,uddot])


# Define shape function for element with one node at each end
h = sym.symarray('h', 2)
h[0] = (1 - x/l)
h[1] = (x/l)
# h[0] = sym.Rational(1,2)*(1-x)
# h[1] = sym.Rational(1,2)*(1+x)


# Compute shape function matrix
H = sym.Matrix([h[0]*sym.eye(2), h[1]*sym.eye(2)]).T
dHdx = H.diff(x)
print('\nH = ')
display(H)


# In[6]:

delta


# ### Needed Matrix Quantities

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
print('\nvp = ')
display(vP)
v1 = H.subs(x,0)*deltadot
# print('\nv1 = ')
# display(v1)
v2 = H.subs(x,l)*deltadot
# print('\nv2 = ')
# display(v2)

# Angular velocites of the cross-sections of the nodes
omega1 = sym.Matrix([omega[0]])
omega2 = sym.Matrix([omega[1]])

# Create 'Spatial' velocities for the nodes (useful for constraint force coefficient matrix)
V1 = sym.Matrix.vstack(omega1,v1)
V2 = sym.Matrix.vstack(omega2,v2)

# Define acceleration of generic point on element
aP = H*deltaddot + H*sym.Matrix.vstack(alpha1_skew,alpha2_skew)*R_interp*s +                H*sym.Matrix.vstack(omega1_skew*omega1_skew,omega2_skew*omega2_skew)*R_interp*s
print('\nap = ')
display(aP)


# ### Compute Partial Velocities

# In[9]:

# generalized speeds need to be omega_i not qdot_i 
gen_speed = sym.Matrix.vstack(sym.Matrix.vstack(sym.Matrix([omega[0]]), sym.Matrix(qdot[1:3])), 
    sym.Matrix.vstack(sym.Matrix([omega[1]]), sym.Matrix(qdot[4:6])))

# Compute partial velocities of point 'P' and the nodes
VrP = [sym.Matrix([sym.diff(v,u) for v in vP]) for u in gen_speed]
# print('\nVrP = ')
# display(VrP)
Vr1 = [sym.Matrix([sym.diff(v,u) for v in V1]) for u in gen_speed]
# print('\nVr1 = ')
# display(Vr1)
Vr2 = [sym.Matrix([sym.diff(v,u) for v in V2]) for u in gen_speed]
# print('\nVr2 = ')
# display(Vr2)


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
# display(kappa)


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

# In[14]:

"""
Note: Sympy bug! Integrating a matrix returns a vector!!!
"""
# Transverse strain energy
kappa_squared = (kappa.T*dHdx.T).dot(dHdx*kappa)
Ut = 1/2*sym.integrate(E*I*kappa_squared, (x,0,l))


# In[15]:

# stiffness matrix in 'N'  
# C = E*A*(R_interp*H*H.T*R_interp.T)
# C = E*A*sym.eye(2)
# Ul = 1/2*sym.simplify(sym.integrate(epsilon.T.dot(C*epsilon), (x,0,l)))
G = E/2.6
C = sym.Matrix([[2/3*E*A, 0],[0, 2/3*5/6*G*A]])
Ul = 1/2*sym.integrate(epsilon.T*R_interp*C*R_interp.T*epsilon, (x,0,l))[0]


# In[16]:

# Compute Total Energy
U = Ul + Ut

# Compute Internal Force Vector
gen_coords = sym.Matrix.vstack(sym.Matrix.vstack(sym.Matrix([theta[0]]), sym.Matrix(q[1:3])), 
    sym.Matrix.vstack(sym.Matrix([theta[1]]), sym.Matrix(q[4:6])))
Qe = sym.Matrix([-U.diff(qi) for qi in q])


# ### Applied and body force vector

# In[17]:

# Applied forces
# Gravity body force
fg = g*sym.Matrix([0,-1])

# Compute beta
beta = sym.Matrix([sym.simplify(sym.integrate(rho*vr.dot(fg),('r_2',-r,r),(x,0,l))) for vr in VrP]) + Qe + RHS


# ### Compute Constraint Force Coefficient Matricies

# In[18]:

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


# In[19]:

# display(Gamma1)
# display(Gamma2)


# In[20]:

# Compute Potential Energy
height = h[0]*delta[1] + h[1]*delta[2]
PE = sym.integrate(rho*g*height,('r_2',-r,r),(x,0,l))
# PE = sym.integrate(rho*g*height,(x,0,l))

# Compute Kinetic Energy
# V = sym.Matrix.vstack(V1,V2)
# KE = 1/2*V.T.dot(M*V)
KE = 1/2*sym.integrate(rho*vP.dot(vP),('r_2',-r,r),(x,0,l))


# ### Output

# In[21]:

pickle.dump( M,      open( "gebf-mass-matrix.dump",   "wb" ) )
pickle.dump( beta,   open( "gebf-force-vector.dump",  "wb" ) )
pickle.dump( Gamma1, open( "gebf-1c-matrix.dump",     "wb" ) )
pickle.dump( Gamma2, open( "gebf-2c-matrix.dump",     "wb" ) )
pickle.dump( U,      open( "gebf-strain-energy.dump", "wb" ) )
pickle.dump( PE,     open( "potential_enrgy",         "wb" ) ) 
pickle.dump( KE,     open( "gebf-KE.dump",            "wb" ) ) 


# In[22]:

M_func      = lambdify((E, A, I, r, rho, l, g, q, omega),      M, "numpy")
beta_func   = lambdify((E, A, I, r, rho, l, g, q, omega),   beta, "numpy")
Gamma1_func = lambdify((E, A, I, r, rho, l, g, q, omega), Gamma1, "numpy")
Gamma2_func = lambdify((E, A, I, r, rho, l, g, q, omega), Gamma2, "numpy")
# M_func      = lambdify((E, A, I, r, rho, l, g, q, theta, omega),      M, "numpy")
# beta_func   = lambdify((E, A, I, r, rho, l, g, q, theta, omega),   beta, "numpy")
# Gamma1_func = lambdify((E, A, I, r, rho, l, g, q, theta, omega), Gamma1, "numpy")
# Gamma2_func = lambdify((E, A, I, r, rho, l, g, q, theta, omega), Gamma2, "numpy")
# U_func    = lambdify((E, A, I, r, rho, l, g, pi, q, theta, omega),    U, "numpy")


# In[23]:

# # Debugging functions to trace source of error 
Qe_func    = lambdify((E, A, I, r, rho, l, g, q, omega), Qe, "numpy")
# Fg_func    = lambdify((E, A, I, r, rho, l, g, q, omega), Fg, "numpy")
Ut_func    = lambdify((E, A, I, r, rho, l, g, q, omega), Ut, "numpy")
Ul_func    = lambdify((E, A, I, r, rho, l, g, q, omega), Ul, "numpy")
# Qe_func    = lambdify((E, A, I, r, rho, l, g, q, theta, omega), Qe, "numpy")
# Fg_func    = lambdify((E, A, I, r, rho, l, g, q, theta, omega), Fg, "numpy")
# Ut_func    = lambdify((E, A, I, r, rho, l, g, q, theta, omega), Ut, "numpy")
# Ul_func    = lambdify((E, A, I, r, rho, l, g, q, theta, omega), Ul, "numpy")


# In[24]:

def test_func(q,omega,g_num):
    M_num      = M_func(     0.7e6, 0.02393, 1.215e-8, 0.02393, 5540, 0.12, g_num, q, omega)
    beta_num   = beta_func(  0.7e6, 0.02393, 1.215e-8, 0.02393, 5540, 0.12, g_num, q, omega)        
    Gamma1_num = Gamma1_func(0.7e6, 0.02393, 1.215e-8, 0.02393, 5540, 0.12, g_num, q, omega)
    Gamma2_num = Gamma2_func(0.7e6, 0.02393, 1.215e-8, 0.02393, 5540, 0.12, g_num, q, omega)

    M11 = np.array(M_num[0:3,0:3])
    M12 = np.array(M_num[0:3,3:6])
    M21 = np.array(M_num[3:6,0:3])
    M22 = np.array(M_num[3:6,3:6])

    # For now use these definitions to cast Fic (constraint forces between GEBF elements) 
    # into generalized constraint forces
    gamma11 = Gamma1_num[0:3,:]
    gamma21 = Gamma1_num[3:6,:]
    gamma12 = Gamma2_num[0:3,:]
    gamma22 = Gamma2_num[3:6,:]

    # partition beta into lambda13 and lambda23
    gamma13 = np.array(beta_num[0:3])
    gamma23 = np.array(beta_num[3:6])

    # Commonly inverted quantities
    iM11 = np.linalg.inv(M11)
    iM22 = np.linalg.inv(M22)
    Chi1 = np.linalg.inv(M11 - M12.dot(iM22.dot(M21)))
    Chi2 = np.linalg.inv(M22 - M21.dot(iM11.dot(M12)))

    # Compute all terms of the two handle equations
    z11 = Chi1.dot(gamma11 - M12.dot(iM22.dot(gamma21)))
    z22 = Chi2.dot(gamma22 - M21.dot(iM11.dot(gamma12)))
    z12 = Chi1.dot(gamma12 - M12.dot(iM22.dot(gamma22)))
    z21 = Chi2.dot(gamma21 - M21.dot(iM11.dot(gamma11)))
    
    z13 = Chi1.dot(gamma13 - M12.dot(iM22.dot(gamma23))).reshape((3,1))
    z23 = Chi2.dot(gamma23 - M21.dot(iM11.dot(gamma13))).reshape((3,1))

    # Cantilever conditions
    # fix node one, D = eye, P = zeros, and let node two be free 
    F1c = np.linalg.inv(z11).dot(-z13)
    A2  = np.dot(z21,F1c) + z23

#     # Free-Free conditions
#     A1 = z13
#     A2 = z23

    return F1c, A2
#     return Gamma2


# In[25]:

g_num = 0

phi = np.linspace(0,2*3.14,200)
delta = np.linspace(0,0.06,200)

q0GEBF = [np.array([-np.pi/2,0,0,0,0,-d]).reshape(6,1) for d,angle in zip(delta,phi)]
# angle = [np.array([-np.pi/2,-np.pi/2]) for val in phi]

# Ul_theta    = np.array([  Ul_func(0.7e6, 0.0018, 1.215e-8, 0.02393, 5540, 0.12, g_num, q0, theta0, np.zeros_like(omega)) for q0,theta0 in zip(q0GEBF, angle)])
# Ut_theta    = np.array([  Ut_func(0.7e6, 0.0018, 1.215e-8, 0.02393, 5540, 0.12, g_num, q0, theta0, np.zeros_like(omega)) for q0,theta0 in zip(q0GEBF, angle)])
# Qe_theta    = np.array([  Qe_func(0.7e6, 0.0018, 1.215e-8, 0.02393, 5540, 0.12, g_num, q0, theta0, np.zeros_like(omega)) for q0,theta0 in zip(q0GEBF, angle)])                                                                  
# beta_theta  = np.array([beta_func(0.7e6, 0.0018, 1.215e-8, 0.02393, 5540, 0.12, g_num, q0, theta0, np.zeros_like(omega)) for q0,theta0 in zip(q0GEBF, angle)])                                                                  
# # Fg_theta    = np.array([  Fg_func(0.7e6, 0.0018, 1.215e-8, 0.02393, 5540, 0.12, g_num, q0, theta0, np.zeros_like(omega)) for q0,theta0 in zip(q0GEBF, angle)])                                                                   
# # M_theta     = np.array([   M_func(0.7e6, 0.0018, 1.215e-8, 0.02393, 5540, 0.12, g_num, q0, theta0, np.zeros_like(omega)) for q0,theta0 in zip(q0GEBF, angle)])                                                                   
# sol_theta   = np.array([test_func(q0, theta0, np.zeros_like(omega), g_num) for q0,theta0 in zip(q0GEBF, angle)])
Ul_theta    = np.array([  Ul_func(0.7e6, 0.02393, 1.215e-8, 0.02393, 5540, 0.12, g_num, q0, np.zeros_like(omega)) for q0 in q0GEBF])
Ut_theta    = np.array([  Ut_func(0.7e6, 0.02393, 1.215e-8, 0.02393, 5540, 0.12, g_num, q0, np.zeros_like(omega)) for q0 in q0GEBF])
Qe_theta    = np.array([  Qe_func(0.7e6, 0.02393, 1.215e-8, 0.02393, 5540, 0.12, g_num, q0, np.zeros_like(omega)) for q0 in q0GEBF])                                                                  
beta_theta  = np.array([beta_func(0.7e6, 0.02393, 1.215e-8, 0.02393, 5540, 0.12, g_num, q0, np.zeros_like(omega)) for q0 in q0GEBF])                                                                  
# Fg_theta    = np.array([  Fg_func(0.7e6, 0.0018, 1.215e-8, 0.02393, 5540, 0.12, g_num, q0, theta0, np.zeros_like(omega)) for q0,theta0 in zip(q0GEBF, angle)])                                                                   
# M_theta     = np.array([   M_func(0.7e6, 0.0018, 1.215e-8, 0.02393, 5540, 0.12, g_num, q0, theta0, np.zeros_like(omega)) for q0,theta0 in zip(q0GEBF, angle)])                                                                   
sol_theta   = np.array([test_func(q0, np.zeros_like(omega), g_num) for q0 in q0GEBF])



# In[26]:

F1c = sol_theta[:,0]
A2 = sol_theta[:,1]


# In[27]:

font = {'family' : 'serif',
        'weight' : 'bold',
        'size'   : 16}

plt.rc('font', **font)


# In[28]:

plt.rc('text', usetex=True)
plt.plot(delta/0.12, Ul_theta)
plt.title('Axial')
plt.ylabel('Strain Energy')
plt.xlabel('$\Delta/l$')
plt.show()


# In[29]:

plt.rc('text', usetex=True)
plt.plot(delta/0.12, Ut_theta)
plt.title('Bending')
plt.ylabel('Strain Energy')
plt.xlabel('$\Delta/l$')
plt.show()


# In[30]:

plt.rc('text', usetex=True)
for i in range(6):
    plt.plot(delta/0.12, beta_theta[:,i])
plt.ylabel('Constraint force')
plt.xlabel('$\Delta/l$')
plt.legend(['1', '2', '3','4','5','6'])
plt.show()


# In[31]:

plt.rc('text', usetex=True)
for i in range(6):
    plt.plot(delta/0.12, Qe_theta[:,i])
plt.ylabel('$Q_{e}$')
plt.xlabel('$\Delta/l$')
plt.legend(['1', '2', '3','4','5','6'])
plt.show()


# In[32]:

plt.rc('text', usetex=True)
for i in range(3):
    plt.plot(delta/0.12, F1c[:,i])
plt.ylabel('Constraint force')
plt.xlabel('$\Delta/l$')
plt.legend(['1', '2', '3'])
plt.show()


# In[33]:

plt.rc('text', usetex=True)
for i in range(3):
    plt.plot(delta/0.12, A2[:,i])
plt.ylabel('Acceleration of Node-2')
plt.xlabel('$\Delta/l$')
plt.legend(['1', '2', '3'])
plt.show()


# In[34]:

plt.rc('text', usetex=True)
plt.plot(delta, F1c[:,2] - (delta.reshape(200,1)*0.02393*0.7e6/0.12))
plt.ylabel('Constraint force')
plt.xlabel('$\Delta/l$')
plt.legend(['1', '2', '3'])
plt.show()


# In[35]:

q_test = np.array([-np.pi/2, 0, 0, 0, 0, -delta[2]])
omega_test = np.zeros_like(omega)

M_num      = M_func(     0.7e6, 0.02393, 1.215e-8, 0.02393, 5540, 0.12, g_num, q_test, omega_test)
beta_num   = beta_func(  0.7e6, 0.02393, 1.215e-8, 0.02393, 5540, 0.12, g_num, q_test, omega_test)        
Gamma1_num = Gamma1_func(0.7e6, 0.02393, 1.215e-8, 0.02393, 5540, 0.12, g_num, q_test, omega_test)
Gamma2_num = Gamma2_func(0.7e6, 0.02393, 1.215e-8, 0.02393, 5540, 0.12, g_num, q_test, omega_test)

M11 = np.array(M_num[0:3,0:3])
M12 = np.array(M_num[0:3,3:6])
M21 = np.array(M_num[3:6,0:3])
M22 = np.array(M_num[3:6,3:6])

# For now use these definitions to cast Fic (constraint forces between GEBF elements) 
# into generalized constraint forces
gamma11 = Gamma1_num[0:3,:]
gamma21 = Gamma1_num[3:6,:]
gamma12 = Gamma2_num[0:3,:]
gamma22 = Gamma2_num[3:6,:]

# partition beta into lambda13 and lambda23
gamma13 = np.array(beta_num[0:3])
gamma23 = np.array(beta_num[3:6])

# Commonly inverted quantities
iM11 = np.linalg.inv(M11)
iM22 = np.linalg.inv(M22)
Chi1 = np.linalg.inv(M11 - M12.dot(iM22.dot(M21)))
Chi2 = np.linalg.inv(M22 - M21.dot(iM11.dot(M12)))

# Compute all terms of the two handle equations
z11 = Chi1.dot(gamma11 - M12.dot(iM22.dot(gamma21)))
z22 = Chi2.dot(gamma22 - M21.dot(iM11.dot(gamma12)))
z12 = Chi1.dot(gamma12 - M12.dot(iM22.dot(gamma22)))
z21 = Chi2.dot(gamma21 - M21.dot(iM11.dot(gamma11)))

z13 = Chi1.dot(gamma13 - M12.dot(iM22.dot(gamma23))).reshape((3,1))
z23 = Chi2.dot(gamma23 - M21.dot(iM11.dot(gamma13))).reshape((3,1))


# In[36]:

np.linalg.inv(z11)


# In[37]:

data = np.hstack((delta.reshape(200,1),
                 F1c[:,2] - (delta.reshape(200,1)*0.02393*0.7e6/0.12)))
np.savetxt('data_def_force.mat', data, fmt='%3.16f')


# In[38]:

# epsilon


# In[39]:

1/12*(0.02393/0.12)**2


# In[40]:

# C = sym.Matrix([[E*A, 0],[0, E*A]])
# sym.simplify(sym.integrate(epsilon.T*R_interp*H*E*A*sym.eye(4)*H.T*R_interp.T*epsilon, (x,0,l)))


# In[ ]:




# In[41]:

sym.nsimplify(sym.simplify(KE))


# In[42]:

sym.nsimplify(sym.simplify(PE))

