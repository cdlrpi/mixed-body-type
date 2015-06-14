# Jeremy Lafin 6/8/2015
# Modified version of MBStructs.py by Mike Hans
# In this version 
#   - ANCF elements and 
#   - GEBF elements have been added
#   ** Everything is 2D planar i.e., Aij = [theta_ddot,a1,a2] **
#   ** This is because GEBF generates Mii = 0  for unused dofs **

import numpy as np
import sympy as sym
import pickle
import MultiBodyFuncts as MBF

from numpy.linalg import inv

# Rigid_Body is the class that defines a rigid body
class Rigid_Body2D:

    def initialize(self, m, l, I):
        self.m = m
        self.l = l
        self.I = I

    def intProps(self):
        # Locate Handles from C.M. in Body-Fixed Framen and rotate to 'N'
        self.r01 = np.dot(self.CBN,np.array([-self.l/2,0]))
        self.r02 = np.dot(self.CBN,np.array([ self.l/2,0]))
        
        # Applied and body forces
        # Only applied or body force is gravity
        Fg = np.array([0,0,self.m*-9.81])
        self.Fa1 = Fg 
        self.Fa2 = Fg 

        # construct shifter matricies for the **PLANAR** problem
        self.S01 = np.eye(3) + np.vstack((np.array([[0, -self.r01[1], self.r01[0]]]),np.zeros((2,3))))
        self.S02 = np.eye(3) + np.vstack((np.array([[0, -self.r02[1], self.r02[0]]]),np.zeros((2,3))))

        # Create the Mass and inverse Mass Matrices	
        self.M = np.vstack(((np.array([self.I,0,0]),np.hstack((np.zeros((2,1)),(self.m*np.eye(2)))))))
        self.Minv = np.linalg.inv(self.M)
        
        # Construct the inverse inertial matricies for the two handle equations
        self.z11 = np.dot(self.S01.T,np.dot(self.Minv,self.S01))
        self.z12 = np.dot(self.S01.T,np.dot(self.Minv,self.S02))
        self.z21 = np.dot(self.S02.T,np.dot(self.Minv,self.S01))
        self.z22 = np.dot(self.S02.T,np.dot(self.Minv,self.S02))

        self.z13 = np.dot(self.S01.T,np.dot(self.Minv,self.Fa1)) - \
                self.omega**2*np.hstack((0,self.r01)).reshape(self.Fa1.shape) 
        self.z23 = np.dot(self.S02.T,np.dot(self.Minv,self.Fa2)) - \
                self.omega**2*np.hstack((0,self.r02)).reshape(self.Fa2.shape) 

# Joint is the class that defines a single 
# kinematic joint in a multibody system		
class Joint():
    def __init__(self,jointType):
        if jointType == 'revolute2D':
            revolute2D.__init__(self)
        elif jointType == 'gebf':
            gebf.__init__(self)

class revolute2D(Joint):
    def __init__(self):
        self.P = np.array([[1], [0], [0]])
        self.D = np.array([[0, 0], [1, 0], [0, 1]])

class gebf(Joint):
    def __init__(self):
        self.P = np.array([[0], [0], [0]])
        self.D = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

class GEBF_Element2D():
    def __init__(self, modulus, inertia, radius, density, length, state):
        """
        Initializes all of the inertial propertias of the element and assigns 
        values to physical and material properties 
        """
        # symbolic system parameters 
        E, I, r, rho, l, = sym.symbols('E I r rho l')

        # Load symbolic mass matrix and substitute material properties
        M = pickle.load( open( "gebf-mass-matrix.dump", "rb" ) ).subs([ \
            (E, modulus), (I, inertia), (r, radius), (rho, density), (l, length)])
        
        # load the body and applied force vector and substitute material properties
        self.beta = pickle.load( open( "gebf-beta.dump", "rb" ) ).subs([ \
            (E, modulus), (I, inertia), (r, radius), (rho, density), (l, length)])
        
        # Substitute State Variables
        # re-make symbols for proper substitution and create the paird list
        q = sym.Matrix(sym.symarray('q',2*8))
        q_sub = [(q, qi) for q, qi in zip(q, state)]
       
        # Substitute state variables
        beta = self.beta.subs(q_sub)
        M = M.subs(q_sub)

        # Form the binary-DCA algebraic quantities
        # Partition mass matrix
        M11 = np.array(M[0:4,0:4])
        M12 = np.array(M[0:4,4:8])
        M21 = np.array(M[4:8,0:4])
        M22 = np.array(M[4:8,4:8])

        # For now use these definitions to cast Fic (constraint forces between GEBF elements) 
        # into generalized constraint forces
        self.gamma11 = np.eye(4)
        self.gamma12 = np.zeros((4,4))
        self.gamma22 = np.eye(4)
        self.gamma21 = np.zeros((4,4))
        
        # partition beta into lambda13 and lambda23
        self.gamma13 = np.array(beta[0:4])
        self.gamma23 = np.array(beta[4:8])


        # Commonly inverted quantities
        iM11 = inv(self.M11)
        iM22 = inv(self.M22)
        Gamma1 = inv(self.M11 - self.M12*iM22*self.M21)
        Gamma2 = inv(self.M22 - self.M21*iM11*self.M12)

        # Compute all terms of the two handle equations
        self.z11 = Gamma1.dot(self.gamma11 - self.M12.dot(iM22.dot(self.gamma21)))
        self.z12 = Gamma1.dot(self.gamma12 - self.M12.dot(iM22.dot(self.gamma22)))
        self.z13 = Gamma1.dot(self.gamma13 - self.M12.dot(iM22.dot(self.gamma23)))

        self.z21 = Gamma2.dot(self.gamma21 - self.M21.dot(iM11.dot(self.gamma11)))
        self.z22 = Gamma2.dot(self.gamma22 - self.M21.dot(iM11.dot(self.gamma12)))
        self.z23 = Gamma2.dot(self.gamma23 - self.M21.dot(iM11.dot(self.gamma13)))

# class ANCF_Element2D():
#     def __init__(self, area, modulus, inertia, density, length, state):
#         """
#         Initializes all of the inertial propertias of the element and 
#         assigns values to physical and material properties 
#         """
# 
#         # re-make symbols for proper substitution
#         e = sym.Matrix(sym.symarray('e',8))
#         # symbolic system parameters 
#         E, I, A, rho, l, = sym.symbols('E I A rho l')
# 
#         # Load symbolic mass matrix
#         M = pickle.load( open( "ancf-mass-matrix.dump", "rb" ) ).subs([(A, area), (l, length), (rho, density)])
#         
#         # load the body and applied force vector (this still has sympolic 'e' values)
#         self.beta = pickle.load( open( "ancf-force-vector.dump", "rb" ) ).subs([(E, modulus), (A, area), \
#                                                                            (I, inertia), (rho, density), (l, length)])
# 
#         # Partition mass matrix
#         self.M11 = np.array(M[0:4,0:4])
#         self.M12 = np.array(M[0:4,4:8])
#         self.M21 = np.array(M[4:8,0:4])
#         self.M22 = np.array(M[4:8,4:8])
# 
#         # For now 
#         self.lambda11 = np.eye(4)
#         self.lambda12 = np.zeros((4,4))
#         self.lambda22 = np.eye(4)
#         self.lambda21 = np.zeros((4,4))
#         
#         # fully numberic (initial) values for body and applied forces
#         e_sub = [(e, ei) for e, ei in zip(e, state)]
#         beta = self.beta.subs(e_sub)
# 
#         # partition beta into lambda13 and lambda23
#         self.lambda13 = np.array(beta[0:4])
#         self.lambda23 = np.array(beta[4:8])
# 
# 
#         # Commonly inverted quantities
#         Gamma = inv(self.M11 - self.M12*inv(self.M22)*self.M21)
#         iM22 = inv(self.M22)
# 
#         # Compute all terms of the two handle equations
#         self.zeta11 = Gamma.dot(self.lambda11 - self.M12.dot(iM22.dot(self.lambda21)))
#         self.zeta12 = Gamma.dot(self.lambda12 - self.M12.dot(iM22.dot(self.lambda22)))
#         self.zeta13 = Gamma.dot(self.lambda13 - self.M12.dot(iM22.dot(self.lambda23)))
# 
#         self.zeta21 = iM22.dot(self.lambda21 - self.M21.dot(self.lambda11))
#         self.zeta22 = iM22.dot(self.lambda22 - self.M21.dot(self.lambda12))
#         self.zeta23 = iM22.dot(self.lambda23 - self.M21.dot(self.lambda13))
# 
#     def Update(self):
#         """
#         This function updated the body-and-applied force vectors as the state variables change.
#         """
#         # re-make symbols for proper substitution
#         e = sym.Matrix(sym.symarray('e',8))
#         e_sub = [(e, ei) for e, ei in zip(e, self.state)]
# 
#         # load the body and applied force vector and make substitutions for 
#         # physical and material properties.
#         beta = self.beta.subs(e_sub)
# 
#         # partition beta into lambda13 and lambda23
#         self.lambda13 = beta[0:4]
#         self.lambda23 = beta[4:8]
# 
