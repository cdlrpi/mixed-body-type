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

class Joint():
    def __init__(self,jointType):
        if jointType == 'revolute2D':
            revolute2D.__init__(self)
        elif jointType == 'fixed':
            fixed.__init__(self)

class revolute2D(Joint):
    def __init__(self):
        self.P = np.array([[1], [0], [0]])
        self.D = np.array([[0, 0], [1, 0], [0, 1]])

class fixed(Joint):
    def __init__(self):
        self.P = np.array([[0], [0], [0]])
        self.D = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

class Body:

    def initialize(self, body_type, *args):
        if body_type == 'rigid':
            Rigid_Body2D.initialize(self, args)
        elif body_type == 'gebf':
            GEBF_Element2D.initialize(self, args)

    def intProps(self, body_type, *args):
        if body_type == 'rigid':
            Rigid_Body2D.intProps(self)
        elif body_type == 'gebf':
            GEBF_Element2D.intProps(self, args)

# Rigid_Body is the class that defines inertial properties of a rigid body
class Rigid_Body2D(Body):
    
    def initialize(self, args):
        self.m = args[0]
        self.l = args[1]
        self.I = args[2]
    
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


class GEBF_Element2D(Body):
    
    def initialize(self, args):
        self.E   = args[0]
        self.A   = args[1]
        self.I   = args[2]
        self.r   = args[3]
        self.rho = args[4]
        self.l   = args[5]

    def intProps(self,args):
       
        # u's are deformations not generalized speeds
        u = args[0]
       
        # re-make symbols for proper substitution
        # no velocity dependent terms 
        q = sym.Matrix(sym.symarray('q',6))
        E, A, I, r, rho, l, = sym.symbols('E A I r rho l')

        # Some cheating going on here with the kinematics and "position" vector 
        # q_e is the "positions" with q0,q3 = sum(q0j), sum(q3j) (j = 1 ... nbodies)
        # self.theta1 and self.theta2 are updated in the kinematic sweep
        q_e = [self.theta1] + u[:2] + [self.theta2] + u[2:]        

        # create the paird list for substitution
        q_sub = [(q, qi) for q, qi in zip(q, q_e)]

        # Load symbolic mass matrix and substitute material properties
        self.M = pickle.load( open( "gebf-mass-matrix.dump", "rb" ) ).\
        subs([(E, self.E), (A, self.A), (r, self.r), (I, self.I), (rho, self.rho), (l, self.l)]).\
        subs(q_sub).evalf()
        
        # load the body and applied force vector and substitute material properties
        self.beta = pickle.load( open( "gebf-force-vector.dump", "rb" ) ).\
        subs([(E, self.E), (A, self.A), (r, self.r), (I, self.I), (rho, self.rho), (l, self.l)]).\
        subs(q_sub).evalf()

        # Form the binary-DCA algebraic quantities
        # Partition mass matrix
        self.M11 = np.array(self.M[0:3,0:3])
        self.M12 = np.array(self.M[0:3,3:6])
        self.M21 = np.array(self.M[3:6,0:3])
        self.M22 = np.array(self.M[3:6,3:6])

        # For now use these definitions to cast Fic (constraint forces between GEBF elements) 
        # into generalized constraint forces
        self.gamma11 = np.eye(3)
        self.gamma12 = np.zeros((3,3))
        self.gamma22 = np.eye(3)
        self.gamma21 = np.zeros((3,3))
        
        # partition beta into lambda13 and lambda23
        self.gamma13 = np.array(self.beta[0:3])
        self.gamma23 = np.array(self.beta[3:6])


        # Commonly inverted quantities
        self.iM11 = inv(self.M11)
        self.iM22 = inv(self.M22)
        self.Gamma1 = inv(self.M11 - self.M12.dot(self.iM22.dot(self.M21)))
        self.Gamma2 = inv(self.M22 - self.M21.dot(self.iM11.dot(self.M12)))

        # Compute all terms of the two handle equations
        self.z11 = self.Gamma1.dot(self.gamma11 - self.M12.dot(self.iM22.dot(self.gamma21)))
        self.z12 = self.Gamma1.dot(self.gamma12 - self.M12.dot(self.iM22.dot(self.gamma22)))
        self.z21 = self.Gamma2.dot(self.gamma21 - self.M21.dot(self.iM11.dot(self.gamma11)))
        self.z22 = self.Gamma2.dot(self.gamma22 - self.M21.dot(self.iM11.dot(self.gamma12)))
        
        
        self.z13 = self.Gamma1.dot(self.gamma13 - self.M12.dot(self.iM22.dot(self.gamma23))).reshape((3,1))
        self.z23 = self.Gamma2.dot(self.gamma23 - self.M21.dot(self.iM11.dot(self.gamma13))).reshape((3,1))

    
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
#        # Load symbolic mass matrix
#        M = pickle.load( open( "ancf-mass-matrix.dump", "rb" ) ).subs([(A, area), (l, length), (rho, density)])
#
#        # load the body and applied force vector (this still has sympolic 'e' values)
#        beta = pickle.load( open( "ancf-force-vector.dump", "rb" ) ).subs([(E, modulus), (A, area), \
#                                                                           (I, inertia), (rho, density), (l, length)])
#        # Partition mass matrix
#        M11 = np.array([M[0:4,0:4]], dtype=np.double).reshape((4,4))
#        M12 = np.array([M[0:4,4:8]], dtype=np.double).reshape((4,4))
#        M21 = np.array([M[4:8,0:4]], dtype=np.double).reshape((4,4))
#        M22 = np.array([M[4:8,4:8]], dtype=np.double).reshape((4,4))
#
#        # For now 
#        lambda11 = np.eye(4)
#        lambda12 = np.zeros((4,4))
#        lambda21 = np.zeros((4,4))
#        lambda22 = np.eye(4)
#
#        # fully numberic (initial) values for body and applied forces
#
#        e_sub = [(e, ei) for e, ei in zip(e, state[:8])]
#        beta = beta.subs(e_sub)
#
#        # partition beta into lambda13 and lambda23
#        lambda13 = np.array([beta[0:4]], dtype=np.double).reshape((4,1))
#        lambda23 = np.array([beta[4:8]], dtype=np.double).reshape((4,1))
#
#        # Commonly inverted quantities
#        iM11 = inv(M11)
#        iM22 = inv(M22)
#
#        # Coefficients
#        Gamma1 = inv(M11 - M12.dot(iM22.dot(M21)))
#        Gamma2 = inv(M22 - M21.dot(iM11.dot(M12)))
#
#
#        # Compute all terms of the two handle equations
#        self.zeta11 = Gamma1.dot(lambda11 - M12.dot(iM22.dot(lambda21)))
#        self.zeta12 = Gamma1.dot(lambda12 - M12.dot(iM22.dot(lambda22)))
#        self.zeta13 = Gamma1.dot(lambda13 - M12.dot(iM22.dot(lambda23)))
#
#        self.zeta21 = Gamma2.dot(lambda21 - M21.dot(iM11.dot(lambda11)))
#        self.zeta22 = Gamma2.dot(lambda22 - M21.dot(iM11.dot(lambda12)))
#        self.zeta23 = Gamma2.dot(lambda23 - M21.dot(iM11.dot(lambda13)))
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
#         beta = beta.subs(e_sub)
