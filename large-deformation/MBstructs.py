# Jeremy Lafin 6/8/2015
# Modified version of MBStructs.py by Mike Hans
# In this version 
#   - ANCF elements and 
#   - GEBF elements have been added
#   - Ridid_Bodies have been modified to use __init__()
#   ** Everything is 2D planar i.e., Aij = [theta_ddot,a1,a2] **
#   ** This is because GEBF generates Mii = 0  for unused dofs **

import numpy as np
import sympy as sym
import pickle
import MultiBodyFuncts as MBF

from numpy.linalg import inv

# Body is the class that defines a single body in a multibody system
# Only really used for list management
class Body(self,body_type,*arg):
    if body_type == 'rigid':
        Rigid_Body(arg)
    elif body_type == 'GEBF':
        GEBF_Element(arg)
    #  ANCF Elements not used in this version


# Rigid_Body is the class that defines a rigid body
class Rigid_Body():
    def __init__(self, (m, l, I)):

        self.r10 = np.dot(np.array([l/2,0,0]),self.C0)
        self.r01 = -self.r10
        self.r12 = np.dot(np.array([l,0,0]),self.C0)
        self.r21 = -self.r12
        self.r20 = np.dot(np.array([-l/2,0,0]),self.C0)
        self.r02 = -self.r20

        # construct shifter matricies
        self.S10 = np.hstack((np.array([[1],[0],[0]]),np.vstack(([-r[1], r[0]],np.eye(2)))))
        self.S01 = np.hstack((np.array([[1],[0],[0]]),np.vstack(([r[1], -r[0]],np.eye(2)))))
        self.S02 = np.hstack((np.array([[1],[0],[0]]),np.vstack(([-r[1], r[0]],np.eye(2)))))
        self.S20 = np.hstack((np.array([[1],[0],[0]]),np.vstack(([r[1], -r[0]],np.eye(2)))))

        # Create the Mass and inverse Mass Matrices	
        self.M = np.hstack(((np.array([[I],[0],[0]]),np.vstack(np.zeros((2,1)),(np.eye(2)))))
        self.Minv=np.linalg.inv(self.M)
        
        # construct the inverse inertial matricies for the two handle equations
        self.z11=np.dot(self.S10,np.dot(self.Minv,self.S01))
        self.z12=np.dot(self.S10,np.dot(self.Minv,self.S02))
        self.z21=np.dot(self.S20,np.dot(self.Minv,self.S01))
        self.z22=np.dot(self.S20,np.dot(self.Minv,self.S02))

        # Determine Initial Body Forces and State Dependent Terms 
        self.i1= np.cross(self.w,self.r10)
        self.i12=np.cross(self.w,self.r20)
        self.i2=np.cross(self.w,self.i1)
        self.i31=self.m*self.i2
        self.i32=self.m*np.cross(self.w,self.i12)
        
        #Total force with gravity included
        self.Fa1[3:]=9.81*self.m*v-self.i31#-self.i5
        self.Fa2[3:]=9.81*self.m*v-self.i32

        self.z13=np.dot(self.S10,np.dot(self.Minv,self.Fa1))
        self.z23=np.dot(self.S20,np.dot(self.Minv,self.Fa2))

    #Function that defines position vectors between any of the
    #Three points on the body (two handles and center of mass)
    def rs(self,r,v):

    #Function that finds the initial zeta values in each timestep
    def zeta(self):

            
    #Function to apply the state dependant forces
    #I think this is where my problem is because I am 
    #not sure if this is correct.
    def Forces(self,v):


# Joint is the class that defines a single 
# kinematic joint in a multibody system		
class Joint:
    # Function to create the P and D matrices
    def __init__(self,v,x):
        """
        (simple joint)
        x = (1 or 0)  - (free or locked joint)
        v =  the vector defining the direction of allowed motion
        """
        if x == 0:
            self.P = np.zeros((6,1))
            self.D = np.eye(6)
        else:
            self.P=np.zeros((6,x))
            self.D=np.zeros((6,6-x))
            j=0
            k=0
            for i in range(0,6):
                if v[i]==1:
                    self.P[i,j]=1
                    j=j+1
                else:
                    self.D[i,k] = 1
                    k=k+1

class ANCF_Element():
    def __init__(self, area, modulus, inertia, density, length, state):
        """
        Initializes all of the inertial propertias of the element and 
        assigns values to physical and material properties 
        """

        # re-make symbols for proper substitution
        e = sym.Matrix(sym.symarray('e',8))
        # symbolic system parameters 
        E, I, A, rho, l, = sym.symbols('E I A rho l')

        # Load symbolic mass matrix
        M = pickle.load( open( "ancf-mass-matrix.dump", "rb" ) ).subs([(A, area), (l, length), (rho, density)])
        
        # load the body and applied force vector (this still has sympolic 'e' values)
        self.beta = pickle.load( open( "ancf-force-vector.dump", "rb" ) ).subs([(E, modulus), (A, area), \
                                                                           (I, inertia), (rho, density), (l, length)])

        # Partition mass matrix
        self.M11 = np.array(M[0:4,0:4])
        self.M12 = np.array(M[0:4,4:8])
        self.M21 = np.array(M[4:8,0:4])
        self.M22 = np.array(M[4:8,4:8])

        # For now 
        self.lambda11 = np.eye(4)
        self.lambda12 = np.zeros((4,4))
        self.lambda22 = np.eye(4)
        self.lambda21 = np.zeros((4,4))
        
        # fully numberic (initial) values for body and applied forces
        e_sub = [(e, ei) for e, ei in zip(e, state)]
        beta = self.beta.subs(e_sub)

        # partition beta into lambda13 and lambda23
        self.lambda13 = np.array(beta[0:4])
        self.lambda23 = np.array(beta[4:8])


        # Commonly inverted quantities
        Gamma = inv(self.M11 - self.M12*inv(self.M22)*self.M21)
        iM22 = inv(self.M22)

        # Compute all terms of the two handle equations
        self.zeta11 = Gamma.dot(self.lambda11 - self.M12.dot(iM22.dot(self.lambda21)))
        self.zeta12 = Gamma.dot(self.lambda12 - self.M12.dot(iM22.dot(self.lambda22)))
        self.zeta13 = Gamma.dot(self.lambda13 - self.M12.dot(iM22.dot(self.lambda23)))

        self.zeta21 = iM22.dot(self.lambda21 - self.M21.dot(self.lambda11))
        self.zeta22 = iM22.dot(self.lambda22 - self.M21.dot(self.lambda12))
        self.zeta23 = iM22.dot(self.lambda23 - self.M21.dot(self.lambda13))

    def Update(self):
        """
        This function updated the body-and-applied force vectors as the state variables change.
        """
        # re-make symbols for proper substitution
        e = sym.Matrix(sym.symarray('e',8))
        e_sub = [(e, ei) for e, ei in zip(e, self.state)]

        # load the body and applied force vector and make substitutions for 
        # physical and material properties.
        beta = self.beta.subs(e_sub)

        # partition beta into lambda13 and lambda23
        self.lambda13 = beta[0:4]
        self.lambda23 = beta[4:8]

class GEBF_Element():
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
        self.M11 = np.array(M[0:4,0:4])
        self.M12 = np.array(M[0:4,4:8])
        self.M21 = np.array(M[4:8,0:4])
        self.M22 = np.array(M[4:8,4:8])

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

        # # Make "full-size" to interface with "recursive" DCA for RigidBody() class
        # self.z11 = np.hstack((np.zeros((6,2)),np.vstack((np.zeros((2,4)),self.z11))))
        # self.z12 = np.hstack((np.zeros((6,2)),np.vstack((np.zeros((2,4)),self.z12))))
        # self.z21 = np.hstack((np.zeros((6,2)),np.vstack((np.zeros((2,4)),self.z21))))
        # self.z22 = np.hstack((np.zeros((6,2)),np.vstack((np.zeros((2,4)),self.z22))))
        # 
        # self.z13 = np.vstack((np.zeros((2,1)),self.z13.reshape(4,1)))
        # self.z23 = np.vstack((np.zeros((2,1)),self.z23.reshape(4,1)))
    # def Update(self):
    #   """
    #   For GEBF everything needs to be updated so just use __init__ 
    #   """
