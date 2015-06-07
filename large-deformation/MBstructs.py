import numpy as np
import sympy as sym
import pickle
import MultiBodyFuncts as MBF

from numpy.linalg import inv

# Body is the class that defines a single body in a multibody system
# Only really used for list management
class Body():
    pass
    # basically used

# Rigid_Body is the class that defines a rigid body
class Rigid_Body():
    def __init__(self, I1, I2, I3, q, r, x, v):
        self.I=np.array(((I1,0,0),(0,I2,0),(0,0,I3)))

        # Make the skew matricies needed 
        self.q = q
        self.r = r
        self.q_skew = MBF.skew(q)
        self.r_skew = MBF.skew(r)

        # # This function is not used in DCA
        # self.wn = np.dot(x, self.Pw)
        # self.W0 = np.transpose(MBF.skew(self.wn))

        #This function is not used in DCA
        self.C0dot=np.dot(self.W0, self.C0)

        #Function to transform the Inertia Matrix with the 
        #direction cosine matrices
        self.I0=np.dot(np.transpose(self.C0),np.dot(self.I,self.C0))
        self.I0T=np.transpose(self.I0)

        #Function to create the Shifter Matrices
        # self.S01=np.zeros((6,6))
        # self.S10=np.zeros((6,6))
        # self.S20=np.zeros((6,6))
        # self.S02=np.zeros((6,6))

        # self.S02[:3,:3]=np.identity(3)
        # self.S02[3:,3:]=np.identity(3)
        # self.S10[:3,:3]=np.identity(3)
        # self.S10[3:,3:]=np.identity(3)
        # self.S20[:3,:3]=np.identity(3)
        # self.S20[3:,3:]=np.identity(3)
        # self.S01[:3,:3]=np.identity(3)
        # self.S01[3:,3:]=np.identity(3)
                        
        self.S02 = MBF.skew(self.r02)
        self.S10 = MBF.skew(self.r10)
        self.S01 = MBF.skew(self.r01)	
        self.S20 = MBF.skew(self.r20)

        # Create the Mass and inverse Mass Matrices	
        self.M=np.zeros((6,6))
        self.M[:3,:3]=self.I0
        self.M[3:,3:]=np.array(((self.m,0,0),(0,self.m,0),(0,0,self.m)))		
        self.Minv=np.linalg.inv(self.M)

    #Function that defines position vectors between any of the
    #Three points on the body (two handles and center of mass)
    def rs(self,r,v):
        self.r10= np.dot(r*v,self.C0)
        self.r01=-1*self.r10
        self.r12= np.dot(self.l*v,self.C0)
        self.r21=-1*self.r12
        self.r20= np.dot((r)*v*-1,self.C0) 
        self.r02=-1*self.r20

    #Function that finds the initial zeta values in each timestep
    def zeta(self):
        self.z11=np.dot(self.S10,np.dot(self.Minv,self.S01))
        self.z12=np.dot(self.S10,np.dot(self.Minv,self.S02))
        self.z13=np.dot(np.dot(self.S10,self.Minv),self.Fa1)
        self.z21=np.dot(self.S20,np.dot(self.Minv,self.S01))
        self.z22=np.dot(self.S20,np.dot(self.Minv,self.S02))
        self.z23=np.dot(self.S20,np.dot(self.Minv,self.Fa2))

            
    #Function to apply the state dependant forces
    #I think this is where my problem is because I am 
    #not sure if this is correct.
    def Forces(self,v):

        self.Fa1=np.zeros((6))
        self.Fa2=np.zeros((6))

        #Force from centripedal motion around the
        #body's first joint
        self.i1= np.cross(self.w,self.r10)
        self.i12=np.cross(self.w,self.r20)
        self.i2=np.cross(self.w,self.i1)
        self.i31=self.m*self.i2
        self.i32=self.m*np.cross(self.w,self.i12)
        
        
        #Total force with gravity included
        self.Fa1[3:]=9.81*self.m*v-self.i31#-self.i5
        self.Fa2[3:]=9.81*self.m*v-self.i32

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
        E, I, A, r, rho, l, = sym.symbols('E I A r rho l')

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
        self.zeta11 = Gamma1.dot(self.gamma11 - self.M12.dot(iM22.dot(self.gamma21)))
        self.zeta12 = Gamma1.dot(self.gamma12 - self.M12.dot(iM22.dot(self.gamma22)))
        self.zeta13 = Gamma1.dot(self.gamma13 - self.M12.dot(iM22.dot(self.gamma23)))

        self.zeta21 = Gamma2.dot(self.gamma21 - self.M21.dot(iM11.dot(self.gamma11)))
        self.zeta22 = Gamma2.dot(self.gamma22 - self.M21.dot(iM11.dot(self.gamma12)))
        self.zeta23 = Gamma2.dot(self.gamma23 - self.M21.dot(iM11.dot(self.gamma13)))

        # Make "full-size" to interface with "standard" DCA for rigid bodies
        self.zeta11 = np.hstack((np.zeros((6,2)),np.vstack((np.zeros((2,4)),self.zeta11))))
        self.zeta12 = np.hstack((np.zeros((6,2)),np.vstack((np.zeros((2,4)),self.zeta12))))
        self.zeta21 = np.hstack((np.zeros((6,2)),np.vstack((np.zeros((2,4)),self.zeta21))))
        self.zeta22 = np.hstack((np.zeros((6,2)),np.vstack((np.zeros((2,4)),self.zeta22))))

        self.zeta13 = np.vstack((np.zeros((2,1)),self.zeta13.reshape(4,1)))
        self.zeta23 = np.vstack((np.zeros((2,1)),self.zeta23.reshape(4,1)))
    # def Update(self):
    #   """
    #   For GEBF everything needs to be updated so just use __init__ 
    #   """
