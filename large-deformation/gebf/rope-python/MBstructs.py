import numpy as np
import sympy as sym
import pickle
import MultiBodyFuncts as MBF

from numpy.linalg import inv

#Body is the class that defines a single body in a multibody system
class Rigid_Body:
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

#Joint is the class that defines a single joint in a multibody system		
class Joint:
    # Function to create the P and D matrices
    def __init__(self,v,x):
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
        Initializes all of the inertial propertias of the element and assigns values to physical and material properties 
        """

        # re-make symbols for proper substitution
        e = sym.Matrix(sym.symarray('e',8))
        # symbolic system parameters 
        E, I, A, rho, l, = sym.symbols('E I A rho l')

        # Load symbolic mass matrix
        M = pickle.load( open( "mass-matrix.dump", "rb" ) ).subs([(A, area), (l, length), (rho, density)])
        
        # load the body and applied force vector (this still has sympolic 'e' values)
        self.beta = pickle.load( open( "force-vector.dump", "rb" ) ).subs([(E, modulus), (A, area), \
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

    def BodyForces(self):
        """
        This function updated the body-and-applied force vectors as the state variables change.
        """
        # re-make symbols for proper substitution
        e = sym.Matrix(sym.symarray('e',8))
        e_sub = [(e, ei) for e, ei in zip(e, self.state)]

        # load the body and applied force vector and make substitutions for physical and material properties.
        beta = self.beta.subs(e_sub)

        # partition beta into lambda13 and lambda23
        self.lambda13 = beta[0:4]
        self.lambda23 = beta[4:8]

class GEBF_Element():
    def __init__(self, area, modulus, inertia, density, length, state):
        """
        Initializes all of the inertial propertias of the element and assigns values to physical and material properties 
        """

        # re-make symbols for proper substitution
        e = sym.Matrix(sym.symarray('e',12))
        # symbolic system parameters 
        E, I, A, rho, l, = sym.symbols('E I A rho l')

        # Load symbolic mass matrix
        M = pickle.load( open( "mass-matrix.dump", "rb" ) ).subs([(A, area), (l, length), (rho, density)])
        
        # load the body and applied force vector (this still has sympolic 'e' values)
        self.beta = pickle.load( open( "force-vector.dump", "rb" ) ).subs([(E, modulus), (A, area), \
                                                                           (I, inertia), (rho, density), (l, length)])

        # Partition mass matrix
        self.M11 = np.array(M[0:6,0:6])
        self.M12 = np.array(M[0:6,6:12])
        self.M21 = np.array(M[6:12,0:6])
        self.M22 = np.array(M[6:12,6:12])

        # For now 
        self.lambda11 = np.eye(6)
        self.lambda12 = np.zeros((6,6))
        self.lambda22 = np.eye(6)
        self.lambda21 = np.zeros((6,6))
        
        # fully numberic (initial) values for body and applied forces
        e_sub = [(e, ei) for e, ei in zip(e, state)]
        beta = self.beta.subs(e_sub)

        # partition beta into lambda13 and lambda23
        self.lambda13 = np.array(beta[0:6])
        self.lambda23 = np.array(beta[6:12])


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

    def BodyForces(self):
        """
        This function updated the body-and-applied force vectors as the state variables change.
        """
        # re-make symbols for proper substitution
        e = sym.Matrix(sym.symarray('e',12))
        e_sub = [(e, ei) for e, ei in zip(e, self.state)]

        # load the body and applied force vector and make substitutions for physical and material properties.
        beta = self.beta.subs(e_sub)

        # partition beta into lambda13 and lambda23
        self.lambda13 = beta[0:6]
        self.lambda23 = beta[6:12]
