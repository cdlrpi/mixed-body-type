# Jeremy Laflin
# Modified version of MultiBodyFuncts.py by Mike Hans
# In this version
#   - Rotation matricies are reduced to 2x2 simple rotation about 3rd axis
#   - Removed skew functions
import numpy as np
#This script will hold some functions that will be used in other scripts
def DCM(theta):
    C = np.array([[np.cos(theta), -np.sin(theta)],\
                  [np.cos(theta), -np.sin(theta)]])

def vb2pts(va,awb,rab):
	return (va+np.cross(awb,rab))

def ab2pts(aa,aab,rab,awb):
	return aa+np.cross(aab,rab)+np.cross(awb,np.cross(awb,rab))
def v1pt(vbar,v):
	return vbar+v
def a1pt(abar,a,awb,v):
	return abar+a+2*np.cross(awb,v)
def genvel(v1,v2,awb,rab):
	return v1+v2+np.cross(awb,rab)
def genacc(a1,a2,aab,rab,awb, v2):
	return a1+a2+np.cross(aab,rab)+np.cross(awb,np.cross(awb,rab))+2*np.cross(awb,v2)

def PendEnergy(Y,bodies):
	val1=0
	val2=0
	height1=0
	temp_V=0
	height2=0
	n=int(len(Y[0,:])/2)
	theta=0
	w=0
	r=np.zeros((2))
	V_cm=np.zeros((2))
	h=0
	PE=np.zeros((len(Y[:,0])))
	KEr=np.zeros_like(PE)
	KEt=np.zeros_like(PE)
	energy=np.zeros((len(Y[:,0]),3))
	for i in range (0,len(Y[:,0])):
		val1=0
		val2=0
		height1=0
		temp_V=0
		height2=0
		theta=0
		w=0
		h=0
		for j in range (0,n):
			val1+=Y[i,j]
			val2+=Y[i,j+n]
			theta=val1
			w=val2
			I=bodies[j].m*(bodies[j].l**2)/12
			if j == 0:
				r[:]=np.array((bodies[j].l*np.cos(theta)/2,bodies[j].l*np.sin(theta)/2))
				V_cm[:]=r[:]*w
				temp_V+=w*r*2
				h=(bodies[0].l/2)-r[0]
				height1+=bodies[0].l
				height2+=r[0]*2
			else:

				r[:]=np.array((bodies[j].l*np.cos(theta)/2,bodies[j].l*np.sin(theta)/2))
				V_cm[:]=temp_V+w*r[:]
				temp_V+=w*r*2
				h=height1+(bodies[j].l/2)-r[0]-height2
				height1+=bodies[j].l
				height2+=r[0]*2
			PE[i]+=bodies[j].m*9.81*h
			KEr[i]+=.5*I*w*w
			KEt[i]+=.5*bodies[j].m*np.dot(V_cm,V_cm)
			
	KE=KEr+KEt
	TE=PE+KE
	energy[:,0]=KE
	energy[:,1]=PE
	energy[:,2]=TE
	return energy
				
