# Jeremy Laflin
# Modified version of MultiBodyFuncts.py by Mike Hans
# In this version
#   - Rotation matricies are reduced to 2x2 simple rotation about 3rd axis
import numpy as np
#This script will hold some functions that will be used in other scripts
def DCM(theta):
    C = np.array([[np.cos(theta), -np.sin(theta)],\
                  [np.sin(theta),  np.cos(theta)]])
    return C
