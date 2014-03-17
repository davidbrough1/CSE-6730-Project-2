import scipy as sp
import numpy as np
from MS import MS


'''
David Brough
This function takes an Eigen microstructure and returns a function that is
microstructure function in the fifith dimension of the array.
(The first four dimension will be used for space and time.)
The microstructure function provides the percentage of a voxel that belongs
to each state space.

For example

Eigen Microstructure
        [[1,2,1]
 MAT =   [1,2,1]
         [1,2,2]]

Microstructure Function
        [[[[[1, 0]]]    <- MAT[0,0] = 1
          [[[0, 1]]]    <- MAT[0,1] = 2
          [[[1, 0]]]]   <- MAT[0,2] = 1
         [[[[1, 0]]]    <- MAT[1,0] = 1
MicroSF = [[[0, 1]]]    <- MAT[1,1] = 2
          [[[1, 0]]]]   <- MAT[1,2] = 1
         [[[[1, 0]]]    <- MAT[2,0] = 1
          [[[0, 1]]]    <- MAT[2,1] = 2
          [[[0, 1]]]]]  <- MAT[2,2] = 2
'''


def MSf(MAT):
    Dim = MAT.shape
    l = len(Dim)
    phases = sp.unique(MAT)
    numPhases = len(phases)
    if l == 2:
        MicroSF = np.zeros([Dim[0],Dim[1],1,1,numPhases])
        for ii in range(numPhases):
            Mask = np.ma.masked_equal(MAT, phases[ii])
            MicroSF[:,:,0,0,ii] = MicroSF[:,:,0,0,ii]*Mask
            MicroSF[:,:,0,0,ii] = MicroSF[:,:,0,0,ii]/phases[ii]
    elif l == 3:
        MicroSF = np.zeros([Dim[0],Dim[1],Dim[2],1,numPhases])
        for ii in range(numPhases):
            Mask = np.ma.masked_equal(MAT, phases[ii])
            MicroSF[:,:,:,0,ii] = MicroSF[:,:,:,0,ii]*Mask
            MicroSF[:,:,:,0,ii] = MicroSF[:,:,:,0,ii]/phases[ii]
    return MicroSF


#TO DO
#Add Fourth Dimension for Time
