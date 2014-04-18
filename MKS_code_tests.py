import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
from MS import MS
from MSf import MSf


#Tests for MS

Mat2D = sp.rand(93,93)
MicroS2D = MS(Mat2D,3,1,20,10,1)
#plt.imshow(MicroS2D)
#plt.show()
Mat3D = sp.rand(95,95,95)
MicroS3D = MS(Mat3D,3,1,20,10,20)
#plt.imshow(MicroS3D[1,:,:])
#plt.show()


#Tests for MSf

MAT2D = MS(sp.rand(3,3),2,1,1,1,1)
#print MAT2D
MSf2D = MSf(MAT2D )
#print MSf2D
MSf3D = MSf(MS(sp.rand(3,3,3),3,1,1,1,1))
#print MSf3D
