import MKS
import numpy as np
import matplotlib.pylab as plt
import pickle as pk
from MSf import MSf
import time


#print MKS.GenDelta(5,False)

MSf1 = MKS.GenDelta(21,False);
MSf2 = MKS.GenDelta(21,True);

filename1 = "deltaM_modeling/21_1_noah2.dat"
filename2 = "deltaM_modeling/21_2_noah2.dat"
#filename1 = "deltaM_modeling/21deltaM_1surroundedBy2_strain.dat"
#filename2 = "deltaM_modeling/21deltaM_2surroundedBy1_strain.dat"
filenamePredict = "outputModeling/41_2phase.dat"

f = open("outputModeling/41_2phase.inp.microstructure")

MS41 = pk.load(f)

f.close()

strain1 = MKS.ABstrains(filename1)
strain2 = MKS.ABstrains(filename2)
strain41 = MKS.ABstrains(filenamePredict)

MSf41 = MSf(MS41)
MSf41 = MSf41[:,:,:,0,:]

#slice = stress2[:,:,8]
#MKS.ShowSlice(slice, slice.min(), slice.max())

Macro = .02


coeff = MKS.GenC(MSf1, MSf2, filename1, filename2, Macro)
coeff41 = MKS.ExpandCoeff(coeff,41)

startTime = time.time()
strainCalc41 = np.real_if_close(MKS.NewResponse(coeff41, Macro, MSf41))
endTime = time.time()
print endTime-startTime

plt.subplot(311)
p1 = plt.imshow(strainCalc41[10,:,:])
plt.colorbar()
plt.subplot(312)
p2 = plt.imshow(strain41[10,:,:])
plt.colorbar()
plt.subplot(313)
p3 = plt.imshow(MS41[10,:,:])
plt.colorbar()
plt.show()
plt.close()
   
