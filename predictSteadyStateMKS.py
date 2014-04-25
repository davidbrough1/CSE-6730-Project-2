import MKS
import numpy as np
import matplotlib.pylab as plt
import pickle as pk
from MSf import MSf
import time


#print MKS.GenDelta(5,False)

MSf1 = MKS.GenDelta(21,False);
MSf2 = MKS.GenDelta(21,True);

#filename1 = "deltaM_modeling/21_1_noah2.dat"
#filename2 = "deltaM_modeling/21_2_noah2.dat"
filename1 = "deltaM_modeling/21deltaM_1surroundedBy2_strain.dat"
filename2 = "deltaM_modeling/21deltaM_2surroundedBy1_strain.dat"
filenamePredict = "outputModeling/42_2phase_strain.dat"

f = open("outputModeling/42_2phase_strain.inp.microstructure")

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

'''
plt.imshow(np.real_if_close(np.fft.ifftn(coeff[10,:,:,0])))
plt.colorbar()
plt.show()
'''

Coeff0 = np.real_if_close(np.fft.ifftn(coeff[:,:,:,0]))
Coeff0 = np.roll(Coeff0,10,axis = 0)
Coeff0 = np.roll(Coeff0,10,axis = 1)
Coeff0 = np.roll(Coeff0,10,axis = 2)
coeff41 = MKS.ExpandCoeff(coeff,MSf41.shape[0])
plt.subplot(311)     
plt.imshow(np.real_if_close(Coeff0[10,:,:]))
plt.colorbar()
plt.title('Calibrated Coefficients')
Coeff410 = np.fft.ifftn(coeff41[:,:,:,0])
Coeff411 = np.fft.ifftn(coeff41[:,:,:,1])
plt.subplot(312)
plt.imshow(np.real_if_close(Coeff410[10,0:21,0:21]))
plt.colorbar()
plt.title('Expanded Coefficients in Small Window')
plt.subplot(313)
plt.imshow(np.real_if_close(Coeff0[10,:,:]-Coeff410[10,0:21,0:21]))
plt.colorbar()
plt.title('Difference between 2 above')
plt.show()
                        
startTime = time.time()
strainCalc41 = np.real_if_close(MKS.NewResponse(coeff41, Macro, MSf41))
endTime = time.time()
print endTime-startTime

plt.subplot(311)
p1 = plt.imshow(strainCalc41[10,:,:])
plt.colorbar()
plt.title('MKS')
plt.subplot(312)
p2 = plt.imshow(strain41[10,:,:])
plt.colorbar()
plt.title('FEM')
plt.subplot(313)
p3 = plt.imshow(MS41[10,:,:])
plt.colorbar()
plt.show()
plt.title('Microstructure')
plt.close()




'''
plt.imshow(strainCalc41[10,:,:])
plt.show()
plt.imshow(MS41[10,:,:])
plt.show()
'''
