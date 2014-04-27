import MKS
reload(MKS)
import matplotlib.pylab as plt
import numpy as np
import Plot
reload(Plot)

#print MKS.GenDelta(5,False)

MSf1 = MKS.GenDelta(21,False);
MSf2 = MKS.GenDelta(21,True);

#filename1 = "deltaM_modeling/21_1_noah2.dat"
#filename2 = "deltaM_modeling/21_2_noah2.dat"
filename1 = "deltaM_modeling/21deltaM_1surroundedBy2_strain.dat"
filename2 = "deltaM_modeling/21deltaM_2surroundedBy1_strain.dat"



stress1 = MKS.ABstrains(filename1)
stress2 = MKS.ABstrains(filename2)

#slice = stress2[:,:,8]
#MKS.ShowSlice(slice, slice.min(), slice.max())

Macro = .02

coeff = MKS.GenC(MSf1, MSf2, filename1, filename2, Macro)

Plot.plotCoeff(coeff)
Plot.compareCoeff(coeff)

'''

stressCalc1 = MKS.NewResponse(coeff, .02, MSf1)
stressCalc1 = np.real_if_close(stressCalc1)

stressCalc2 = MKS.NewResponse(coeff, .02, MSf2)
stressCalc2 = np.real_if_close(stressCalc2)
#print stressCalc[8,8,8]
#print stress1[8,8,8]

#plt.subplot(411)
plt.imshow(stressCalc1[11,:,:])
plt.title('MKS')
plt.colorbar()
plt.show()
#plt.subplot(412)
plt.imshow(stress1[11,:,:])
plt.title('FEM')
plt.colorbar()
plt.show()
#plt.subplot(413)
plt.imshow(stress1[11,:,:]-stressCalc1[11,:,:])
plt.title('Error')
plt.colorbar()
plt.show()

#plt.subplot(421)
plt.imshow(stressCalc2[11,:,:])
plt.title('MKS')
plt.colorbar()
plt.show()
#plt.subplot(422)
plt.imshow(stress2[11,:,:])
plt.title('FEM')
plt.colorbar()
plt.show()
#plt.subplot(423)
plt.imshow(stress2[11,:,:]-stressCalc2[11,:,:])
plt.title('Error')
plt.colorbar()
plt.show()
'''