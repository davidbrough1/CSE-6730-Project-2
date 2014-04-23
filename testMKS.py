import MKS
reload(MKS)
import MSf as MSf
import matplotlib.pylab as plt
import numpy as np

#print MKS.GenDelta(5,False)

MSf1 = MKS.GenDelta(21,False);
MSf2 = MKS.GenDelta(21,True);
'''
print 'MSf1,0'
print MSf1[9:12,9:12,9:12,0]
print 'MSf1,1'
print MSf1[9:12,9:12,9:12,1]
print ' '
print 'MSf2,0'
print MSf2[9:12,9:12,9:12,0]
print 'MSf2,1'
print MSf2[9:12,9:12,9:12,1]
'''

stress1 = MKS.ABstrains("deltaM_modeling/21deltaM_1surroundedBy2.dat")
stress2 = MKS.ABstrains("deltaM_modeling/21deltaM_2surroundedBy1.dat")
stress41 = MKS.ABstrains("outputModeling/41_2phase.dat")

"""
Noah's Code Test
"""

stressNoah1 = stress1.reshape((1,1,21**3)
stressNoah2 = stress2.reshape((1,1,21**3)




#slice = stress2[:,:,8]
#MKS.ShowSlice(slice, slice.min(), slice.max())

Macro = .02

fcoeff = MKS.GenC(MSf1, MSf2, "deltaM_modeling/21deltaM_1surroundedBy2.dat", "deltaM_modeling/21deltaM_2surroundedBy1.dat", Macro)


coeff = np.zeros(fcoeff.shape,dtype = 'float')
for ii in range(coeff.shape[3]):
    coeff[:,:,:,ii] = np.real_if_close(np.fft.ifftn(fcoeff[:,:,:,ii]))
    coeff[:,:,:,ii] = np.roll(coeff[:,:,:,ii], 10, axis = 0)
    coeff[:,:,:,ii] = np.roll(coeff[:,:,:,ii], 10, axis = 1)
    coeff[:,:,:,ii] = np.roll(coeff[:,:,:,ii], 10, axis = 2)
'''
for ii in range(2):
    Coef11 = coeff[11,:,:,ii]
    C1 = plt.imshow(Coef11)
    plt.colorbar()
    plt.title('Coef 11 H='+str(ii+1))
    plt.show()
    Coef22 = coeff[:,11,:,ii]
    plt.imshow(Coef22)
    plt.colorbar()
    plt.title('Coef 22 H='+str(ii+1))
    plt.show()
    Coef33 = coeff[:,:,11,ii]
    plt.imshow(Coef33)
    plt.colorbar()
    plt.title('Coef 33 H='+str(ii+1))
    plt.show()
'''
MSf41 = MSf.MSf(stress41)
stressCalc = MKS.NewResponse(fcoeff, .02, MSf2)
stressCalc = np.real_if_close(stressCalc)

plt.subplot(311)
plt.imshow(stressCalc[8,:,:])
plt.colorbar()
#plt.show()
plt.subplot(321)
plt.imshow(stress1[8,:,:])
plt.colorbar()
#plt.show()
aveS = np.average(stress1)
plt.subplot(331)
#plt.imshow((stress1[8,:,:]-stressCalc[8,:,:])/(aveS*stress1.shape[0]**3))
plt.imshow((stress1[8,:,:]-stressCalc[8,:,:]))
plt.colorbar()
plt.show()

