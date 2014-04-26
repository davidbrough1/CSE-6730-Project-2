import MKS as MKS
import matplotlib.pylab as plt
import numpy as np

MSf211 = MKS.GenDelta(21,False);
MSf212 = MKS.GenDelta(21,True);
MSf411 = MKS.GenDelta(41,False);
MSf412 = MKS.GenDelta(41,True);

filename211 = "deltaM_modeling/21deltaM_1surroundedBy2.dat"
filename212 = "deltaM_modeling/21deltaM_2surroundedBy1.dat"
filename411 = "deltaM_modeling/41deltaM_1surroundedBy2.dat"
filename412 = "deltaM_modeling/41deltaM_2surroundedBy1.dat"

strain211 = MKS.ABstrains(filename211)
strain212 = MKS.ABstrains(filename212)
strain411 = MKS.ABstrains(filename411)
strain411 = MKS.ABstrains(filename412)

Macro = .02

fcoeff21 = MKS.GenC(MSf211, MSf212, filename211, filename212, Macro)
fcoeff41 = MKS.GenC(MSf411, MSf412, filename411, filename412, Macro)

fexpandcoeff21 = fcoeff21
fcoeff21exp = MKS.ExpandCoeff(fexpandcoeff21,41)

coeff21 = np.ndarray(shape = fcoeff21exp.shape,dtype = 'float')

coeff41 = np.ndarray(shape = fcoeff41.shape,dtype = 'float')

coeff21[:,:,:,0] = np.real_if_close(np.fft.ifftn(fcoeff21exp[:,:,:,0]))
coeff21[:,:,:,1] = np.real_if_close(np.fft.ifftn(fcoeff21exp[:,:,:,1]))

coeff41[:,:,:,0] = np.real_if_close(np.fft.ifftn(fcoeff41[:,:,:,0]))
coeff41[:,:,:,1] = np.real_if_close(np.fft.ifftn(fcoeff41[:,:,:,1]))

coeff21[:,:,:,0] = np.roll(coeff21[:,:,:,0],10,axis = 0)
coeff21[:,:,:,0] = np.roll(coeff21[:,:,:,0],10,axis = 1)
coeff21[:,:,:,0] = np.roll(coeff21[:,:,:,0],10,axis = 2)
coeff21[:,:,:,1] = np.roll(coeff21[:,:,:,1],10,axis = 0)
coeff21[:,:,:,1] = np.roll(coeff21[:,:,:,1],10,axis = 1)
coeff21[:,:,:,1] = np.roll(coeff21[:,:,:,1],10,axis = 2)


coeff41[:,:,:,0] = np.roll(coeff41[:,:,:,0],20,axis = 0)
coeff41[:,:,:,0] = np.roll(coeff41[:,:,:,0],20,axis = 1)
coeff41[:,:,:,0] = np.roll(coeff41[:,:,:,0],20,axis = 2)
coeff41[:,:,:,1] = np.roll(coeff41[:,:,:,1],20,axis = 0)
coeff41[:,:,:,1] = np.roll(coeff41[:,:,:,1],20,axis = 1)
coeff41[:,:,:,1] = np.roll(coeff41[:,:,:,1],20,axis = 2)


plt.subplot(311)
plt.imshow(coeff21[20,:,:,0])
plt.colorbar()
plt.title('21 H = 0')
plt.subplot(312)
plt.imshow(coeff41[20,:,:,0])
plt.colorbar()
plt.title('41 H = 0')
plt.subplot(313)
plt.imshow(coeff41[20,:,:,0]-coeff21[20,:,:,0])
plt.colorbar()
plt.title('Difference')
plt.show()

'''
plt.subplot(311)
plt.imshow(coeff21[20,:,:,1])
plt.colorbar()
plt.title('21 H = 1')
plt.subplot(312)
plt.imshow(coeff41[20,:,:,1])
plt.colorbar()
plt.title('41 H = 1')
plt.subplot(313)
plt.imshow(coeff41[20,:,:,1]-coeff21[20,:,:,1])
plt.colorbar()
plt.title('Difference')
plt.show()
'''

print ' '
print '21 H = 1'
print coeff21[0,0,0,1]
print coeff21[0,0,0,1]*41**3
print ' '
print '41 H = 1'
print coeff41[0,0,0,1]
print coeff41[0,0,0,1]*41**3

coef21 = np.ndarray(shape = fcoeff21.shape,dtype = 'float')

coef21[:,:,:,0] = np.real_if_close(np.fft.ifftn(fcoeff21[:,:,:,0]))
coef21[:,:,:,1] = np.real_if_close(np.fft.ifftn(fcoeff21[:,:,:,1]))
coef21[:,:,:,0] = np.roll(coef21[:,:,:,0],10,axis = 0)
coef21[:,:,:,0] = np.roll(coef21[:,:,:,0],10,axis = 1)
coef21[:,:,:,0] = np.roll(coef21[:,:,:,0],10,axis = 2)

plt.subplot(311)
plt.imshow(coef21[10,:,:,0])
plt.colorbar()
plt.title('21 H = 0')
plt.subplot(312)
plt.imshow(coeff41[20,10:31,10:31,0])
plt.colorbar()
plt.title('41 H = 0')
plt.subplot(313)
plt.imshow(coeff41[20,10:31,10:31,0]-coef21[10,:,:,0])
plt.colorbar()
plt.title('Difference')
plt.show()

'''
plt.subplot(311)
plt.imshow(coef21[10,:,:,1])
plt.colorbar()
plt.title('21 H = 1')
plt.subplot(312)
plt.imshow(coeff41[20,10:31,10:31,1])
plt.colorbar()
plt.title('41 H = 1')
plt.subplot(313)
plt.imshow(coeff41[20,10:31,10:31,1]-coef21[20,:,:,1])
plt.colorbar()
plt.title('Difference')
plt.show()
'''

print ' '
print '21 H = 1'
print coef21[0,0,0,1]
print coef21[0,0,0,1]*21**3
print ' '
print '41 H = 1'
print coeff41[0,0,0,1]
print coeff41[0,0,0,1]*41**3



