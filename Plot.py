import matplotlib.pylab as plt
import numpy as np
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm

#plots coefficients in the spatial domain after taking in frequency domain coefficients
def plotCoeff(fcoeff):
	
	coeff = spatialCoeff(fcoeff)
	#plot a surface graph of the influence coefficients for the middle slice
	fig = plt.figure()
	ax = fig.gca(projection='3d')

	dim_len = coeff.shape[0]

	Z = coeff[dim_len/2,:,:,0]
	X = np.zeros((dim_len,dim_len))
	Y = np.zeros((dim_len,dim_len))
	for i in range(dim_len):
		for j in range(dim_len):
			X[i,j] = i
			Y[i,j] = j
	ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=0.3,cmap=cm.coolwarm)
	cset = ax.contour(X, Y, Z, zdir='z', offset=Z.min(), cmap=cm.coolwarm)
	cset = ax.contour(X, Y, Z, zdir='x', offset=0, cmap=cm.coolwarm)
	cset = ax.contour(X, Y, Z, zdir='y', offset=0, cmap=cm.coolwarm)

	ax.set_xlabel('X')
	ax.set_xlim(0, dim_len)
	ax.set_ylabel('Y')
	ax.set_ylim(0, dim_len)
	ax.set_zlabel('Z')
	ax.set_zlim(Z.min(), Z.max())

	plt.show()
	
def compareCoeff(fcoeff):
	coeff = spatialCoeff(fcoeff)
	coeff = coeff[:,:,:,0]
	dim_len = coeff.shape[0]
	
	print 'Coefficient min, max:'
	fcoeff = fcoeff[:,:,:,0]
	abs = np.absolute(fcoeff)
	print abs.min()
	print abs.max()
	
	#compare symmetry of coefficients across all axis
	dim1 = coeff[:,dim_len/2,dim_len/2]
	dim2 = coeff[dim_len/2,:,dim_len/2]
	dim3 = coeff[dim_len/2,dim_len/2,:]
	print 'dimension coefficients'
	print dim1
	print dim2
	print dim3
		
	
	
	
def spatialCoeff(fcoeff):
	coeff = np.zeros(fcoeff.shape,dtype = 'float')
	dim_len = coeff.shape[0]
	for ii in range(coeff.shape[3]):
		coeff[:,:,:,ii] = np.real_if_close(np.fft.ifftn(fcoeff[:,:,:,ii]))
		coeff[:,:,:,ii] = np.roll(coeff[:,:,:,ii], dim_len/2, axis = 0)
		coeff[:,:,:,ii] = np.roll(coeff[:,:,:,ii], dim_len/2, axis = 1)
		coeff[:,:,:,ii] = np.roll(coeff[:,:,:,ii], dim_len/2, axis = 2)
	return coeff