import numpy as np
from numpy import *
import matplotlib.pylab as plt

#generate a delta microstructure function, returns a 4d array, 3d spatial 1d microstructure
def GenDelta(side_len,swap_phases):
	MSf_1 = zeros((side_len,side_len,side_len,1))
	MSf_2 = ones((side_len,side_len,side_len,1))
	middle = side_len/2
	MSf = np.concatenate((MSf_1,MSf_2),axis=3)
	MSf[middle,middle,middle,0] = 1
	MSf[middle,middle,middle,1] = 0
	if swap_phases:
		MSf = MSf[:,:,:,::-1]
	return MSf

#reads in the Abaqus output file and returns the desired strains for all elements in a 3d array
def ABstrains(ABout):
	#'17_2phase.dat' 
	f = open(ABout, 'r')
	cur_line = f.readline()
	while(cur_line != "" and "NUMBER OF ELEMENTS IS" not in cur_line):
		cur_line = f.readline()
	#get number of elements to loop through
	num_elements = int((cur_line.split())[4])


	#find start of element data
	while(cur_line != "" and "THE FOLLOWING TABLE IS PRINTED FOR ALL ELEMENTS "\
		  "WITH TYPE C3D8 AT THE INTEGRATION POINTS" not in cur_line):
		cur_line = f.readline()
	for i in range(0,5):
		cur_line = f.readline()
	
	#create storage of element values, only strain in 11 direction currently
	#numStrains are the numeric values for 11, 22, 33, 12, 13, 23 strains
	#                         respectively  0,  1, 2 , 3,  4,  5
	num_strains = range(0,1)
	dim_length = ceil(num_elements ** (1/3.0))
	layer_size = dim_length ** 2
	strains = zeros( (dim_length,dim_length,dim_length,size(num_strains)) )
	integration_point_strains = zeros( (8,size(num_strains)) )
	#1st dim = x, 2nd dim = y, 3rd dim = z
	for i in range(0,num_elements):
		z = i%dim_length
		y = (i%layer_size)/dim_length
		x = i/layer_size
		#average strains over the 8 integration points for c3d8
		for k in range(0,8):
			split_str = cur_line.split()
			for j in num_strains:
				integration_point_strains[k,j] = float(split_str[2+j])
			cur_line = f.readline()
		strains[x,y,z,:] = integration_point_strains.mean(axis=0)

	#for single strain case reshape strains to only 3d array, will need to be changed for time rate to be nxnxnxt in this case 
	strains = strains[:,:,:,0]
	return strains
	
#MicroSF are the Microstructure functions, ABout are the ABAQUS output file names, Macro is the imposed macro strain in the 11 direction
def GenC(MicroSF_1, MicroSF_2, ABout1, ABout2, Macro):
	#process ABAQUS output files and return 11 strains for each element
	strains_1 = ABstrains(ABout1)
	strains_2 = ABstrains(ABout2)
	'''
	plt.imshow(strains_2[11,:,:])
	plt.colorbar()
	plt.show()
	plt.imshow(strains_2[:,11,:])
	plt.colorbar()
	plt.show()
	plt.imshow(strains_2[:,:,11])
	plt.colorbar()
	plt.show()
	
	plt.imshow(MicroSF_1[10,:,:,1])
	plt.colorbar()
	plt.show()
	plt.imshow(MicroSF_1[:,10,:,1])
	plt.colorbar()
	plt.show()
	plt.imshow(MicroSF_1[:,:,10,1])
	plt.colorbar()
	plt.show()
	plt.imshow(MicroSF_2[10,:,:,1])
	plt.colorbar()
	plt.show()
	plt.imshow(MicroSF_2[:,10,:,1])
	plt.colorbar()
	plt.show()
	plt.imshow(MicroSF_2[:,:,10,1])
	plt.colorbar()
	plt.show()
	'''
	
	#calculate DFT space responses and microstructure functions
	response_1_k = np.fft.fftn(strains_1)
	micro_1_k = zeros(MicroSF_1.shape,dtype=complex)
	micro_1_k[:,:,:,0] = np.fft.fftn(MicroSF_1[:,:,:,0])
	micro_1_k[:,:,:,1] = np.fft.fftn(MicroSF_1[:,:,:,1])
	#print (micro_1_k[0,0,0,0]+micro_1_k[0,0,0,1])
	response_2_k = np.fft.fftn(strains_2)
	micro_2_k = zeros(MicroSF_2.shape,dtype=complex)
	micro_2_k[:,:,:,0] = np.fft.fftn(MicroSF_2[:,:,:,0])
	micro_2_k[:,:,:,1] = np.fft.fftn(MicroSF_2[:,:,:,1])
	
	'''
        

        plt.imshow(MicroSF_2[10,:,:,0])
        plt.colorbar()
        plt.show()
	plt.imshow(MicroSF_2[10,:,:,1])
        plt.colorbar()
        plt.show()
	
	M1S0 = np.real_if_close(np.fft.ifftn(micro_1_k[:,:,:,0]))
	M1S1 = np.real_if_close(np.fft.ifftn(micro_1_k[:,:,:,1]))
	M1S0 = np.roll(M1S0,10, axis = 0)
        M1S0 = np.roll(M1S0,10, axis = 1)
        M1S0 = np.roll(M1S0,10, axis = 2)
        plt.imshow(M1S0[10,:,:])
        plt.colorbar()
        plt.show()
        M1S1 = np.roll(M1S1,10, axis = 0)
        M1S1 = np.roll(M1S1,10, axis = 1)
        M1S1 = np.roll(M1S1,10, axis = 2)
        plt.imshow(M1S1[10,:,:])
        plt.colorbar()
        plt.show()
        
	M2S0 = np.real_if_close(np.fft.ifftn(micro_2_k[:,:,:,0]))
	M2S1 = np.real_if_close(np.fft.ifftn(micro_2_k[:,:,:,1]))
	M2S0 = np.roll(M2S0,10, axis = 0)
        M2S0 = np.roll(M2S0,10, axis = 1)
        M2S0 = np.roll(M2S0,10, axis = 2)
        plt.imshow(M2S0[10,:,:])
        plt.colorbar()
        plt.show()
        M2S1 = np.roll(M2S1,10, axis = 0)
        M2S1 = np.roll(M2S1,10, axis = 1)
        M2S1 = np.roll(M2S1,10, axis = 2)
        plt.imshow(M2S1[10,:,:])
        plt.colorbar()
        plt.show()
        
	
	print M2S0[10:13,10:13,10:13]
	print M2S1[10:12,10:12,10:12]
	print (micro_2_k[0,0,0,0]+micro_2_k[0,0,0,1])
        '''
	
	#divide responses by inputs to begin calculating coefficients
	#response_1_k[:,:,:] = [x/Macro for x in response_1_k] 
	#response_2_k[:,:,:] = [x/Macro for x in response_2_k] 
	
	#explicitly invert matrix and solve for coefficients
	dim_len = response_1_k.shape[0]
	#coeff = zeros( (dim_len,dim_len,dim_len, 2) , dtype=complex)
	#det = zeros( (dim_len,dim_len,dim_len), dtype=complex )
	#det[:,:,:] = ( micro_1_k[:,:,:,0]*micro_2_k[:,:,:,1] - micro_2_k[:,:,:,0]*micro_1_k[:,:,:,1])
	#print det
	#coeff[:,:,:,0] = ( response_1_k*micro_2_k[:,:,:,1] - response_2_k*micro_1_k[:,:,:,1] )/det
	#coeff[:,:,:,1] = ( response_2_k*micro_1_k[:,:,:,0] - response_1_k*micro_2_k[:,:,:,0] )/det
	
	#attempt to track error
	#errorLocs = np.where(det == 0)
	#print errorLocs[0][2]
	#print errorLocs[1][2]
	#print errorLocs[2][2]
	
	#for blah in range(17):
	#	print micro_1_k[0,0,blah,0]
	#	print micro_1_k[0,0,blah,1]
	micro_1_k[:,:,:,1] = [x*1.5 for x in micro_1_k[:,:,:,1]]
	micro_2_k[:,:,:,1] = [x*1.5 for x in micro_2_k[:,:,:,1]]
	coeff = zeros( (dim_len,dim_len,dim_len,2) , dtype=complex)

	#attempt to perform least squares solution
	for i in range(dim_len):
		for j in range(dim_len):
			for k in range(dim_len):
# 				if(i==0 and k==0 and j==0):
# 					break
# 				x1 = micro_1_k[i,j,k,0] 
# 				x2 = micro_2_k[i,j,k,0]
# 				y1 = response_1_k[i,j,k]
# 				y2 = response_2_k[i,j,k]
# 				y_vec = zeros((2,1), dtype=complex)
# 				y_vec[0] = y1
# 				y_vec[1] = y2
# 				trans_inv = zeros((2,2), dtype=complex)
# 				trans_inv[0,0] = x1**2+x2**2
# 				trans_inv[0,1] = -x1-x2
# 				trans_inv[1,0] = -x1-x2
# 				trans_inv[1,1] = 2
# 				det = (x1**2-2*x1*x2+x2**2)
# 				trans_inv[:,:] = [x/det for x in trans_inv] 
# 				x_trans = zeros((2,2), dtype=complex)
# 				x_trans[0,0] = 1
# 				x_trans[0,1] = 1
# 				x_trans[1,0] = x1
# 				x_trans[1,1] = x2
# 				crap = np.dot(trans_inv,x_trans)
# 				temp = np.dot(crap,y_vec)
# 				arg = temp[1]
# 				coeff[i,j,k] = arg[0]
	
				#u,v,w = i,j,k
				#The following 2 quantities are used in the normal equation during regression.
				MM = np.zeros((2,2)) + 0j*np.zeros((2,2))# Microstucture matrix (M' M*) 
				PM = np.zeros((2,1)) + 0j*np.zeros((2,1))# Property matrix (P M*) 

				mSQc = np.conjugate(micro_1_k[i,j,k,:])  # Conjugate of FFT of Microstructure     
				mSQt = np.mat(micro_1_k[i,j,k,:]).T      # Transpose of FFT of Microstructure

				MM = MM + np.outer(mSQt,mSQc)        # Calculate MM
				PM = PM + (response_1_k[i,j,k] * mSQc) # Calculate PM
				
				mSQc = np.conjugate(micro_2_k[i,j,k,:])  # Conjugate of FFT of Microstructure     
				mSQt = np.mat(micro_2_k[i,j,k,:]).T      # Transpose of FFT of Microstructure

				MM = MM + np.outer(mSQt,mSQc)        # Calculate MM
				PM = PM + (response_2_k[i,j,k] * mSQc) # Calculate PM
	
				#if ((i==0 and j ==0) and (k==0 or k==1)):
                                p = independent_columns(MM, .001)

				calred = MM[p,:][:,p]      # Linearly independent columns of MM
				resred = PM[p,0].conj().T  # Linearly independent columns of PM 

				coeff[i,j,k,p] = np.linalg.solve(calred, resred)   

	
	#coeff[:,:,:] = response_1_k/micro_1_k[:,:,:,0]
	
	return coeff

#takes in conj of DFT coefficients, output new matrix with conj of DFT coefficients such that it is ready to use in the DFT space
def ExpandCoeff(coeff, new_side_len):
	old_side_len = coeff.shape[0]
	if(old_side_len > new_side_len):
		return coeff
		
	#perform complex conjugate
	coeff = np.conj(coeff)
	#coefficients to the spatial format instead of DFT
	coeff[:,:,:,0] = np.fft.ifftn(coeff[:,:,:,0])
	coeff[:,:,:,1] = np.fft.ifftn(coeff[:,:,:,1])
	
	coeff_centered = np.zeros(coeff.shape,dtype = 'float')
	for ii in range(coeff.shape[3]):
		coeff_centered[:,:,:,ii] = np.real_if_close(np.fft.ifftn(coeff[:,:,:,ii]))
    	coeff_centered[:,:,:,ii] = np.roll(coeff[:,:,:,ii], 10, axis = 0)
    	coeff_centered[:,:,:,ii] = np.roll(coeff[:,:,:,ii], 10, axis = 1)
    	coeff_centered[:,:,:,ii] = np.roll(coeff[:,:,:,ii], 10, axis = 2)
	
	#both old side length and new should be odd so centering is easy
	new_coeff = zeros((new_side_len,new_side_len,new_side_len,2))
	offset = (new_side_len-old_side_len)/2
	last_position = offset+old_side_len
	new_coeff[offset:last_position,offset:last_position,offset:last_position,:] = coeff_centered
	
	#convert coeff to DFT space
	new_coeff[:,:,:,0] = np.fft.fftn(new_coeff[:,:,:,0])
	new_coeff[:,:,:,1] = np.fft.fftn(new_coeff[:,:,:,1])
	#perform complex conjugate
	new_coeff = np.conj(new_coeff)
	
	return new_coeff

#solves for new response given inputs
#macro is the imposed macro strain, coeff are the conj DFT coefficients, MSf is the microstructure function
#returns the spatial (not DFT) response
def NewResponse(coeff, macro, MSf):
	response = zeros(coeff.shape[0:3])
	
	MSf_DFT = zeros(MSf.shape,dtype=complex)
	MSf_DFT[:,:,:,0] = np.fft.fftn(MSf[:,:,:,0])
	MSf_DFT[:,:,:,1] = np.fft.fftn(MSf[:,:,:,1])
	
	response = lin_sum = np.sum(np.conjugate(coeff) * MSf_DFT[:,:,:,:], 3)
	response = np.fft.ifftn(response)
	return response


#graphically displays the grayscale picture of a 2D array slice
#min and max vals allow for consistent scaling between two output arrays
def ShowSlice(vals, min_val, max_val):
	import matplotlib.pyplot as plt

	# tell imshow about color map so that only set colors are used
	img = plt.imshow(vals,cmap = plt.get_cmap('gray'),vmin=min_val, vmax=max_val)

	plt.show()

"""
Noah Paulson, 3/27/2014
"""
def independent_columns(A, tol = 1e-05):
    """
    This function is from: http://stackoverflow.com/q/13312498

    Returns an array composed of independent columns of A.

    Note that answer may not be unique; this function returns one of many
    possible answers.
    """
    Q, R = np.linalg.qr(A)
    independent = np.where(np.abs(R.diagonal()) > tol)[0]
    #return A[:, independent]
    return independent
