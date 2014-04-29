import numpy as np
import re
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

	currentStep = 1
	keepProcessing = True
	strains = None
      # determine if this is a static or visco analysis
	visco = False
	searchStaticOrVisco = "S T E P       {}     (S T A T I C   A N A L Y S I S|V I S C O   A N A L Y S I S)".format(currentStep)
	m = re.search(searchStaticOrVisco, cur_line)
	while(cur_line != "" and not m):
		cur_line = f.readline()
		m = re.search(searchStaticOrVisco, cur_line)	
	
	if (cur_line == ""):
		return None
      
	if m.group(1) == "V I S C O   A N A L Y S I S":
		visco = True
	
	while (keepProcessing):
		# Read strains for the next step, for static analysis, there should be only one
		stepStrains = ReadStrainsStep(currentStep, f, num_elements, visco)
		if (stepStrains == None):
			keepProcessing = False
		else:
			if strains == None:
				strains = stepStrains
			else:
				strains = np.concatenate((strains, stepStrains), axis=3)
		currentStep = currentStep + 1
	return strains

def ReadStrainsStep(currentStep, f, num_elements, visco = False):
	
	cur_line = f.readline()
	
	# get time increments and time period for step if visco
	timeIncrement = 0
	timePeriod = 0
	numberOfTimeSteps = 1
	if visco:
		while(cur_line != "" and "TIME INCREMENT IS" not in cur_line):
			cur_line = f.readline()
		if cur_line == "":
			return None
		timeIncrement = float((cur_line.split())[3])
		while(cur_line != "" and "TIME PERIOD IS" not in cur_line):
			cur_line = f.readline()
		timePeriod = float((cur_line.split())[3])
		numberOfTimeSteps = timePeriod/timeIncrement
		
	currentTimeStep = 0						

	#create storage of element values, only strain in 11 direction currently
	#numStrains are the numeric values for 11, 22, 33, 12, 13, 23 strains
	#                         respectively  0,  1, 2 , 3,  4,  5
	num_strains = range(0,1)
	dim_length = ceil(num_elements ** (1/3.0))
	if visco:
		strains = zeros( (dim_length,dim_length,dim_length,numberOfTimeSteps,size(num_strains)) )
	else:
		strains = zeros( (dim_length,dim_length,dim_length,size(num_strains)) )		
	
	# read strains for the next timestep, if static only one timestep to read
	if (visco):		
		continueReading = True
		while (continueReading):
			continueReading = ReadStrainsTimestep(currentStep, currentTimeStep, f, num_elements, strains, True)
			currentTimeStep = currentTimeStep + 1
	else:
		continueReading = ReadStrainsTimestep(currentStep, currentTimeStep, f, num_elements, strains)
		if (not continueReading):
			return None
	
	# for single strain case reshape strains to only 3d array, will need to be changed for time rate to be nxnxnxt in this case 
	# Including time for the visco state
	if (visco):
		strains = strains[:, :, :, :, 0]
	else:
		strains = strains[:, :, :, 0]
		
	return strains

def ReadStrainsTimestep(currentStep, timestep, f, num_elements, strains, visco=False):
	searchStaticOrVisco = "S T E P       {}     (S T A T I C   A N A L Y S I S|V I S C O   A N A L Y S I S)".format(currentStep + 1)
	
	num_strains = range(0,1)
	dim_length = ceil(num_elements ** (1/3.0))
	layer_size = dim_length ** 2
	
	cur_line = f.readline()
	
	#find start of element data
	while(cur_line != "" and "THE FOLLOWING TABLE IS PRINTED FOR ALL ELEMENTS "\
		  "WITH TYPE C3D8 AT THE INTEGRATION POINTS" not in cur_line):
		cur_line = f.readline()
		# reach end of file or the next step in the analysis
		if cur_line == "" or re.search(searchStaticOrVisco, cur_line):
			return False
	
	for i in range(0,5):
		cur_line = f.readline()
		
	integration_point_strains = zeros( (8,size(num_strains)) )
	#1st dim = x, 2nd dim = y, 3rd dim = z
	for i in range(0,num_elements):
		y = i%dim_length
		z = (i%layer_size)/dim_length
		x = i/layer_size
		#average strains over the 8 integration points for c3d8
		for k in range(0,8):
			split_str = cur_line.split()
			for j in num_strains:
				integration_point_strains[k,j] = float(split_str[2+j])
			cur_line = f.readline()
		if (visco):
			strains[x,y,z,timestep,:] = integration_point_strains.mean(axis=0)
		else:
			strains[x,y,z,:] = integration_point_strains.mean(axis=0)
	return True
	
	
#MicroSF are the Microstructure functions, ABout are the ABAQUS output file names, Macro is the imposed macro strain in the 11 direction
def GenC(MicroSF_1, MicroSF_2, ABout1, ABout2, Macro):
	#process ABAQUS output files and return 11 strains for each element
	strains_1 = ABstrains(ABout1)
	strains_2 = ABstrains(ABout2)
	strains_1 = strains_1/np.mean(strains_1)
        strains_2 = strains_2/np.mean(strains_2)	
	#calculate DFT space responses and microstructure functions
	response_1_k = np.fft.fftn(strains_1)
	micro_1_k = zeros(MicroSF_1.shape,dtype=complex)
	micro_1_k[:,:,:,0] = np.fft.fftn(MicroSF_1[:,:,:,0])
	micro_1_k[:,:,:,1] = np.fft.fftn(MicroSF_1[:,:,:,1])
	
	response_2_k = np.fft.fftn(strains_2)
	micro_2_k = zeros(MicroSF_2.shape,dtype=complex)
	micro_2_k[:,:,:,0] = np.fft.fftn(MicroSF_2[:,:,:,0])
	micro_2_k[:,:,:,1] = np.fft.fftn(MicroSF_2[:,:,:,1])

	#explicitly invert matrix and solve for coefficients
	dim_len = response_1_k.shape[0]
		
	coeff = zeros( (dim_len,dim_len,dim_len,2) , dtype=complex)

	#attempt to perform least squares solution
	for i in range(dim_len):
		for j in range(dim_len):
			for k in range(dim_len):
	
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
	
				if ((i==0 and j ==0) and (k==0 or k==1)):
					p = independent_columns(MM, .001)

				calred = MM[p,:][:,p]      # Linearly independent columns of MM
				resred = PM[p,0].conj().T  # Linearly independent columns of PM 

				coeff[i,j,k,p] = np.linalg.solve(calred, resred)   

	
	#coeff[:,:,:] = response_1_k/micro_1_k[:,:,:,0]
	
	return coeff

#takes in conj of DFT coefficients, output new matrix with conj of DFT coefficients such that it is ready to use in the DFT space

#Takes the influence coefficients and zero padds them to fit a larger microstructure
def ExpandCoeff(coeff, new_side_len):
        coeff[:,:,:,0] = np.real_if_close(np.fft.ifftn(coeff[:,:,:,0]))
        coeff[:,:,:,0] = np.fft.fftshift(coeff[:,:,:,0])
        new_coeff = np.ndarray(shape = (new_side_len,new_side_len,new_side_len,coeff.shape[3]),dtype = 'complex128')
        new_coeff[new_side_len/2-10:new_side_len/2+11,new_side_len/2-10:new_side_len/2+11,\
                  new_side_len/2-10:new_side_len/2+11,0]=coeff[:,:,:,0]
        new_coeff[:,:,:,0] = np.fft.ifftshift(new_coeff[:,:,:,0])
        new_coeff[:,:,:,0] = np.fft.fftn(new_coeff[:,:,:,0])
        new_coeff[0,0,0,1] = new_coeff[0,0,0,0]
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
	
	return np.real_if_close(response)


#graphically displays the grayscale picture of a 2D array slice
#min and max vals allow for consistent scaling between two output arrays
def ShowSlice(vals, min_val, max_val):
	import matplotlib.pyplot as plt

	# tell imshow about color map so that only set colors are used
	img = plt.imshow(vals,cmap = plt.get_cmap('gray'),vmin=min_val, vmax=max_val)

	plt.show()


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
