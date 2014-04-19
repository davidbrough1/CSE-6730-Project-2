from numpy import *

#MicroSF are the Microstructure functions, ABout are the ABAQUS output file names, Macro is the imposed macro strain in the 11 direction
def GenC(MicroSF_1, MicroSF_2 ABout1, ABout2, Macro)
	#process ABAQUS output files and return 11 strains for each element
	strains_1 = ABstrains(ABout1)
	strains_2 = ABstrains(ABout2)
	
	#calculate DFT space responses and microstructure functions
	response_1_k = np.fft.fftn(strains_1, strains_1.size, [0,1,2])
	micro_1_k = np.fft.fftn(MicroSF_1, MicroSF_1.size, [0,1,2])
	
	response_2_k = np.fft.fftn(strains_2, strains_2.size, [0,1,2])
	micro_2_k = np.fft.fftn(MicroSF_2, MicroSF_2.size, [0,1,2])
	
	#divide responses by inputs to begin calculating coefficients
	response_1_k[:,:,:] = [x/Macro for x in response_1_k] 
	response_2_k[:,:,:] = [x/Macro for x in response_2_k] 
	
	#explicitly invert matrix and solve for coefficients
	dim_len = response_1_k.shape[0]
	coeff = zeros( (dim_len,dim_len,dim_len, 2) )
	det = zeros( (dim_len,dim_len,dim_len) )
	det[:,:,:] = 1/( micro_1_k[:,:,:,0]*micro_2_k[:,:,:,1] - micro_2_k[:,:,:,0]*micro_1_k[:,:,:,1])
	coeff[:,:,:,0] = ( response_1_k*micro_2_k[:,:,:,1] - response_2_k*micro_1_k[:,:,:,1] )/det
	coeff[:,:,:,1] = ( response_2_k*micro_1_k[:,:,:,0] - response_1_k*micro_2_k[:,:,:,0] )/det
	
	coeff = np.fft.ifftn(coeff, coeff.size, [0,1,2])
	
	return coeff