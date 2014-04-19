import numpy
def GenDelta(side_len,swap_phases)
	MSf_1 = zeros((side_len,side_len,side_len,1))
	MSf_2 = ones((side_len,side_len,side_len,1))
	middle = ceil(side_len/2.0)
	if swap_phases:
		MSf = np.concatenate((MSf_1,MSf_2),axis=3)
		MSf[middle,middle,middle,0] = 1
		MSf[middle,middle,middle,1] = 0
	else:
		MSf = np.concatenate((MSf_2,MSf_1),axis=3)
		MSf[middle,middle,middle,0] = 0
		MSf[middle,middle,middle,1] = 1
	return MSF