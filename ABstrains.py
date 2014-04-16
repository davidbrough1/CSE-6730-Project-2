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