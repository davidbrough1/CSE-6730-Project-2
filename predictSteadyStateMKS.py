import MKS
import numpy as np
import matplotlib as plt

import MKS
reload(MKS)
import matplotlib.pylab as plt
import numpy as np

#print MKS.GenDelta(5,False)

MSf1 = MKS.GenDelta(21,False);
MSf2 = MKS.GenDelta(21,True);

#filename1 = "deltaM_modeling/21_1_noah2.dat"
#filename2 = "deltaM_modeling/21_2_noah2.dat"
filename1 = "deltaM_modeling/21deltaM_1surroundedBy2_strain.dat"
filename2 = "deltaM_modeling/21deltaM_2surroundedBy1_strain.dat"
filenamePredict = "outputModeling/41_2phase.dat"



stress1 = MKS.ABstrains(filename1)
stress2 = MKS.ABstrains(filename2)

#slice = stress2[:,:,8]
#MKS.ShowSlice(slice, slice.min(), slice.max())

Macro = .02

coeff = MKS.GenC(MSf1, MSf2, filename1, filename2, Macro)
