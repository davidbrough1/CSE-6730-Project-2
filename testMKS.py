import MKS
reload(MKS)

#print MKS.GenDelta(5,False)

MSf1 = MKS.GenDelta(21,False);
MSf2 = MKS.GenDelta(21,True);

stress1 = MKS.ABstrains("deltaM_modeling/21_1_noah2.dat")
stress2 = MKS.ABstrains("deltaM_modeling/21_2_noah2.dat")

#slice = stress2[:,:,8]
#MKS.ShowSlice(slice, slice.min(), slice.max())

Macro = .02

coeff = MKS.GenC(MSf1, MSf2, "deltaM_modeling/21_1_noah2.dat", "deltaM_modeling/21_2_noah2.dat", Macro)

stressCalc = MKS.NewResponse(coeff, .02, MSf1)

print stressCalc[8,8,8]
print stress1[8,8,8]
