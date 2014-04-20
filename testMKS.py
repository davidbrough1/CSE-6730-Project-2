import MKS
reload(MKS)

#print MKS.GenDelta(5,False)

MSf1 = MKS.GenDelta(17,False);
MSf2 = MKS.GenDelta(17,True);

stress1 = MKS.ABstrains("deltaM_modeling/17deltaM_1surroundedBy2.dat")
stress2 = MKS.ABstrains("deltaM_modeling/17deltaM_2surroundedBy1.dat")

#slice = stress2[:,:,8]
#MKS.ShowSlice(slice, slice.min(), slice.max())

Macro = .02

coeff = MKS.GenC(MSf1, MSf2, "deltaM_modeling/17deltaM_1surroundedBy2.dat", "deltaM_modeling/17deltaM_2surroundedBy1.dat", Macro)

stressCalc = MKS.NewResponse(coeff, .02, MSf1)

print stressCalc[0,0,0]
print stress1[0,0,0]