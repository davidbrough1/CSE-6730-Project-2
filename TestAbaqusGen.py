import scipy as sp
import numpy as np
import math
from MS import MS
from MSf import MSf
from AbaqusGen import generateAbaqusInp
from AbaqusGen import getEleNumber
from AbaqusGen import getIndices

ms = MS(sp.rand(31,31,31),2,1,1,1,1)

generateAbaqusInp('23_2phase.inp', MS(sp.rand(23,23,23),2,1,1,1,1))
generateAbaqusInp('21_2phase.inp', MS(sp.rand(21,21,21),2,1,1,1,1))
generateAbaqusInp('17_2phase.inp', MS(sp.rand(17,17,17),2,1,1,1,1))

#generateAbaqusInp('81_2phase.inp', MS(sp.rand(81,81,81),2,1,1,1,1))

generateAbaqusInp('21_2phase_strain.inp', MS(sp.rand(21,21,21),2,1,1,1,1))
generateAbaqusInp('41_2phase_strain.inp', MS(sp.rand(41,41,41),2,1,1,1,1))
generateAbaqusInp('42_2phase_strain.inp', MS(sp.rand(42,42,42),2,1,1,1,1))
generateAbaqusInp('63_2phase.inp', MS(sp.rand(63,63,63),2,1,1,1,1))

#generateAbaqusInp('101_2phase.inp', MS(sp.rand(101,101,101),2,1,1,1,1))

def generateDeltaM(sideLen, viscoelastic=False):
    deltaM_1surroundedBy2 = sp.ones((sideLen, sideLen, sideLen))*2
    middle = sideLen//2
    deltaM_1surroundedBy2[middle, middle, middle] = 1
    if (viscoelastic):
        generateAbaqusInp('{}deltaM_1surroundedBy2_viscoelastic.inp'.format(sideLen), deltaM_1surroundedBy2, True)
    else:
        generateAbaqusInp('{}deltaM_1surroundedBy2.inp'.format(sideLen), deltaM_1surroundedBy2,True)
    deltaM_2surroundedBy1 = sp.ones((sideLen, sideLen, sideLen))
    deltaM_2surroundedBy1[middle, middle, middle] = 2
    if (viscoelastic):
        generateAbaqusInp('{}deltaM_2surroundedBy1_viscoelastic.inp'.format(sideLen), deltaM_2surroundedBy1, True)
    else:
        generateAbaqusInp('{}deltaM_2surroundedBy1.inp'.format(sideLen), deltaM_1surroundedBy2)
    
generateDeltaM(17)
generateDeltaM(21)
generateDeltaM(23)

generateDeltaM(41)

generateDeltaM(15,True)
generateDeltaM(17,True)
generateDeltaM(21,True)
generateDeltaM(23,True)

generateAbaqusInp('15_2phase_viscoelastic.inp', MS(sp.rand(15,15,15),2,1,1,1,1), True)
generateAbaqusInp('21_2phase_viscoelastic.inp', MS(sp.rand(21,21,21),2,1,1,1,1), True)
generateAbaqusInp('31_2phase_viscoelastic.inp', MS(sp.rand(31,31,31),2,1,1,1,1), True)
#generateAbaqusInp('41_2phase_viscoelastic.inp', MS(sp.rand(41,41,41),2,1,1,1,1), True)
