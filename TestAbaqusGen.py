import scipy as sp
import numpy as np
from MS import MS
from MSf import MSf
from AbaqusGen import generateAbaqusInp
from AbaqusGen import getEleNumber

print 'Ele Number for 0, 0, 0: {}'.format(getEleNumber(0, 0, 0, 21, 21, 21, 1))
print 'Ele Number for 0, 20, 0: {}'.format(getEleNumber(0, 20, 0, 21, 21, 21, 1))
print 'Ele Number for 1, 0, 0: {}'.format(getEleNumber(1, 0, 0, 21, 21, 21, 1))
print 'Ele Number for 0, 0, 1: {}'.format(getEleNumber(0, 0, 1, 21, 21, 21, 1))
print 'Ele Number for 20, 20, 20: {}'.format(getEleNumber(20, 20, 20, 21, 21, 21, 1))

ms = MS(sp.rand(21,21,21),2,1,1,1,1)

generateAbaqusInp('testFile.inp', ms)