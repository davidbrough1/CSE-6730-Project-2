from AbaqusGen import generateAbaqusInp
from MS import MS
import numpy as np
import matplotlib.pylab as plt

A = np.ndarray(shape = (63,63,63),dtype = 'float')
Ms = MS(A,2,1,1,1,1)
plt.imshow(Ms[31,:,:])
plt.show()
generateAbaqusInp('63_2phase',Ms,viscoelastic=False)

