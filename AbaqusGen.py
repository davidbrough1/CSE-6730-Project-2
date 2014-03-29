import scipy as sp
import numpy as np
from matplotlib.cbook import flatten

def generateAbaqusIn(inputFileName, MicroSF):
    f = open(inputFileName, 'w')
    nl = "\n"
    headerLines = '*Preprint, echo=NO, model=No, history=NO, contact=NO', nl, '*Heading', nl
    headerLines += '****************************************************', nl
    headerLines += '*node', nl
                  
    f.writelines(headerLines)
    
    Length_X = 2
    Length_Y = 2
    Length_Z = 2
    
    intx = 10
    inty = 10
    intz = 10
    
    Fine_int = 1
    
    blockNo = 1
    NodeNo1 = 1
    ElementNo1 = 1
    x0 = 0
    y0 = 0
    z0 = 0
    
    setA = NodeNo1
    setC = setA + (inty*Fine_int)
    setB = setA + (inty*Fine_int + 1) * (intz*Fine_int)
    setD = setA + (inty*Fine_int + 1) * (intz*Fine_int) + (inty*Fine_int)
    
    LayerNodes = (inty*Fine_int + 1) * (intz*Fine_int) + (inty*Fine_int) + 1
    
    setE = (intx*Fine_int)*LayerNodes + setA
    setG = setE + (inty*Fine_int)
    setF = setE + (inty*Fine_int + 1)*(intz*Fine_int)
    setH = setE + (inty*Fine_int + 1)*(intz*Fine_int) + (inty*Fine_int)
    
    Total_Nodes = setH

    modelCoordinates = [[setA, x0, y0, z0], [setB, x0, y0, z0 + Length_Z],
                        [setC, x0, y0 + Length_Y, z0],[setD, x0, y0 + Length_Y, z0 + Length_Z],
                        [setE, x0 +Length_X, y0, z0], [setF, x0 + Length_X, y0, z0],
                        [setG, x0 + Length_X, y0 + Length_Y, z0],[setH, x0 + Length_X, y0 + Length_Y, z0 + Length_Z]]   
    vertexLineFormat = '{}, {}, {}, {} \n'
    for i in range(len(modelCoordinates)):
        f.writelines(vertexLineFormat.format(*modelCoordinates[i]))
    
    nSetNames = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    nSetFormat = '*NSET, NSET={} \n {} \n'
    f.writelines(('*****************************************************', nl))
    for i in range(len(nSetNames)):
        f.writelines(nSetFormat.format(nSetNames[i], modelCoordinates[i][0]))
    
    # Generate nodes for main region
    # Below: [x y] = x for interval & y for numbering
    nFillDict = {}
    setAB_f = [intz, (inty*Fine_int + 1)*Fine_int]
    nFillDict['AB'] = setAB_f
    setCD_f = [intz, (inty*Fine_int + 1)*Fine_int]
    nFillDict['CD'] = setCD_f
    setABCD_f = [inty, Fine_int]
    nFillDict['ABCD'] = setABCD_f
    
    nFillFormat = '*Nfill, Nset={}\n{}, {}, {}, {}\n'
    f.writelines(('****************************************************',nl))
        
    setEF_f = [intz, (inty*Fine_int + 1)*Fine_int]
    nFillDict['EF'] = setEF_f
    setGH_f = [intz, (inty*Fine_int + 1)*Fine_int]
    nFillDict['GH'] = setGH_f
    setEFGH_f = [inty, Fine_int]
    nFillDict['EFGH'] = setEFGH_f
    
    setABCDEFGH_f = [intx, LayerNodes*Fine_int]
    nFillDict['ABCDEFGH'] = setABCDEFGH_f
    
#    for i in range(1, np.log2(len(nSetNames))):
#        windowSize = 2**i
#        for j in range(len(nSetNames)/windowSize):
#            nodeSet = []
#            for k in range(windowSize):
#                nodeSet += nSetNames[j*windowSize + k]
    f.write(nFillFormat.format('AB','A','B', setAB_f[0], setAB_f[1]))
    f.write(nFillFormat.format('CD','C','D', setCD_f[0], setCD_f[1]))
    f.write(nFillFormat.format('ABCD','AB','CD', setABCD_f[0], setABCD_f[1]))
    f.write(nFillFormat.format('EF','E','F', setEF_f[0], setEF_f[1]))
    f.write(nFillFormat.format('GH','G','H', setGH_f[0], setGH_f[1]))
    f.write(nFillFormat.format('EFGH','EF','GH', setEFGH_f[0], setEFGH_f[1]))
    f.write(nFillFormat.format('ABCDEFGH','ABCD','EFGH', setABCDEFGH_f[0], setABCDEFGH_f[1]))
    
    f.writelines(('****************************************************', nl))    
    f.write(nFillFormat.format('BD','B','D', setABCD_f[0], setABCD_f[1]))
    f.write(nFillFormat.format('FH','F','H', setABCD_f[0], setABCD_f[1]))
    f.write(nFillFormat.format('BDFH','BD','FH', setABCDEFGH_f[0], setABCDEFGH_f[1]))    
    
    f.writelines(('*****************************************************', nl))
    f.writelines(('** Bottom Face Node Set', nl))
    f.write(nFillFormat.format('AC','A','C', setABCD_f[0], setABCD_f[1]))
    f.write(nFillFormat.format('EG','E','G', setABCD_f[0], setABCD_f[1]))
    f.write(nFillFormat.format('ACEG','AC','EG', setABCDEFGH_f[0], setABCDEFGH_f[1]))  
    
    # Define Element#1
    Cube_1=setA
    Cube_2=setA+Fine_int
    Cube_3=setA+(inty*Fine_int+2)*Fine_int
    Cube_4=setA+(inty*Fine_int+1)*Fine_int
    Cube_5=(1*Fine_int)*LayerNodes+setA
    Cube_6=(1*Fine_int)*LayerNodes+setA+Fine_int
    Cube_7=Cube_5+(inty*Fine_int+2)*Fine_int
    Cube_8=Cube_5+(inty*Fine_int+1)*Fine_int
    
    Element1 = [ElementNo1, Cube_1, Cube_2, Cube_3, Cube_4, Cube_5, Cube_6, Cube_7, Cube_8]
    
    f.writelines(('*****************************************************', nl))
    f.writelines(('*ELEMENT, TYPE=C3D8',nl))
    f.writelines((','.join(map(str,Element1)), nl))
    
    
    # Generate Elements for main region
    # Below: [x y z] = x no of elements including master
    #                  y increment in node number
    #                  z increment in element number

    Row_ELM_f=[inty, Fine_int, Fine_int]
    Layer_ELM_f=[intz, setAB_f, (inty*Fine_int)*Fine_int]
    Block_ELM_f=[intx, setABCDEFGH_f, (inty*Fine_int)*(intz*Fine_int)*Fine_int]
    
    ELGEN=[ElementNo1, Row_ELM_f, Layer_ELM_f, Block_ELM_f]
    
    f.writelines(('*****************************************************', nl))
    f.writelines(('*ELGEN, elset=allel',nl))
    f.writelines((','.join(map(str,flatten(ELGEN))), nl))
    
    Total_Elements=intx*inty*intz
    
    
