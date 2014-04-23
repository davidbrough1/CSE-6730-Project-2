# -*- coding: utf-8 -*-
"""
3D, Isotropic, 2st Order MKS

This script calibrates against reference datasets and plots the FE and MKS
response for a validation microstructure.

Noah Paulson, 4/8/2014
"""

import os
import time
import numpy as np
import matplotlib.pyplot as plt
import MKS as MKS

def independent_columns(A, tol = 1e-05):
    """
    This function is from: http://stackoverflow.com/q/13312498

    Returns an array composed of independent columns of A.

    Note that answer may not be unique; this function returns one of many
    possible answers.
    """
    Q, R = np.linalg.qr(A)
    independent = np.where(np.abs(R.diagonal()) > tol)[0]
    #return A[:, independent]
    return independent

def res_red(filename = "21_1_noah.dat", el = 21, slc = 10):
    """
    This function reads the E11 values from a .dat file and reorganizes the
    data into a el x el x el array with the correct organization
    
    It will also plot a certain x-slice in the dataset if called within this
    script
    """
    f = open(filename, "r")

    linelist = f.readlines()

    # finds a location several lines above the start of the data
    # linelist[n] reads the entire line at location n
    for n in range(1000):
        if 'THE FOLLOWING TABLE' in linelist[n]:
            break

    # line0 is the index of first line of the data
    line0 = n + 5;      

    e11 = np.zeros((21**3,8))
    c = 0

    # this series of loops generates a 9261x8 dataset of E11s (element x integration point) 
    for ii in range(21**3):
        for jj in range(8):
            e11pre = linelist[line0 + c].split()[2]
            c = c + 1
            e11[ii,jj] = float(e11pre)
    
    f.close()    
    
    # here we average all 8 integration points in each element cell
    e11cond = np.flipud(np.mean(e11, axis=1))
    # element 4630 is at the centroid of a 21x21x21 dataset
    #print e11cond[4630]

    # here we reshape the data from a 9261 length vector to a 21x21x21 3D matrix
    pre_e11mat = np.reshape(e11cond, [21,21,21])
    e11mat = np.swapaxes(pre_e11mat,1,2)

    

    # here we generate an image and scalebar for a particular slice of the 3D matrix
    if __name__ == '__main__':
        
        plt.clf()        
                
        plt.imshow(e11mat[slc,:,:], origin='lower', interpolation='none',
                   cmap='jet')
        
        plt.colorbar()
        
    else:
        return e11mat

def pha_loc(filename = "msf.txt", el = 21, ns = 3):   
    """
    Opens file with microstructure info where is column represents a single
    microstructure. Converts each microstructure column vector into a 
    cube data structure with indexing which matches that of the 
    Finite Element structure
    
    Inputs: 'filename', 'ns'== number of microstructures

    Output: 21x21x21xns array where ns is the microstructure index    
    
    """
    
    f = open(filename, "r")

    linelist = f.readlines()
    
    pre_micr1 = np.zeros((21**3,ns))    
    
    for ii in range(21**3):
        for jj in range(ns):
            pre_micr1[ii,jj] = linelist[ii].split()[jj]

    f.close()
            
    # element 4630 is at the centroid of a 21x21x21 dataset
    #print e11cond[4630]

    # here we reshape the data from a 9261 length vector to a 21x21x21 3D matrix
        
    pre_micr2 = np.zeros((21,21,21,ns))
    micr = pre_micr2
    
    for jj in range(ns):    
        pre_micr2 = np.reshape(np.flipud(pre_micr1[:,jj]), [21,21,21])
        micr[:,:,:,jj] = np.swapaxes(pre_micr2,1,2)
        
    return micr

def remzer(r_ini):  
    """
    This function shrinks a vector by removing the zeros.
    NOTE: Find a better way to do this    
    """
    c = 0
    r = np.zeros(len(r_ini))
    for ii in range(len(r_ini)):
        if r_ini[ii] != 0:
            
            r[c] = r_ini[ii]        
            
            c = c + 1
    
    return np.trim_zeros(r)

el = 21
# ns refers to calibration datasets + VALIDATION DATASET
ns = 2
# H is the number of conformations of location and phase
H = 2

### THE MICROSTRUCTURE FUNCTION ###

## import delta and random microstructures
start = time.time()

#micr = pha_loc("cont5_micr_old.txt", el, ns)
micr = np.ndarray(shape =(21,21,21,2),dtype = 'float')
micr[:,:,:,0] = np.zeros((21,21,21),dtype = 'float')
micr[10,10,10,0] = 1
micr[:,:,:,1] = np.ones((21,21,21),dtype = 'float')
micr[10,10,10,1] = 0

end = time.time()
timeE = np.round((end - start),2)
print "Import microstructures: %s seconds" %timeE

# Black cells = 1, White cells = 0
# The black delta has a high stiffness cell surrounded by low stiffness cells,
# The white delta has a low stiffness cell surrounded by high stiffness cells

## microstructure functions
start = time.time()

m = np.zeros([el,el,el,ns,2])
m[:,:,:,:,0] = (micr == 0)
m[:,:,:,:,1] = (micr == 1)
m = m.astype(int)

end = time.time()
timeE = np.round((end - start),2)
print "Microstructure function generation: %s seconds" %timeE
    
## Microstructure functions in frequency space
start = time.time()

M = np.zeros((el, el, el, ns, H)) + 0j*np.zeros((el, el, el, ns, H))
for n in xrange(ns):
    for h in xrange(H):
        M[:,:,:,n,h] = np.fft.fftn(m[:,:,:,n,h])

end = time.time()
timeE = np.round((end - start),2)
print "Convert microstructure function to frequency space: %s seconds" %timeE


### FINITE ELEMENT RESPONSES ###

## responses of the black and white delta microstructure and a random
## microstructure.

# if read_dat == 1 the program will reload all of the .dat files and save them
# to FE_results.npy
read_dat = 1
if read_dat == 1:
    start = time.time()    
    
    #os.chdir("C:\mks_data\iso_2ph_5cont_dat")
    filename1 = "deltaM_modeling/21deltaM_1surroundedBy2_strain.dat"
    filename2 = "deltaM_modeling/21deltaM_2surroundedBy1_strain.dat"
    resp = np.zeros((el,el,el,ns))
    resp[:,:,:,0] = res_red(filename1)
    resp[:,:,:,1] = res_red(filename2)
    #print "%s is loaded" %filename 
    
    #os.chdir("C:/Users/nhpnp3/Documents/GitHub/MKS_repository/MKS_2nd")
    #np.save('FE_results',resp)    
    
    end = time.time()
    timeE = np.round((end - start),1)
    #print "Import FE results: %s seconds" %np.round((end - start),1)

# if read_dat == 0 the script will simply reload the results from a previously
# saved FE_results.npy
else:
    resp = np.load('FE_results.npy')
    print "FE results loaded"    

## responses in frequency space
start = time.time()

resp_fft = np.zeros((el,el,el,ns)) + 0j*np.zeros((el,el,el,ns))
for n in xrange(ns):
    resp_fft[:,:,:,n] = np.fft.fftn(resp[:,:,:,n])
    
end = time.time()
timeE = np.round((end - start),3)
print "Convert FE results to frequency space: "
print "%s seconds" %np.round((end - start),3)

### CALIBRATION OF INFLUENCE COEFFICIENTS ###

start = time.time()

#specinfc = np.zeros((el**3,H)) + 0j*np.zeros((el**3,H))
specinfc = np.zeros((el**3,H),dtype = 'complex128')
for k in xrange(el**3):
    
    [u,v,w] = np.unravel_index(k,[el,el,el]) #Index in a vector

    #The following 2 quantities are used in the normal equation during regression.
    MM = np.zeros((H,H)) + 0j*np.zeros((H,H))# Microstucture matrix (M' M*) 
    PM = np.zeros((H,1)) + 0j*np.zeros((H,1))# Property matrix (P M*) 
    
    for n in xrange(ns-1):
        mSQc = np.conjugate(M[u,v,w,n,:])  # Conjugate of FFT of Microstructure     
        mSQt = np.mat(M[u,v,w,n,:]).T      # Transpose of FFT of Microstructure
        
        MM = MM + np.outer(mSQt,mSQc)        # Calculate MM
        PM = PM + (resp_fft[u,v,w,n] * mSQc) # Calculate PM
        
    if k < 2:
        p = independent_columns(MM, .001)
    
    calred = MM[p,:][:,p]      # Linearly independent columns of MM
    resred = PM[p,0].conj().T  # Linearly independent columns of PM 

    ## for examining variables at desired frequency    
#    if k == 7:
#        print calred
#        print resred
#        print np.linalg.solve(calred,resred)
#        print specinfc[k, range(len(p))]

    if k % 500 == 0:
        print "frequency completed: %s" % k
    
    specinfc[k, p] = np.linalg.solve(calred, resred)   

end = time.time()
timeE = np.round((end - start),3)
print "Calibration: %s seconds" %np.round((end - start),1)

### VALIDATION WITH RANDOM ARRANGEMENT ###

## vectorize the frequency-space microstructure function for the validation
## dataset
lin_M = np.zeros((el**3,H)) + 0j*np.zeros((el**3,H))
for h in xrange(H):
    lin_M[:,h] = np.reshape(M[:,:,:,-1,h],el**3)


## find the frequency-space response of the validation microstructure
## and convert to the real space
lin_sum = np.sum(np.conjugate(specinfc) * lin_M, 1)
mks_F = np.reshape(lin_sum,[21,21,21])
mks_R = np.fft.ifftn(mks_F).real

np.save('MKS_1stOrd_resp',mks_R) 

### MEAN ABSOLUTE STRAIN ERROR (MASE) ###

MASE = 0
avgE = np.average(resp[:,:,:,-1])
print "The average strain is %s" %avgE

for k in xrange(el**3):
    
    [u,v,w] = np.unravel_index(k,[el,el,el])
    MASE = MASE + ((np.abs(resp[u,v,w,-1] - mks_R[u,v,w]))/(avgE * el**3))

print "The mean absolute strain error (MASE) is %s%%" %(MASE*100)

plt.imshow(resp[:,:,:,-1])
plt.colorbar()
plt.show()
plt.imshow(mks_R)
plt.colorbar()
plt.show()


## VISUALIZATION OF MKS VS. FEM ###

plt.close()
#fig = plt.figure()

## pick a slice perpendicular to the x-direction
slc = 10

## find the min and max of both datasets (needed to scale both images the same) 
dmin = np.amin([resp[slc,:,:,-1],mks_R[slc,:,:]])
dmax = np.amax([resp[slc,:,:,-1],mks_R[slc,:,:]])

plt.subplot(221)
ax = plt.imshow(mks_R[slc,:,:], origin='lower', interpolation='none',
    cmap='jet', vmin=dmin, vmax=dmax)
plt.colorbar(ax)
plt.title('MKS approximation, E11')

plt.subplot(222)
ax = plt.imshow(resp[slc,:,:,-1], origin='lower', interpolation='none',
    cmap='jet', vmin=dmin, vmax=dmax)
plt.colorbar(ax)
plt.title('FE response, E11')   

plt.subplot(212)

feb = remzer(np.reshape(resp[:,:,:,-1]*m[:,:,:,-1,1],el**3))
few = remzer(np.reshape(resp[:,:,:,-1]*m[:,:,:,-1,0],el**3))
mksb = remzer(np.reshape(mks_R[:,:,:]*m[:,:,:,-1,1],el**3))
mksw = remzer(np.reshape(mks_R[:,:,:]*m[:,:,:,-1,0],el**3))

bn = 40

n, bins, patches = plt.hist(feb, bins = bn, histtype = 'step', hold = True,
                            range = (dmin, dmax), color = 'white')
bincenters = 0.5*(bins[1:]+bins[:-1])
febp, = plt.plot(bincenters,n,'k', linestyle = '--', lw = 1.5)

n, bins, patches = plt.hist(few, bins = bn, histtype = 'step', hold = True,
                            range = (dmin, dmax), color = 'white')
fewp, = plt.plot(bincenters,n,'k')

n, bins, patches = plt.hist(mksb, bins = bn, histtype = 'step', hold = True,
                            range = (dmin, dmax), color = 'white')
mksbp, = plt.plot(bincenters,n,'r', linestyle = '--', lw = 1.5)

n, bins, patches = plt.hist(mksw, bins = bn, histtype = 'step', hold = True,
                            range = (dmin, dmax), color = 'white')
mkswp, = plt.plot(bincenters,n,'r')

plt.grid(True)

plt.legend([febp,fewp,mksbp,mkswp], ["FE - stiff phase", 
           "FE - compliant phase", "MKS - stiff phase",
           "MKS - compliant phase"])

plt.xlabel("Strain")
plt.ylabel("Frequency")
plt.title("Frequency comparison of FE and MKS")
plt.show()

#
#del MM, PM, M, lin_M, specinfc, lin_sum, mks_F, m, micr, pm
#del febp, fewp, mksbp, mkswp, ax
#
#plt.subplot(221)
#ax = plt.imshow(micr[slice,:,:,0], origin='lower', interpolation='none',
#    cmap='binary')
#plt.colorbar(ax)
#plt.title('Black delta microstructure')
#
#plt.subplot(222)
#ax = plt.imshow(micr[slice,:,:,1], origin='lower', interpolation='none',
#    cmap='binary')
#plt.colorbar(ax)
#plt.title('White delta microstructure')
#
#plt.subplot(223)
#ax = plt.imshow(micr[slice,:,:,2], origin='lower', interpolation='none',
#    cmap='binary')
#plt.colorbar(ax)
#plt.title('Validation microstructure')
#
#resp_fft_lin = np.reshape(resp_fft[:,:,:,-1],el**3).real
#freq = range(el**3)
#plt.plot(freq,resp_fft_lin,'b')

