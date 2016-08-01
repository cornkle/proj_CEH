#!/bin/env python

from scipy import ndimage

#Find maxima of the 2D wavelet spectrum

#TIR
maxpowerTIR = (powerTIR == ndimage.maximum_filter(powerTIR,size=(5,5,5),
           mode='constant',cval=np.amax(powerTIR)+1))
maxpowerTIR = maxpowerTIR.astype(int)
powerTIRthresh = np.repeat(50*(scales2d/1.)**.5,powerTIR.shape[1],axis=1) #NOTICE the scale-dependant threshold for power
powerTIRthresh = np.repeat(powerTIRthresh,powerTIR.shape[2],axis=2)
powerTIRbkg = (powerTIR < powerTIRthresh).astype(int) #Apply the threshold
powerTIRpeaks = maxpowerTIR - powerTIRbkg  #Remove all smaller than the threshold
zpks,ypks,xpks = np.argwhere(powerTIRpeaks == 1).T
