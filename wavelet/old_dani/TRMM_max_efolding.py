#!/bin/env python

import numpy as np
#import matplotlib.pyplot as plt
from netCDF4 import Dataset
from wavelet import twod as w2d
import os
from scipy import ndimage

def locate(pattern, root_path):
    # pattern: a pattern specifying the end of file (file type you are looking for, eg.  ".nc")
    # root_path: full path to directory that should be searched
    strg=[]
    flist = os.listdir(root_path)
    flist.sort()
    for ffile in flist:
        if ffile.endswith(pattern):
            strg.append(os.path.join(root_path, ffile))
    return strg

dt = 5. #grid spacing in km

pcpMax = []
tmpmin = []
hourMax = []
xpcp = []
ypcp = []
xrad1 = []
xrad2 = []
yrad1 = []
yrad2 = []
#Read Conni's files
ccfile = locate(".nc","/users/global/cornkle/TRMMfiles")
for ffile in ccfile:
    fh = Dataset(ffile, mode='r')
    pcp = fh.variables['p'][:]
    lat = fh.variables['lat'][:]
    lon = fh.variables['lon'][:]
    tmp = fh.variables['t'][:]
    hour = int(ffile[-11:-9]) + int(ffile[-8:-6])/60.
    pcp[np.isnan(pcp)] = 0.0
#    pcp[pcp<0.1] = 0.0
    if (np.max(pcp) < 0.1):
        continue
    print(ffile)
    pcpMax = np.append(pcpMax, np.amax(pcp))
    hourMax = np.append(hourMax, hour)
    locpcp = np.unravel_index(np.argmax(pcp),pcp.shape)
    ypcptmp = locpcp[0]
    xpcptmp = locpcp[1]
    ypcp = np.append(ypcp,locpcp[0])
    xpcp = np.append(xpcp,locpcp[1])
    for i in range(1,30):
        if (xpcptmp-i < 0):
            xrad1t = -1
            break
        if pcp[ypcptmp,xpcptmp-i]<np.exp(-1)*np.amax(pcp):
            xrad1t = i
            break
    for i in range(1,30):
        if (xpcptmp+i > pcp.shape[1]-1):
            xrad2t = 999
            break
        if pcp[ypcptmp,xpcptmp+i]<np.exp(-1)*np.amax(pcp):
            xrad2t = i
            break
    for j in range(1,30):
        if (ypcptmp-j < 0):
            yrad1t = -1
            break
        if pcp[ypcptmp-j,xpcptmp]<np.exp(-1)*np.amax(pcp):
            yrad1t = j
            break
    for j in range(1,30):
        if (ypcptmp+j > pcp.shape[0]-1):
            yrad2t = 999
            break
        if pcp[ypcptmp+j,xpcptmp]<np.exp(-1)*np.amax(pcp):
            yrad2t = j
            break
    xrad1 = np.append(xrad1,xrad1t)
    xrad2 = np.append(xrad2,xrad2t)
    yrad1 = np.append(yrad1,yrad1t)
    yrad2 = np.append(yrad2,yrad2t)
    if ( (xrad1t>-1) & (xrad2t<900) & (yrad1t>-1) & (yrad2t<900) ):
        tmpmin = np.append(tmpmin, np.amin(tmp[ypcptmp-yrad1t:ypcptmp+yrad2t, 
                                               xpcptmp-xrad1t:xpcptmp+xrad2t]))
    else:
        tmpmin = np.append(tmpmin, np.nan)
            
saveFile = "/users/global/cornkle/pcp_max_efolding01_tmp"
np.savez(saveFile, pcpMax=pcpMax, tmpmin=tmpmin, hourMax=hourMax, 
         xpcp=xpcp, ypcp=ypcp, xrad1=xrad1, xrad2=xrad2, 
         yrad1=yrad1, yrad2=yrad2)
