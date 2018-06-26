#!/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
            
readFile = "/users/global/cornkle/pcp_max_efolding01_tmp.npz"
wav = np.load(readFile)
xrad1 = (wav['xrad1'])
xrad2 = (wav['xrad2'])
yrad1 = (wav['yrad1'])
yrad2 = (wav['yrad2'])
pcpMax = (wav['pcpMax'])
tmpmin = (wav['tmpmin'])

ffile = "/users/global/cornkle/rubbish/tp_99th.nc"
fh = Dataset(ffile, mode='r')
pcpC = fh.variables['p'][:]
tmpC = fh.variables['t'][:]

ind = (xrad1>-1)&(xrad2<998)&(yrad1>-1)&(yrad2<998)&(pcpMax>0.)

pcpMaxind = pcpMax[ind]
tmpminind = tmpmin[ind]

#Plot:
#T-pcp scatter
plt.figure()
plt.scatter(tmpminind[~np.isnan(tmpminind)],pcpMaxind[~np.isnan(tmpminind)])
#plt.scatter(-tmpC,pcpC,color='r')
plt.show()

