#!/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from netCDF4 import Dataset

dt = 5. #grid spacing in km
            
readFile = "/users/global/danbel/wavelet/Conni/pcp_max_efolding01_tmp.npz"
wav = np.load(readFile)
xpcp = (wav['xpcp'])
ypcp = (wav['ypcp'])
xrad1 = (wav['xrad1'])
xrad2 = (wav['xrad2'])
yrad1 = (wav['yrad1'])
yrad2 = (wav['yrad2'])
houror = (wav['hourMax'])
pcpMax = (wav['pcpMax'])
tmpmin = (wav['tmpmin'])

ffile = "/users/global/cornkle/tp_99th.nc"
fh = Dataset(ffile, mode='r')
pcpC = fh.variables['p'][:]
tmpC = fh.variables['t'][:]

ind = (xrad1>-1)&(xrad2<998)&(yrad1>-1)&(yrad2<998)&(pcpMax>0.)

pcpMaxind = pcpMax[ind]
tmpminind = tmpmin[ind]
hour = houror[ind]

hr1 = 0
hr2 = 24
#stratify for hours
indhr = (hour>=hr1) & (hour<=hr2)

#Size of systems
Lx = xrad2[ind] + xrad1[ind]
Ly = yrad2[ind] + yrad1[ind]
Lx = Lx * dt
Ly = Ly * dt

#Histogram
#bins = np.unique(Ly[indhr])[:-2] - 0.1
bins = np.arange(9.9,39.9,5)
bins = np.append(bins,bins[-1]+0.5)
#hist, bin_edges = np.histogram(Lscale, bins=bins)

#Bin average
bin_means, bin_edges, binnumber = stats.binned_statistic(
#            hour[indhr], pcpMaxind[indhr], bins=24, statistic=stat)
            Lx[indhr], pcpMaxind[indhr], bins=bins, statistic='mean')
#            hour, pcpMax, bins=24, statistic=stat)
bin_meansy, bin_edgesy, binnumbery = stats.binned_statistic(
#            hour[indhr], pcpMaxind[indhr], bins=24, statistic=stat)
            Ly[indhr], pcpMaxind[indhr], bins=bins, statistic='mean')
#            hour, pcpMax, bins=24, statistic=stat)

#Confidence Interval
bin_std, bin_edges, binnumber = stats.binned_statistic(
            Lx[indhr], pcpMaxind[indhr], bins=bins, statistic=np.std)
bin_count, bin_edges, binnumber = stats.binned_statistic(
            Lx[indhr], pcpMaxind[indhr], bins=bins, statistic='count')
R = stats.t.interval(0.95,bin_count-1,loc=bin_means,scale=bin_std/np.sqrt(bin_count))
bin_stdy, bin_edgesy, binnumbery = stats.binned_statistic(
            Ly[indhr], pcpMaxind[indhr], bins=bins, statistic=np.std)
bin_county, bin_edgesy, binnumbery = stats.binned_statistic(
            Ly[indhr], pcpMaxind[indhr], bins=bins, statistic='count')
Ry = stats.t.interval(0.95,bin_county-1,loc=bin_meansy,scale=bin_stdy/np.sqrt(bin_county))

#Plot:
#T-pcp scatter
plt.scatter(tmpminind[~np.isnan(tmpminind)],pcpMaxind[~np.isnan(tmpminind)])
plt.scatter(-tmpC,pcpC,color='r')
plt.show()
#bin-averaged
#plt.scatter(Ly[indhr],pcpMaxind[indhr])
#plt.scatter(hour[indhr],pcpMaxind[indhr])
plt.errorbar(bin_edges[:-1]-0.2-5, bin_means, yerr=R[1]-bin_means, fmt='-o',
             lw=2, ms=10, label='Lx')
plt.errorbar(bin_edgesy[:-1]+0.4-5, bin_meansy, yerr=Ry[1]-bin_meansy, fmt='-d',
             color='r', lw=2, ms=10, label='Ly')
#plt.plot(bin_edges[:-1],bin_means,'r',linewidth=2)
plt.xlabel('L (km)',fontsize=20)
#plt.xlabel('Time')
plt.ylabel('P intensity (mm/h)',fontsize=20)
plt.xlim([0,35])
plt.legend(loc='upper left')
plt.tick_params(axis='both', which='major', labelsize=18)
#plt.ylim([10,60])
plt.show()
