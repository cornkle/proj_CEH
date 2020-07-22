# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 10:15:40 2016

@author: cornkle
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib
import ipdb
import pandas as pd
from collections import OrderedDict
import salem
from utils import u_met, u_parallelise, u_arrays as ua, constants as cnst, u_darrays
from scipy.interpolate import griddata
import multiprocessing
import glob
from scipy import ndimage
from utils import u_statistics as ustats
import salem
from metpy import calc
from metpy.units import units


import pickle as pkl


matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)



dic = pkl.load(open(cnst.network_data + "figs/LSTA/paper/saves/fig01.p", "rb"))

extent = (dic['lsta0'].shape[1]-1)/2
xlen = dic['lsta0'].shape[1]
ylen = dic['lsta0'].shape[0]

xv, yv = np.meshgrid(np.arange(ylen), np.arange(xlen))
st=30
xquiv = xv[4::st, 4::st]
yquiv = yv[4::st, 4::st]

u = (dic['u925'])[4::st, 4::st]
v = (dic['v925'])[4::st, 4::st]
levels = list(np.arange(-2.7,-0.1,0.3)) + list(np.arange(0.3,2.71,0.3))
qlevels = list(np.round(np.arange(-0.55, -0.01, 0.05),2)) + list(np.round(np.arange(0.05, 0.56, 0.05),2))
divlevels = list(np.arange(-0.8, -0.01, 0.1)) + list(np.arange(0.1, 0.805, 0.1))

f = plt.figure(figsize=(15.8,4), dpi=250)


ax1 = f.add_subplot(131)

plt.contourf(dic['lsta0'], extend='both', cmap='RdBu', levels=levels)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=12)
cbar.ax.set_ylabel(r'%', size=12)

cb=plt.contourf(dic['lsta0_sig'], colors='none', hatches='.',
             levels=[0.5, 1])
for i, collection in enumerate(cb.collections):
    collection.set_edgecolor('peachpuff')
for collection in cb.collections:
    collection.set_linewidth(0.)


lev = [-0.6, -0.5, -0.4,-0.3, -0.2, 0, 0.2,0.3, 0.4, 0.5,0.6]

contours = plt.contour(dic['t'], extend='both', #(dic['t']-dic['tclim'])
                       levels=lev,
                       colors='k',linestyles = 'solid', linewidths = 1)

plt.clabel(contours, inline=True, fontsize=11, fmt='%1.1f')

plt.axvline(x=200, linestyle='dashed', color='dimgrey', linewidth=1.2)
plt.axhline(y=200, linestyle='dashed', color='dimgrey', linewidth=1.2)
plt.plot(extent, extent, marker='o', color='dimgrey')
ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -200) * 3, dtype=int))
ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) -200) * 3, dtype=int))
ax1.set_xlabel('km', fontsize=12)
ax1.set_ylabel('km', fontsize=12)


ax1 = f.add_subplot(133)

plt.contourf(dic['q'], extend='both',  cmap='RdBu',levels=qlevels) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=12)
cbar.ax.set_ylabel(r'g kg$^{-1}$', size=12)

cb=plt.contourf(dic['q_sig'], colors='none', hatches='.',
             levels=[0.5, 1])
for i, collection in enumerate(cb.collections):
    collection.set_edgecolor('peachpuff')
for collection in cb.collections:
    collection.set_linewidth(0.)

contours = plt.contour(dic['shear'] , extend='both',levels=np.arange(-15.25,-12.75,0.5), colors='k', linestyles='solid', linewidths=1) #np.arange(-15,-10,0.5)

plt.clabel(contours, inline=True, fontsize=11, fmt='%1.2f')
plt.axvline(x=200, linestyle='dashed', color='dimgrey', linewidth=1.2)
plt.axhline(y=200, linestyle='dashed', color='dimgrey', linewidth=1.2)
plt.plot(extent, extent, marker='o', color='dimgrey')

ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -200) * 3, dtype=int))
ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) -200) * 3, dtype=int))
ax1.set_xlabel('km', fontsize=12)

ax1 = f.add_subplot(132)
div = ndimage.gaussian_filter(((dic['div'])/ dic['cnte2'])*100, 3, mode='nearest')
vals = np.percentile(div,[9,91])
#ipdb.set_trace()
plt.contourf(div, extend='both',  cmap='RdBu', levels=divlevels) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
mask =  (div <= vals[0])
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=12)
cbar.ax.set_ylabel(r'10$^{-3}$s$^{-1}$', size=11)

# plt.contour(div, extend='both', colors='k',
#              levels=divlevels, linewidths=0.05)

cb = plt.contourf(np.arange(mask.shape[1])[::15], np.arange(mask.shape[0])[::15], mask[::15,::15], colors='none', hatches='.',
             levels=[0.5, 1])
for i, collection in enumerate(cb.collections):
    collection.set_edgecolor('peachpuff')
for collection in cb.collections:
    collection.set_linewidth(0.)

v925 = ndimage.gaussian_filter(dic['v925_orig'] / dic['cnte'], 3, mode='nearest')
plt.contour(v925, extend='both', colors='k', linewidths=5,
                       levels=[-50, 0.06, 50])

#itd = ndimage.gaussian_filter(dic3['itd'] / dic3['cnte'], 8, mode='nearest')

#
#itd = xr.DataArray(np.array((dic4['itd']/dic4['cnte'])[1::,1::]).astype(float))
itd = xr.DataArray(np.array((dic['itd'] / dic['cnte4'])).astype(float))
# cnt = xr.DataArray(np.array(dic3['cnte'][1::, 1::]).astype(float))
#
#itd = salem.reduce(itd, factor=8,how=np.sum)
# cnt = salem.reduce(cnt, factor=16,how=np.max)

itd = ndimage.gaussian_filter((itd)*100, 28, mode='nearest')

#ipdb.set_trace()
#cs = plt.contour(np.arange(0,400, 8),np.arange(0,400,8),itd, extend='both', colors='r', linewidths=2) #
cs = plt.contour( itd, extend='both', colors='k', linewidths=0.5, levels=[24.9,30,35], linestyles='dashed')
cs = plt.contourf(itd, colors='lightgrey',levels=[24.9, 31, 35], alpha=0.27)
plt.clabel(cs, inline=True, fontsize=11, fmt='%1.0f')

plt.axvline(x=200, linestyle='dashed', color='dimgrey', linewidth=1.2)
plt.axhline(y=200, linestyle='dashed', color='dimgrey', linewidth=1.2)
plt.plot(extent, extent, marker='o', color='dimgrey')
# contours = plt.contour((dic['shear'] / dic['cntp']), extend='both',levels=np.arange(-17,-12,0.5), cmap='viridis') #np.arange(-15,-10,0.5)
# plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
qu = ax1.quiver(xquiv, yquiv, u, v, scale=15, headwidth=5)
qk = plt.quiverkey(qu, 0.9, 0.03,1, r'1 m s$^{-1}$',
                   labelpos='E', coordinates='figure')

ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -200) * 3, dtype=int))
ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) -200) * 3, dtype=int))
ax1.set_xlabel('km', fontsize=12)
#ax1.set_ylabel('km')
#plt.title('shading: divergence, vectors: 925hPa wind anomaly', fontsize=9)

plt.tight_layout()
text = ['a', 'b', 'c']
plt.annotate(text[0], xy=(0.006, 0.92), xytext=(0, 4),xycoords=('figure fraction', 'figure fraction'),
             textcoords='offset points', fontweight='bold', fontname='Ubuntu', fontsize=16)
plt.annotate(text[1], xy=(0.33, 0.92), xytext=(0, 4), xycoords=('figure fraction', 'figure fraction'),
             textcoords='offset points', fontweight='bold', fontname='Ubuntu', fontsize=16)
plt.annotate(text[2], xy=(0.66, 0.92), xytext=(0, 4),  xycoords=('figure fraction', 'figure fraction'),
             textcoords='offset points', fontweight='bold', fontname='Ubuntu', fontsize=16)


plt.savefig(cnst.network_data + "figs/LSTA/paper/plots/fig01.png")
plt.close()
