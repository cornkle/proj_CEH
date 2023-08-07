#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 09:37:30 2019

@author: abamba
"""

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import glob
import xarray as xr
import pdb
import numpy.ma as ma


#from netCDF4 import Dataset
myfile = '/media/abamba/AMMA2050/AMMAGAP/ammagapdata/vegfrac/'
myfile2 = '/media/abamba/AMMA2050/AMMAGAP/ammagaptest/temp/'


veg = glob.glob(myfile+'qrparm.veg.frac_modis_inv_tewnf.nc')
lonlat = glob.glob(myfile2+'tewnfa.pb20110501.nc_STASH_m01s03i236_T_1_5m.nc')

dsv = xr.open_mfdataset(veg, decode_times=False)#, concat_dim='TH1_MN')
dsl = xr.open_mfdataset(lonlat, decode_times=False)#, concat_dim='TH1_MN')

lat = dsl.latitude_t
lon = dsl.longitude_t

#veg = dsv['field1391'].mean(dim=['pseudo'])#.mean(dim=['TH1_MN'])

lst = np.where((dsv['field1391'] > 0) & (dsv['field1391'] < 1))

pdb.set_trace()

#veg[0,:,:].values.max()
#veg[0,:,:].values.min()

f=plt.figure(figsize=(12, 9), dpi=80)
ax = f.add_subplot(111, projection=ccrs.Mercator())
plt.contourf(lon, lat, lst[0,:,:], extend='both', cmap='jet', transform=ccrs.PlateCarree())#, levels=np.arange(0,1,0.1)) 
ax.coastlines()
xl = ax.gridlines(draw_labels=True);
xl.xlabels_top = True 
xl.ylabels_right = True
plt.colorbar()
plt.title('January')


plt.savefig('vegmean.png')
#plt.yticks(y)
plt.tight_layout()
plt.show()


#pdb.set_trace()


