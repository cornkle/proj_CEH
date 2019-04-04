#!/usr/bin/env python
# coding: utf-8

# In[53]:


import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from utils import u_grid
import ipdb

file = '/prj/vera/semval/hansen_forest_2018_0.05deg/Hansen_GFC_2000-2017-v1.5_treecover_WAfrica_250m_cdomerge_time.nc'

da = xr.open_dataset(file)
da = da['treecover']#.sel(LON=slice(-8,-5), LAT=slice(6,8))


# In[65]:

#
# plt.figure(figsize=(9,7))
# plt.imshow(da[-1], origin='lower')
# plt.colorbar()
#

# In[66]:


da.values[np.isnan(da.values)]=0


#ipdb.set_trace()

grid = u_grid.make(da['LON'].values, da['LAT'].values, 5000)
outt = grid.lookup_transform(da, return_lut=False, method=np.nanmean)

da_new = xr.DataArray(outt, coords={'time': da.time, 'lat': grid.ll_coordinates[1][:,0], 'lon': grid.ll_coordinates[0][0,:]},
                        dims=['time','lat', 'lon'])


da_new.to_netcdf('/prj/vera/cornkle/treefrac_5km.nc')

#
# plt.figure(figsize=(9,7))
# plt.imshow(outt[0], origin='lower')
# plt.colorbar()
#
#
# plt.figure(figsize=(9,7))
# plt.imshow(outt[-1], origin='lower')
# plt.colorbar()
#
#
# plt.figure(figsize=(9,7))
# plt.imshow(outt[0]-outt[-1], origin='lower', cmap='BrBG', vmin=-50, vmax=50)
# plt.colorbar()
#
# plt.figure(figsize=(9,7))
# plt.plot(grid.ll_coordinates[0][0,:],outt[0,25,:], label='forested')   # zonal crossection through deforestation
# plt.plot(grid.ll_coordinates[0][0,:],outt[-1,25,:], label='deforested')
# plt.legend()
# plt.ylabel('Forest cover fraction')
# plt.ylabel('Longitudes')




