
import numpy as np
from scipy.ndimage import gaussian_filter, uniform_filter
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from netCDF4 import Dataset
from scipy.interpolate import griddata

#read the hurtt change fractions
fdfrac = "dfrac_primary_land_2010_1950.nc"
 
ncid = Dataset(fdfrac)
dfrac = ncid.variables["absolute_dfrac"][:,:]
seamask = ncid.variables["lsm"][:,:]

vlon = ncid.variables["longitude"][:]
vlat = ncid.variables["latitude"][:]
mlon, mlat = np.meshgrid(vlon, vlat)

seamask = seamask[350:420,620:800]
a = dfrac[350:420,620:800]

#get an array for interpolation
a2 = a.copy()
a2[seamask==0]=np.nan

points = np.where(np.isfinite(a2))
inter = np.where(np.isnan(a2))

#interpolate over sea from land points
a2[inter] = griddata(points, np.ravel(a2[points]), inter, method = 'nearest')

#apply gaussian filter to interpolated image
a3 = gaussian_filter(a2,sigma = 0.8 )
a3[seamask==0]=np.nan
print('gaussian interpolated: ',np.nanmean(a3),np.nansum(a3), np.nanmax(a3), np.nanmin(a3))

#apply gaussian filter to original image
a4 = gaussian_filter(a, sigma = 0.8)
a4[seamask==0]=np.nan
print('gaussian no interpolation: ',np.nanmean(a4),np.nansum(a4), np.nanmax(a4), np.nanmin(a4))

#mask out sea on original image
a[seamask==0]=np.nan
print('original: ',np.nanmean(a),np.nansum(a), np.nanmax(a), np.nanmin(a))

f = plt.figure(figsize=(10,6))
ax=f.add_subplot(321,title = "original")
plt.imshow(a, origin='lower', cmap='jet', vmin=-0.7)
plt.colorbar()

ax=f.add_subplot(322, title = "gaussian sigma = 0.8, Interpolated (nearest)")
plt.imshow(a3, origin='lower', cmap='jet', vmin=-0.7)
plt.colorbar()

ax=f.add_subplot(323, title = "gaussian sigma = 0.8, No interpolation")
plt.imshow(a4, origin='lower', cmap='jet', vmin=-0.7)
plt.colorbar()

ax=f.add_subplot(324, title = "Interpolated gaussian - gaussian")
plt.imshow(a3-a4, origin='lower', cmap='viridis_r', vmin = -0.03)
plt.colorbar()

ax=f.add_subplot(325, title = "original - interpolated gaussian")
plt.imshow(a-a3, origin='lower', cmap='RdBu', vmin=-0.05, vmax=0.05)
plt.colorbar()

plt.tight_layout()
plt.show()

