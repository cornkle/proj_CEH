import xarray as xr
import numpy as np
import ipdb
import matplotlib.pyplot as plt
from utils import u_plot

cube = xr.open_dataarray('/users/global/cornkle/MCSfiles/blob_map_June.nc')
cube = cube[0:200, : , :]
outmax = np.zeros_like(cube[0,:,:])
outmin = np.zeros_like(cube[0,:,:])

cvals = cube.values.flatten()

y ,x = np.meshgrid(range(len(cube.lat)), range(len(cube.lon)))

cube_flat = cvals[np.isfinite(cvals)]
scales = np.unique(cube_flat)

weightso = np.ones_like(cube_flat) / float(len(cube_flat))

histo, ho = np.histogram(cube_flat, bins=np.arange(10,201,10), weights=weightso, range=(10,200))

mid = ho[1::]-5
#plt.imshow(cube.values[5,:,:])
for yy, xx in zip(y.flatten(),x.flatten()):

    dat = cube.values[:, yy, xx]
    dat = dat[np.isfinite(dat)]

    weights = np.ones_like(dat) / float(len(dat))
    hist, h = np.histogram(dat, bins=np.arange(10,201,10), weights=weights, range=(10,200))

    hplot = hist - histo

    try:
       pos =  np.nanargmax(hplot)
       ismax = mid[pos]
       if (hplot[pos] <= 0) or (hist[pos] == 0) or (np.abs(hplot[pos]) < 3*np.std(hist)):
           ismax = np.nan

    except ValueError:
       ismax = np.nan

    try:
        pos = np.nanargmin(hplot)
        ismin = mid[pos]
        if (hplot[pos] >= 0) or (hist[pos] == 0) or (np.abs(hplot[pos]) < 3*np.std(hist)):
            ismin = np.nan

    except ValueError:
        ismin = np.nan

    # if np.nansum(dat) > 10000:
    #     ipdb.set_trace()
    outmax[yy,xx] = ismax
    outmin[yy, xx] = ismin

damax = xr.DataArray(outmax, coords=[cube.lat, cube.lon, ], dims=['lat', 'lon'])
damin = xr.DataArray(outmin, coords=[cube.lat, cube.lon, ], dims=['lat', 'lon'])

u_plot.quick_map(damax, save = '/users/global/cornkle/C_paper/wavelet/figs/paper/map_max_June.png')

u_plot.quick_map(damin, save = '/users/global/cornkle/C_paper/wavelet/figs/paper/map_min_June.png')

plt.close('all')
