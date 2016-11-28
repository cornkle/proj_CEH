import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from utils import u_arrays as ua
import cartopy
import cartopy.crs as ccrs
import salem

import pyproj
import multiprocessing
from functools import partial


def trmm_map():
    files = ua.locate(".nc", '/users/global/cornkle/TRMMfiles')
    pool = multiprocessing.Pool(processes=4)

    files=files[0:5]

    # make a salem grid
    proj = pyproj.Proj('+proj=merc +lat_0=0. +lon_0=0.')
    # Transform lon, lats to the mercator projection
    x, y = pyproj.transform(salem.wgs84, proj,[-18.5, 35], [2.5,25])
    # take the min and max
    xmax, xmin = np.max(x), np.min(x)
    ymax, ymin = np.max(y), np.min(y)
    # Count the number of pixels
    dx = 5000
    nx, r = divmod(xmax - xmin, dx)
    ny, r = divmod(ymax - ymin, dx)
    # Here one could add + 1 to be sure that the last pixel is always included
    grid = salem.Grid(nxny=(nx, ny), dxdy=(dx, dx), ll_corner=(xmin, ymin), proj=proj)

    lon, lat = grid.ll_coordinates
    xi, yi = grid.ij_coordinates

    trmm_arr = np.zeros_like(lon)


    func = partial(file_loop, grid)

    res = pool.map(func, files)
    pool.close()
    #pool.join()
    res = [x for x in res if x is not None]
    all = [item for sublist in res for item in sublist]

    for id, a in enumerate(all):
       print('Doing', id , 'from', len(all))
       trmm_arr[a[1], a[0]] += 1


    fig = plt.figure(figsize=(9, 7), dpi=400)
    ax = fig.add_subplot(111, projection=ccrs.PlateCarree())

    plt.imshow(trmm_arr, extent=(lon.min(), lon.max(), lat.min(), lat.max()), transform=ccrs.PlateCarree(),
               cmap='hot_r')
    # plt.contourf(lon, lat, po30_arr, transform=ccrs.PlateCarree(), cmap='viridis_r', levels=np.arange(0,3,0.5) )
    ax.coastlines()
    plt.colorbar()
    plt.title('Pixel count precip>30')
    ax.add_feature(cartopy.feature.BORDERS, linestyle='--')
    plt.show()


def file_loop(grid, f):
    tlon = []
    tlat = []
    ff = xr.open_dataset(f)
    print('Doing file: ' + f)

    # if np.sum(ff['t']) == 0:
    #     # print('no t, continue')
    #     return
    loc = np.where(np.isfinite(ff['p'].values))
    tlon.append(ff['lon'].values[loc].tolist())
    tlat.append(ff['lat'].values[loc].tolist())

    tlon = [item for sublist in tlon for item in sublist]
    tlat = [item for sublist in tlat for item in sublist]

    xm, ym = grid.transform(tlon, tlat, crs=salem.wgs84)

    ym = np.round(ym)
    xm = np.round(xm)

    tuple = list(zip(xm,ym))

    ff.close()

    return tuple
