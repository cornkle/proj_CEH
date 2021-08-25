import numpy as np
import xarray as xr
import glob
from wavelet import util
import os
from utils import u_grid, u_interpolate as u_int, constants as cnst



dummy = xr.open_dataset(glob.glob('/prj/vera/cores/bigDomain/*.nc')[0])
grid = dummy.salem.grid

filepath = glob.glob('/prj/swift/SEVIRI_LST/data_anom_wrt_historic_clim_withmask/*/*.nc')

dlst = xr.open_dataset(filepath[0])
inds, weights, shape = u_int.interpolation_weights_grid(dlst['lon'].values, dlst['lat'].values, grid)

for f in filepath:

    outpath = f.replace('data_anom_wrt_historic_clim_withmask', 'data_anom_wrt_historic_clim_withmask_wavelet')
    outpath = outpath.replace('withHistClim', 'withHistClim_wavelet')

    if os.path.isfile(outpath):
        continue


    ds = xr.Dataset()
    dat = xr.open_dataset(f)

    try:
        lsta = u_int.interpolate_data(dat['lsta'].values, inds, weights, shape)
    except IndexError:
        print('Interpolation problem, continue')
        continue

    lon, lat = grid.ll_coordinates

    nbslot = u_int.interpolate_data(dat['NbSlot'].values, inds, weights, shape)

    lsta[np.isnan(lsta)] = 0
    dlsta = xr.DataArray((np.round(lsta, 2)*100).astype(np.int16), coords={
        'lat': lat[:, 0],
        'lon': lon[0, :]}, dims=['lat', 'lon'])

    dslot = xr.DataArray(np.round(nbslot,0).astype(np.int8), coords={
        'lat': lat[:, 0],
        'lon': lon[0, :]}, dims=['lat', 'lon'])

    ds['lsta'] = dlsta
    ds['NbSlot'] = dslot

    lsw = lsta.copy()

    lsw[lsw > 1000] = 0
    dic = util.applyHat_pure(lsw, dataset='NOWCAST')

    wav = xr.DataArray(np.array(np.round(dic['coeffs'], 1) * 10).astype(np.int16),
                      coords={'scales': dic['scales'], 'lat': lat[:,0],
                                         'lon': lon[0,:]},
                      dims=['scales', 'lat', 'lon'])
    ds['wavelet'] = wav

    # for sliced in ds['wavelet'].values:
    #     sliced[np.isnan(lsta)] = np.nan
    pos = np.where(lsta==0)
    pos = np.where(lsta==0)
    ds['wavelet'].values[:, pos[0], pos[1]] = 0

    ds = ds.isel(lon=slice(0,1410), lat=slice(0,594))
    comp = dict(zlib=True, complevel=5)
    enc = {var: comp for var in ds.data_vars}
    savefile = f.replace('netcdf', 'netcdf_onCores')
    ds.to_netcdf(path=savefile, mode='w', encoding=enc, format='NETCDF4')

    del lsw