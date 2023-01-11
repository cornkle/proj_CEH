import numpy as np
import xarray as xr
import glob
from wavelet import util
import os
from utils import u_grid, u_interpolate as u_int, constants as cnst
import ipdb


#dummy = xr.open_dataset(glob.glob('/prj/vera/cores_bigDomain/*.nc')[0])
#grid = dummy.salem.grid

filepath = glob.glob('/media/ck/Elements/Africa/WestAfrica/CP4/CP4fut/EF/*.nc')

#dlst = xr.open_dataset(filepath[0])
#inds, weights, shape = u_int.interpolation_weights_grid(dlst['lon'].values, dlst['lat'].values, grid)

for f in filepath:

    outpath = f.replace('EF', 'EF_wav')

    if os.path.isfile(outpath):
        continue


    ds = xr.Dataset()
    dat = xr.open_dataarray(f)
    dat = dat.sel(time=dat.time['time.hour']==12).squeeze()

    lsw = dat.copy()
   # ipdb.set_trace()
    lsw.values[(lsw.values<0.05) | (lsw.values>=1)] = 0
    pos = np.where((lsw.values<0.05) | (lsw.values>=1))

    dic = util.applyHat_pure(lsw, dataset='NOWCAST')

    wav = xr.DataArray(np.array(np.round(dic['coeffs'], 1) * 10).astype(np.int16),
                      coords={'scales': dic['scales'], 'latitude':dat.latitude,
                                         'longitude':dat.longitude},
                      dims=['scales','latitude','longitude'])
    ds['wavelet'] = wav

    # for sliced in ds['wavelet'].values:
    #     sliced[np.isnan(lsta)] = np.nan
    ds['wavelet'].values[:, pos[0], pos[1]] = 0

    comp = dict(zlib=True, complevel=5)
    enc = {var: comp for var in ds.data_vars}

    ds.to_netcdf(path=outpath, mode='w', encoding=enc, format='NETCDF4')

    del lsw
