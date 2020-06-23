import numpy as np
import xarray as xr
import glob
from wavelet import util
import os


filepath = glob.glob('/prj/swift/SEVIRI_LST/data_anom_wrt_historic_clim_withmask/*/*.nc')
for ff in filepath:

    outpath = ff.replace('data_anom_wrt_historic_clim_withmask', 'data_anom_wrt_historic_clim_withmask_wavelet')
    outpath = outpath.replace('withHistClim', 'withHistClim_wavelet')

    if os.path.isfile(outpath):
        continue

    dat = xr.open_dataarray(ff)

    dat.values[np.isnan(dat.values)] = 0
    dat.values[dat.values > 1000] = 0
    dic = util.applyHat_pure(dat, dataset='NOWCAST')

    da = xr.DataArray(np.array(np.round(dic['coeffs'], 2) * 100).astype(int),
                      coords={'scales': dic['scales'], 'phony_dim_0': dat.phony_dim_0, 'phony_dim_1': dat.phony_dim_1},
                      dims=['scales', 'phony_dim_0', 'phony_dim_1'])
    da.name = 'wavelet_coeff'

    if not os.path.isdir(os.path.dirname(outpath)):
        os.makedirs(os.path.dirname(outpath))

    comp = dict(zlib=True, complevel=5)
    encoding = {'wavelet_coeff': comp}

    try:
        da.to_netcdf(outpath, format='NETCDF4', encoding=encoding)
    except OSError:
        print('Did not find ' + outpath)
        print('Out directory not found')
