import numpy as np
import xarray as xr
from utils import constants as cnst, u_grid, u_interpolate as u_int
import os
import salem
import glob
import ipdb

cmip_models = ['MPI-ESM1-2-LR', 'CESM2', 'IPSL-CM6A-LR', 'CNRM-CM6-1', 'EC-Earth3']#, 'EC-Earth3']
experiments = {'historical' : (1980,2007), 'ssp585': (2080,2101), 'amip-lfmip-pdLC_hist' : (1980,2007), 'amip-lfmip-pdLC_fut' : (2080,2101)}

VARS=['tas', 'ua', 'va', 'prw', 'hfls', 'hfss', 'pr', 'psl','tasmax', 'tasmin'] # tasmax tasmin, missing in CESM2 historical
print(cnst.lmcs_drive + 'CMIP6/LS3MIP/slices_30y/amip-lfmip-pdLC_hist/tas_30y_*_' + cmip_models[0] + '.nc')
dummy = xr.open_dataset(glob.glob(cnst.lmcs_drive + 'CMIP6/LS3MIP/slices_30y/amip-lfmip-pdLC_hist/tas_30y_*_' + cmip_models[0] + '.nc')[0])['tas']
#grid = dummy.salem.grid
# for cm in cmip_models[1::]:

#

for MODEL in cmip_models:
    print(MODEL)
    for vv in VARS:
        print(vv)
        for exp in experiments:
            print(exp)

            path = cnst.lmcs_drive+'CMIP6/LS3MIP/slices_30y/'+exp+'/'+vv+'_*_'+MODEL+'*.nc'
            oopath = glob.glob(path)
            for opath in oopath: 
                outpath = opath.replace('slices_30y', 'slices_30y_regrid')

                # if os.path.isfile(outpath):
                #     print('File exists, continue')
                #     continue

                ds = xr.open_dataset(opath)
                # inds, weights, shape = u_int.interpolation_weights(ds.lon, ds.lat, dummy.lon, dummy.lat)
                # dint = u_int.interpolate_data(ds, inds, weights, shape)
                ds = dummy.salem.transform(ds[vv])

                comp = dict(zlib=True, complevel=5)
                encoding = {vv: comp}

                ds.to_netcdf(outpath, mode='w', encoding=encoding, format='NETCDF4')
                del ds

