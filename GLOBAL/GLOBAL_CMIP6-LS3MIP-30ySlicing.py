import numpy as np
import xarray as xr
from utils import constants as cnst
import os


cmip_models = ['MPI-ESM1-2-LR', 'CESM2', 'IPSL-CM6A-LR', 'CNRM-CM6-1']#, 'EC-Earth3']
experiments = {'historical' : (1980,2007), 'ssp585': (2080,2101), 'amip-lfmip-pdLC_hist' : (1980,2007), 'amip-lfmip-pdLC_fut' : (2080,2101)}

VARS=['tas', 'ua', 'va', 'prw', 'hfls', 'hfss', 'pr', 'psl'] # tasmax tasmin, missing in CESM2 historical
#MODEL = cmip_models[1]
MPI = {}
for MODEL in cmip_models:
    print(MODEL)
    for vv in VARS:
        print(vv)
        for exp in experiments:
            print(exp)
            if "amip" in exp:
                exp_read = 'amip-lfmip-pdLC'
            else:
                exp_read = exp
            path = cnst.lmcs_drive+'CMIP6/LS3MIP/'+exp_read+'/'+vv+'_*_'+MODEL+'*.nc'
            outpath = path.replace(vv+'_*_'+MODEL+'*', '30y/'+vv+'_30y_'+str(experiments[exp][0])+'-'+str(experiments[exp][1])+'_'+MODEL)

            if os.path.isfile(outpath):
                print('File exists, continue')
                continue

            ds = xr.open_mfdataset(path, use_cftime=True)
            ds = ds.sel(time=(ds['time.year']>=experiments[exp][0]) & (ds['time.year']<=experiments[exp][1]))
            #ipdb.set_trace()
            uni = np.unique(ds['time.year'])
            print(exp, 'year slice', uni[0],uni[-1])
            ds = ds.mean('time')
            comp = dict(zlib=True, complevel=5)
            encoding = {var: comp for var in ds.data_vars}


            ds.to_netcdf(outpath, mode='w', encoding=encoding, format='NETCDF4')
            del ds

