import numpy as np
import xarray as xr
from utils import constants as cnst
import os
import sys


cmip_models = ['MPI-ESM1-2-LR', 'CESM2', 'IPSL-CM6A-LR', 'CNRM-CM6-1', 'EC-Earth3']
experiments = {'amip-lfmip-rmLC_hist' : (1980,2007), 'amip-lfmip-rmLC_fut' : (2080,2101), 'historical' : (1980,2007), 'ssp585': (2080,2101), 'amip-lfmip-pdLC_hist' : (1980,2007), 'amip-lfmip-pdLC_fut' : (2080,2101)}
mtag = sys.argv[1]
VARS=['tas', 'ua', 'va', 'prw', 'hfls', 'hfss', 'pr', 'psl', 'tasmax', 'tasmin'] # tasmax tasmin, missing in CESM2 historical
#MODEL = cmip_models[1]
MPI = {}
for MODEL in cmip_models:
    print(MODEL)
    for vv in VARS:
        print(vv)
        for exp in experiments:
            print(exp)
            if "rmLC" in exp:
                exp_read = 'amip-lfmip-rmLC'
            elif "pdLC" in exp:
                exp_read = 'amip-lfmip-pdLC'
            else:
                exp_read = exp
            path = cnst.lmcs_drive+'CMIP6/LS3MIP/'+exp_read+'/'+vv+'_*_'+MODEL+'*.nc'
            opath = cnst.lmcs_drive+'CMIP6/LS3MIP/slices_30y/'+exp_read+'/'
            outpath = opath + vv+'_30y_'+mtag+'_'+str(experiments[exp][0])+'-'+str(experiments[exp][1])+'_'+MODEL+'.nc'

            if os.path.isfile(outpath):
                print('File exists, continue')
                continue
            try:
                ds = xr.open_mfdataset(path, use_cftime=True)
            except:
                print('FIle not found, continue', path)
                continue
            ds = ds.sel(time=(ds['time.year']>=experiments[exp][0]) & (ds['time.year']<=experiments[exp][1]))
            #ipdb.set_trace()
            if mtag == 'DJF':
                ds = ds.sel(time=(ds['time.month']>=12) | (ds['time.month']<=2))
            if mtag == 'JJA':
               ds = ds.sel(time=(ds['time.month']>=6) & (ds['time.month']<=8))
            if mtag == 'MAM':
               ds = ds.sel(time=(ds['time.month']>=3) & (ds['time.month']<=5))
            if mtag == 'SON':
               ds = ds.sel(time=(ds['time.month']>=9) & (ds['time.month']<=11))
            uni = np.unique(ds['time.year'])
            print(exp, 'year slice', uni[0],uni[-1])
            print('months', np.unique(ds['time.month']))
            ds = ds.mean('time')
            comp = dict(zlib=True, complevel=5)
            encoding = {var: comp for var in ds.data_vars}


            ds.to_netcdf(outpath, mode='w', encoding=encoding, format='NETCDF4')
            print('Written', outpath)
            del ds

