import xarray as xr
import pandas as pd
import numpy as np
import datetime
import multiprocessing
import glob
import os
from utils import constants as cnst
import pdb

def run():

    in_path =  '/home/users/cornkle/impala/shared/CP4A/ncfiles/4km/ANCILS/'

    out = '/home/users/cornkle/CP4home/EDW_files/'
    
    for ip in glob.glob(in_path+'*.nc'):    
        
        ds = xr.open_dataset(ip, decode_times=False)
        try:
           ds = ds.sel(latitude=slice(-36,39), longitude=slice(-18+360,52+360))                    # [-18+360, 52+360, -36,39]
        except:
           ds = ds.sel(rlat=slice(-36,39), rlon=slice(-18+360,52+360))
        ff = ip 
        fbase = os.path.basename(ip)
        fname = fbase
        if os.path.isfile(out+fname):
           print(out+fname, 'File exists, next!')
           continue
        try:      
           ds = ds.assign_coords(longitude=ds.longitude.values-360)
        except:
           ds = ds.assign_coords(rlon=ds.rlon.values-360)
        #pdb.set_trace()
        ds = ds.squeeze()
        comp = dict(zlib=True, complevel=5)
        encoding = {var: comp for var in ds.data_vars}

        ds.to_netcdf(out+fname, mode='w', encoding=encoding,
                    format='NETCDF4')

        print('Saved', out+fname)
        del ds








