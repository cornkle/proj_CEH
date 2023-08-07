import xarray as xr
import pandas as pd
import numpy as np
import datetime
import multiprocessing
import glob
import os
from utils import constants as cnst
import pdb

def run(var):

    in_path = '/home/users/cornkle/runscript/out_EDW/'+var
    fut_path ='/home/users/cornkle/runscript/fut_out_EDW/'+var

    out = '/home/users/cornkle/CP4home/EDW_files/'
    
    for ip in [in_path, fut_path]:    
        
        ds = xr.open_mfdataset(ip+'/*.nc')
        ff = glob.glob(ip+'/*.nc')[0]  
        fbase = os.path.basename(ff)[0:-20]
        fname = fbase+'mean.nc'
        if os.path.isfile(out+fname):
           print(out+fname, 'File exists, next!')
           continue
       
        ds = ds.mean('time')
        comp = dict(zlib=True, complevel=5)
        encoding = {var: comp for var in ds.data_vars}

        ds.to_netcdf(out+fname, mode='w', encoding=encoding,
                    format='NETCDF4')

        print('Saved', out+fname)
        del ds








