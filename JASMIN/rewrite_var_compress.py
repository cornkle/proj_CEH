import numpy as np
import xarray as xr
import os
import glob
import ipdb

def rewrite(file, varname):
    
    with xr.open_dataset(file) as da:
        da = xr.open_dataset(file)
        ipdb.set_trace()
        #da = xr.open_dataarray(file)
        da.name = varname
        #da.values = (np.round(da.values-273.15, decimals=2)*100).astype(int)
        da_write = da.copy()
    encoding = {da_write.name: {'complevel': 5, 'zlib': True}}
    filen = file.replace('.nc', '_compressed.nc')
    if not os.path.isfile(filen):
        da_write.to_netcdf(filen, format='NETCDF4', encoding=encoding)
        os.remove(file)


files = glob.glob('/home/users/cornkle/CP4home/CP4fut/t_pl/*.nc')
for f in files:
    rewrite(f)
