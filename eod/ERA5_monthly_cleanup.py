import glob
import xarray as xr
import ipdb
from utils import constants as cnst

def synop_day():
    files = glob.glob('/prj/AMMA2050/ERA5/monthly/synoptic/surface/*.nc')

    for f in files:
        ds = xr.open_dataset(f)
        ds = ds.where(ds['time.day']==1, drop=True)

        out = f.replace('synoptic', 'synoptic_clean')
        #ipdb.set_trace()
        comp = dict(zlib=True, complevel=5)
        encoding = {var: comp for var in ds.data_vars}
        ds.to_netcdf(out, mode='w', encoding=encoding, format='NETCDF4')



def rewrite():

    files = glob.glob(cnst.ERA5+'monthly/synop_selfmade/CLIM_2006-2010_new/ERA5_2006-2010_CLIM_*.nc')

    for f in files:
        ds = xr.open_dataset(f, ).load()
        #ds = ds.where(ds['time.day'] == 1, drop=True)
        ds = ds.rename({'lat': 'latitude', 'lon': 'longitude'})
        out = f.replace('.nc', '_rw.nc')
        #ipdb.set_trace()
        # ipdb.set_trace()
        #comp = dict(zlib=True, complevel=5)
        #encoding = {var: comp for var in ds.data_vars}
        ds.to_netcdf(out, mode='w', format='NETCDF4')