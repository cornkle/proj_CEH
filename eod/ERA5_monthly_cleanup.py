import glob
import xarray as xr
import ipdb

files = glob.glob('/prj/AMMA2050/ERA5/monthly/synoptic/surface/*.nc')

for f in files:
    ds = xr.open_dataset(f)
    ds = ds.where(ds['time.day']==1, drop=True)

    out = f.replace('synoptic', 'synoptic_clean')
    #ipdb.set_trace()
    comp = dict(zlib=True, complevel=5)
    encoding = {var: comp for var in ds.data_vars}
    ds.to_netcdf(out, mode='w', encoding=encoding, format='NETCDF4')