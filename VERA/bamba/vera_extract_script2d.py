import xarray as xr
import pdb
import numpy as np
import glob
import os

from matplotlib.pyplot import grid

from utils import u_interpolate as uint


#veg = '/users/global/cornkle/w2018_bamba/qrparm.cci.4km.nc'

### 2d vars , xmh*.pc*.nc files
folder = '/scratch/ssf/'
out = '/users/global/cornkle/w2018_bamba/regrid_raw/'

atmo_bstream = '.pb'
flux_cstream = '.pc'
current = 'current'
past = 'past'
timex = [current, past]
streamx = [atmo_bstream, flux_cstream]

units_cstream = {'sw_net': 'W m-2',
         'sw_in': 'W m-2',
         'lw_net': 'W m-2',
         'lw_in': 'W m-2',
         'SH': 'W m-2',
         'LH': 'W m-2',
         'ConvRain': 'kg m-2s-1',
         'TotalRain': 'kg m-2s-1',
         'LargeScaleRain': 'kg m-2s-1' }
units_bstream = {
         'lst' : 'K' ,
         'q2' : 'kg kg-1',
         'T2' : 'K',
         'slp' : 'Pa',
         'v10' : 'm s-1',
         'u10' : 'm s-1',
         }

flist = glob.glob(folder+'current*.pb*.nc')

dummy = xr.open_dataset(flist[0])

box = [0, 690, 0, 300 ]
dummy = dummy.isel(grid_latitude_t=slice(box[2],box[3]), grid_longitude_t=slice(box[0],box[1]),
                   grid_latitude_uv=slice(box[2],box[3]), grid_longitude_uv=slice(box[0],box[1]))

lats = dummy.latitude_t[:, 0].values
lons = dummy.longitude_t[0, :].values
latmin = dummy.latitude_t.min().values
latmax = dummy.latitude_t.max().values
lonmin = dummy.longitude_t.min().values
lonmax = dummy.longitude_t.max().values
dist = np.round(np.float(np.mean(lats[1::] - lats[0:-1])), decimals=4)

lat_regular = np.arange(latmin+10*dist, latmax - 10*dist, dist)
lon_regular = np.arange(lonmin+10*dist, lonmax - 10*dist, dist)

# interpolation weights for new grid
inds, weights, shape = uint.interpolation_weights(dummy.longitude_t.values, dummy.latitude_t.values, lon_regular, lat_regular)
inds_uv, weights_uv, shr = uint.interpolation_weights(dummy.longitude_uv.values, dummy.latitude_uv.values, lon_regular, lat_regular)

for tt in timex:
    flist = glob.glob(folder + tt+'*.nc')
    dates = []
    for f in flist:
        strs = f.split('.')
        dates.append(strs[1][2:10])
    dates = np.unique(dates)

    for d in dates:

        flist = []
        ds = xr.Dataset(coords={'lat': lat_regular, 'lon': lon_regular})
        for xs in streamx:
            flist.extend(glob.glob(folder + tt + xs + d+ '.nc'))
            flist.extend(glob.glob(folder + tt + xs + d + '_12.nc'))

        for f in flist:

            fname = f.split(os.sep)[-1]

            if '.pc' in fname:
                varsdat = xr.open_mfdataset(f[0:-6]+'*.nc', concat_dim='TH1_MN')

            else:
                varsdat = xr.open_dataset(f)


            if 'time' not in ds.keys():
                if '.pc' in fname:
                    time = varsdat.TH1_MN.values
                else:
                    time = varsdat.TH1.values
                ds['time'] = time
            units = units_cstream
            if atmo_bstream in fname:
                units = units_bstream
                #varsdat = varsdat.sel(P_ECMWF=slice(975,600))


            if atmo_bstream in fname:
                varsdat = varsdat[list(units.keys())]
                varsdat = varsdat.isel(grid_longitude_t=slice(box[0], box[1]), grid_longitude_uv=slice(box[0], box[1]),
                                       grid_latitude_t=slice(box[2], box[3]), grid_latitude_uv=slice(box[2], box[3]))

            if flux_cstream in fname:
                varsdat = varsdat.isel(grid_longitude_t=slice(box[0], box[1]), grid_latitude_t=slice(box[2], box[3]))


            vkeys = varsdat.keys()
            for key in vkeys:

                if ('longitude' in key) | ('latitude' in key) | ('TH1' in key) | ('height' in key) :
                    continue
                print('Doing variable', key)

                data = varsdat[key].squeeze().values
                if data.ndim == 4:
                    data = np.squeeze(data[0,:,:,:])
                if (key == 'u10') | (key == 'v10'):
                    regridded = uint.interpolate_data(data, inds_uv, weights_uv, shape)
                else:
                    regridded = uint.interpolate_data(data, inds, weights, shape)
                unit = units[key]
                ds[key] = (('time', 'lat', 'lon'), regridded)
                ds[key].attrs['unit'] = unit

            if current in fname:
                fname = fname.replace(current, 'current')
            if past in fname:
                fname = fname.replace(past, 'past')


        comp = dict(zlib=True, complevel=5)
        encoding = {var: comp for var in ds.data_vars}
        ds.to_netcdf(out+tt+'_2d_'+d+'.nc', format='NETCDF4', encoding=encoding)
        ds.close()



