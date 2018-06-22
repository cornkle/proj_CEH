import xarray as xr
import pdb
import numpy as np
import glob
import os
from utils import u_interpolate as uint


#veg = '/users/global/cornkle/w2018_bamba/qrparm.cci.4km.nc'

### 2d vars , xmh*.pc*.nc files
folder = '/scratch/ssf/'
out = '/users/global/cornkle/w2018_bamba/clean_files_raw/'

atmo_bstream = '.pb'
flux_cstream = '.pc'
current = 'current'
past = 'past'
timex = [current, past]
streamx = [atmo_bstream, flux_cstream]

dict_cstream = {
    'STASH_m01s01i201_2': 'sw_net',
    'STASH_m01s01i235_2': 'sw_in',
    'STASH_m01s02i201_2': 'lw_net',
    'STASH_m01s02i207_2': 'lw_in',
    'STASH_m01s03i217_2': 'SH',
    'STASH_m01s03i234_2': 'LH',
    'STASH_m01s05i205_2': 'ConvRain',
    'STASH_m01s05i216_2': 'TotalRain',
    'STASH_m01s04i203_2': 'LargeScaleRain',
}

dict_bstream = {
    'STASH_m01s00i024': 'lst',
    'STASH_m01s03i237': 'q2',
    'STASH_m01s03i236': 'T2',
    'STASH_m01s16i222': 'slp',
    'STASH_m01s03i226': 'v10',
    'STASH_m01s03i225': 'u10'
}

units = {'sw_net': 'W m-2',
         'sw_in': 'W m-2',
         'lw_net': 'W m-2',
         'lw_in': 'W m-2',
         'SH': 'W m-2',
         'LH': 'W m-2',
         'ConvRain': 'kg m-2s-1',
         'TotalRain': 'kg m-2s-1',
         'LargeScaleRain': 'kg m-2s-1',
         'lst' : 'K' ,
         'q2' : 'kg kg-1',
         'T2' : 'K',
         'slp' : 'Pa',
         'v10' : 'm s-1',
         'u10' : 'm s-1',

         }

box = [30, 730, 10, 75 ]  # x1, x2, y1, y2 Tai park box

flist = glob.glob(folder+'current*.nc')

dummy = xr.open_dataset(flist[0])
#dummy = dummy.isel(grid_longitude_t=slice(box[0], box[1]), grid_latitude_t=slice(box[2], box[3]))

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
inds, weights, shape = uint.interpolation_weights(lons, lats, lon_regular, lat_regular)

for tt in timex:
    flist = glob.glob(folder + tt+'*.nc')
    dates = []
    for f in flist:
        dates.append(f[-14:-3])
    dates = np.unique(dates)

    for d in dates:

        flist = []
        ds = xr.Dataset(coords={'lat': lat_regular, 'lon': lon_regular})
        for xs in streamx:
            flist.extend(glob.glob(folder + tt + xs + d+ '.nc'))

        for f in flist:

            fname = f.split(os.sep)[-1]
            varsdat = xr.open_dataset(f)
            #varsdat = varsdat.sel(P_ECMWF=slice(975,600))
            time = varsdat.TH1_MN.values
            if 'time' not in ds.keys():
                ds['time'] = time

            # if atmo_bstream in fname:
            #     varsdat = varsdat.isel(grid_longitude_t=slice(box[0], box[1]), grid_longitude_uv=slice(box[0], box[1]),
            #                            grid_latitude_t=slice(box[2], box[3]), grid_latitude_uv=slice(box[2], box[3]))
            #     # varsdat = varsdat[list(dict_bstream.keys())]
            #     # varsdat.rename(dict_bstream, inplace=True)
            # if flux_cstream in fname:
            #     varsdat = varsdat.isel(grid_longitude_t=slice(box[0], box[1]), grid_latitude_t=slice(box[2], box[3]))
                # varsdat = varsdat[list(dict_cstream.keys())]
                # varsdat.rename(dict_cstream, inplace=True)

            vkeys = varsdat.keys()
            for key in vkeys:

                # if ('longitude' in key) | ('latitude' in key) | ('TH1' in key) :
                #     continue
                # print('Doing variable', key)
                # data = varsdat[key].values
                # regridded = uint.interpolate_data(data, inds, weights, shape)
                unit = units[key]
                # ds[key] = (('time', 'lat', 'lon'), regridded)
                varsdat[key].attrs['unit'] = unit

            if current in fname:
                fname = fname.replace(current, 'current')
            if past in fname:
                fname = fname.replace(past, 'past')

            varsdat.to_netcdf(out+fname)
            varsdat.close()



