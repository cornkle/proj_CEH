import xarray as xr
import glob
import os
import itertools


#veg = '/users/global/cornkle/w2018_bamba/qrparm.cci.4km.nc'

### 2d vars , xmh*.pc*.nc files
folder = '/scratch/ssf/'
out =  '/users/global/cornkle/w2018_bamba/linked/small/'

atmo_bstream = '.pb'
flux_cstream = '.pc'
current = 'xmhkga'
past = 'xmhkha'
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
    'STASH_m01s03i225': 'u10',
    'STASH_m01s30i203': 'w_pl',
    'STASH_m01s30i201': 'u_pl',
    'STASH_m01s30i202': 'v_pl',
    'STASH_m01s30i204': 'T_pl',
    'STASH_m01s30i205': 'q_pl',
    'STASH_m01s30i206': 'rh_pl',
    'STASH_m01s30i208': 'omega_pl'
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
         'w_pl': 'm s-1',
         'u_pl': 'm s-1',
         'v_pl': 'm s-1',
         'T_pl': 'K',
         'q_pl': 'kg kg-1',
         'rh_pl': '%',
         'omega_pl': 'Pa s-1'

         }

box = [150, 900, 93, 400]
flist = []

for tx, sx in itertools.product(timex, streamx):
    flist.extend(glob.glob(folder +tx + sx+ '*.nc'))

for f in flist:

    fname = f.split(os.sep)[-1]
    varsdat = xr.open_dataset(f)

    time = varsdat.TH1_MN.values

    if atmo_bstream in fname:
        varsdat = varsdat[list(dict_bstream.keys())]
        varsdat = varsdat.isel(grid_longitude_t=slice(box[0], box[1]), grid_longitude_uv=slice(box[0], box[1]),
                               grid_latitude_t=slice(box[2], box[3]), grid_latitude_uv=slice(box[2], box[3]))
        varsdat.rename(dict_bstream, inplace=True)
    if flux_cstream in fname:
        varsdat = varsdat[list(dict_cstream.keys())]
        varsdat = varsdat.isel(grid_longitude_t=slice(box[0], box[1]), grid_latitude_t=slice(box[2], box[3]))
        varsdat.rename(dict_cstream, inplace=True)

    vkeys = varsdat.keys()
    for key in vkeys:
        try:
            unit = units[key]
        except KeyError:
            continue
        varsdat[key].attrs['unit'] = unit

    if current in fname:
        fname = fname.replace(current, 'current')
    if past in fname:
        fname = fname.replace(past, 'past')

    varsdat.to_netcdf(out+fname)
    varsdat.close()



