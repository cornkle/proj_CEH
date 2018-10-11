import xarray as xr
import glob
import os
import itertools
import salem
import numpy as np
import pdb
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy
from utils import u_interpolate

#veg = '/users/global/cornkle/w2018_bamba/qrparm.cci.4km.nc'

### 2d vars , xmh*.pc*.nc files
vera_folder = '/users/global/cornkle/figs/VERA/bamba/w2018_bamba/'  #'/home/users/cornkle/linked_vera/' #
ancils = vera_folder + 'ancils/'   # where to find veg_past, veg_current and topography
#in_folder = vera_folder + 'intest/'
out =  vera_folder + 'outtest/' #'/work/scratch/cornkle/vera_out/' ##
dummy_grid = vera_folder + 'dummy_grid.nc'   # this dummy grid is the big vera grid (ancils have no coord info)

atmo_bstream = '.pb'
flux_cstream = '.pc'
current = ['xmhkg', 'xmhkk', 'xmhkm', 'xmhko']
past = ['xmhkj', 'xmhkl', 'xmhkn', 'xmhkp']
timex = current + past
streamx = [atmo_bstream, flux_cstream]

dict_cstream = {
    'STASH_m01s01i201_2': 'sw_net',
    'STASH_m01s01i235_2': 'sw_in',
    'STASH_m01s02i201_2': 'lw_net',
    'STASH_m01s02i207_2': 'lw_in',
    'STASH_m01s03i217_2': 'SH',
    'STASH_m01s03i234_2': 'LH',
    #'STASH_m01s05i205_2': 'ConvRain',
    'STASH_m01s05i216_2': 'TotalRain',
    #'STASH_m01s04i203_2': 'LargeScaleRain',

}

dict_bstream = {
    'STASH_m01s00i024': 'lst',
    'STASH_m01s03i237': 'q2',
    'STASH_m01s03i236': 'T2',
    'STASH_m01s16i222': 'slp',
    'STASH_m01s03i226_2': 'v10',
    'STASH_m01s03i225_2': 'u10',
    'STASH_m01s30i203': 'w_pl',
    'STASH_m01s30i201': 'u_pl',
    'STASH_m01s30i202': 'v_pl',
    'STASH_m01s30i204': 'T_pl',
    'STASH_m01s30i205': 'q_pl',
    'STASH_m01s30i206': 'rh_pl',
    'STASH_m01s30i208': 'omega_pl',
}

units = {'sw_net': 'W m-2',
         'sw_in': 'W m-2',
         'lw_net': 'W m-2',
         'lw_in': 'W m-2',
         'SH': 'W m-2',
         'LH': 'W m-2',
         'ConvRain': 'mm h-1', #'kg m-2s-1',
         'TotalRain': 'mm h-1',
         'LargeScaleRain': 'mm h-1',
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

box = [150, 900, 20, 400]
wind_box = [box[0]-1,box[1],box[2]-1,box[3]]



def unstagger_metum(var):

    if var.ndim == 1:
        unstag_ndim = var[0:-1] + (var[1::] - var[0:-1]) / 2

    if var.ndim == 2:
        unstag = var[:, 0:-1]  + (var[:, 1::] - var[:, 0:-1]) / 2 # xdir
        unstag_ndim = unstag[0:-1, :]  + (unstag[1::, :] - unstag[0:-1, :]) / 2 # ydir

    if var.ndim == 3:
        unstag = var[:, :, 0:-1]  + (var[:, :, 1::] - var[:, :, 0:-1]) / 2 # xdir
        unstag_ndim = unstag[:, 0:-1, :]  + (unstag[:, 1::, :] - unstag[:, 0:-1, :]) / 2 # ydir

    if var.ndim == 4:
        unstag = var[:, :, :, 0:-1]  + (var[:, :, :, 1::] - var[:, :, :, 0:-1]) / 2 # xdir
        unstag_ndim = unstag[:, :, 0:-1, :]  + (unstag[:, :, 1::, :] - unstag[:, :, 0:-1, :]) / 2 # ydir

    return unstag_ndim

def create_regularGrid(lat2d, lon2d):

    mid_start_lon = lon2d[:,0].mean()
    mid_end_lon = lon2d[:,-1].mean()
    mid_start_lat = lat2d[0,:].mean()
    mid_end_lat = lat2d[-1,:].mean()

    #dist = np.round(np.float(np.mean(lats[1::] - lats[0:-1])), decimals=4)
    pdb.set_trace()
    lat_reg = np.linspace(mid_start_lat, mid_end_lat, lon2d.shape[0])
    lon_reg = np.linspace(mid_start_lon, mid_end_lon, lon2d.shape[1])

    return lon_reg, lat_reg


def create_da(data, time, lat1d, lon1d, lat2d, lon2d, plevels=None):
    pdb.set_trace()
    if plevels is not None:

        da = xr.DataArray(data,
                          coords={'time': time, 'false_latitude': lat1d,
                                  'false_longitude': lon1d,
                                  'true_latitude': (
                                      ['false_latitude', 'false_longitude'], lat2d),
                                  'true_longitude': (
                                      ['false_latitude', 'false_longitude'], lon2d),
                                  'pressure_levels': plevels},
                          dims=['time', 'pressure_levels', 'false_latitude', 'false_longitude'])
    else:

        da = xr.DataArray(data,
                          coords={'time': time, 'false_latitude': lat1d,
                                  'false_longitude': lon1d,
                                  'true_latitude': (
                                      ['false_latitude', 'false_longitude'], lat2d),
                                  'true_longitude': (
                                      ['false_latitude', 'false_longitude'], lon2d)},
                          dims=['time', 'false_latitude', 'false_longitude'])
    return da





def run_datafiles():
    print('start script')
    flist = []

    for tx, sx in itertools.product(timex, streamx):

        run_name = tx+'a'
        time_dir = ''
        if tx in current:
            time_dir = 'veg_current'
        if tx in past:
            time_dir = 'veg_1950'
        print(vera_folder + time_dir + os.sep +tx + os.sep + run_name + sx)
        flist.extend(glob.glob(vera_folder + time_dir + os.sep +tx + os.sep + run_name + sx+ '*.nc'))
    print('Created file list')
    #
    # dummy = xr.open_dataarray(dummy_grid)
    # lats = dummy.latitude_t.values
    # lons = dummy.longitude_t.values
    #
    # start_lon =

    for f in flist:
        
        print('Doing ',f)
        fname = f.split(os.sep)[-1]

        run_name = fname.split('.')[0][0:-1]

        varsdat = salem.open_metum_dataset(f)
        newdat = xr.Dataset(attrs=varsdat.attrs)
        rot_da = varsdat['rotated_latitude_longitude']
        time = varsdat.TH1_MN.values

        if atmo_bstream in fname:
            varsdat = varsdat[list(dict_bstream.keys())]
            varsdat = varsdat.isel(grid_longitude_t=slice(box[0], box[1]), grid_longitude_uv=slice(wind_box[0], wind_box[1]),
                                   grid_latitude_t=slice(box[2], box[3]), grid_latitude_uv=slice(wind_box[2], wind_box[3]))
            #varsdat.rename(dict_bstream, inplace=True)

            for var in varsdat.data_vars:

                darr = varsdat[var]
                data = darr.values.squeeze()
                plevels=None

                if ('P_ECMWF' in darr.coords) | ('height_10m' in darr.coords):

                    lat_t = varsdat.grid_latitude_t.values
                    lon_t = varsdat.grid_longitude_t.values

                    lat_uv = varsdat.grid_latitude_uv.values
                    lon_uv = varsdat.grid_longitude_uv.values

                    lat_ustag = unstagger_metum(lat_uv)
                    lon_ustag = unstagger_metum(lon_uv)


                    np.testing.assert_allclose(lat_t, lat_ustag, atol=1e-5)
                    np.testing.assert_allclose(lon_t, lon_ustag, atol=1e-5)

                    data = unstagger_metum(data)

                if ('P_ECMWF' in darr.coords):

                        data = data[:,1:-1,:,:]
                        nans = np.where(data==0)
                        data[nans] = np.nan
                        plevels=varsdat.P_ECMWF.values[1:-1]

                da = create_da(data, time, varsdat.grid_latitude_t.values, varsdat.grid_longitude_t.values,
                               varsdat.latitude_t.values, varsdat.longitude_t.values, plevels=plevels)

                newdat[dict_bstream[var]] = da


        if flux_cstream in fname:

            varsdat = varsdat[list(dict_cstream.keys())]
            varsdat = varsdat.isel(grid_longitude_t=slice(box[0], box[1]), grid_latitude_t=slice(box[2], box[3]))
            #varsdat.rename(dict_cstream, inplace=True)
            c = 1
            for var in varsdat.data_vars:
                if 'Rain' in var:
                    c = 3600
                data = varsdat[var].values*c
                da = create_da(data, time, varsdat.grid_latitude_t.values, varsdat.grid_longitude_t.values,
                               varsdat.latitude_t.values, varsdat.longitude_t.values)

                newdat[dict_cstream[var]] = da

        vkeys = newdat.keys()
        for key in vkeys:
            try:
                unit = units[key]
            except KeyError:
                continue
            newdat[key].attrs['unit'] = unit

        newdat['rotated_latitude_longitude'] = rot_da

        if run_name in current:
            fname = fname.replace(run_name, 'current_'+run_name)
        if run_name in past:
            fname = fname.replace(run_name, 'past_'+run_name)

        comp = dict(zlib=True, complevel=5)
        encoding = {var: comp for var in newdat.data_vars}
        newdat.to_netcdf(path=out+fname, mode='w', encoding=encoding, format='NETCDF4')

        varsdat.close()
        newdat.close()


def create_ancils():


    dummy = xr.open_dataarray(dummy_grid)

    ds = xr.Dataset(attrs=dummy.attrs)
    dummy = dummy.isel(grid_longitude_t=slice(box[0], box[1]), grid_latitude_t=slice(box[2], box[3]))

    files = glob.glob(ancils+'*.nc')

    for f in files:
        varsdat = xr.open_dataset(f, decode_times=False)

        if 'pseudo' in varsdat.keys():

            varsdat = varsdat.isel(rlon=slice(box[0], box[1]), rlat=slice(box[2], box[3]))

            data = varsdat['field1391'].values[0, 0, :,:].squeeze() # time, plant type, y, x

            if 'past' in f:
                var = 'veg_past'
            elif 'current' in f:
                var = 'veg_current'
            else:
                print('Ancils not found')
                return
            ds[var] = xr.DataArray(data, coords={'false_latitude': dummy.grid_latitude_t.values,
                                          'false_longitude': dummy.grid_longitude_t.values,
                                          'true_latitude': (
                                          ['false_latitude', 'false_longitude'], dummy.latitude_t.values),
                                          'true_longitude': (
                                          ['false_latitude', 'false_longitude'], dummy.longitude_t.values)},
                                  dims=['false_latitude', 'false_longitude'])

        if 'ht' in varsdat.keys():

            varsdat = varsdat.isel(rlon=slice(box[0], box[1]), rlat=slice(box[2], box[3]))
            data = varsdat['ht'].values[0, 0, :, :].squeeze()  # time, plant type, y, x

            ds['topo'] = xr.DataArray(data, coords={'false_latitude': dummy.grid_latitude_t.values,
                                                 'false_longitude': dummy.grid_longitude_t.values,
                                                 'true_latitude': (
                                                     ['false_latitude', 'false_longitude'], dummy.latitude_t.values),
                                                 'true_longitude': (
                                                     ['false_latitude', 'false_longitude'], dummy.longitude_t.values)},
                                   dims=['false_latitude', 'false_longitude'])

    ds.to_netcdf(out+'ancils_vera.nc')

# xarray.DataArray(var.data, var.coords, var.dims, var.attrs)



