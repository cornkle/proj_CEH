# -*- coding: utf-8 -*-


import numpy as np
import xarray as xr
import ipdb
import os
from utils import u_grid, u_interpolate as u_int
import glob
import CCores.cores as cores
import datetime
import multiprocessing


def interpolate(da):

    data_resolution = 3  # in km
    # make salem grid
    grid = u_grid.make(np.arange(-19, 0), np.arange(4, 20), data_resolution * 1000)
    dlon = da['lon_2d'].squeeze().values.T
    dlat = da['lat_2d'].squeeze().values.T
    inds, weights, shape = u_int.interpolation_weights_grid(dlon, dlat, grid)

    data = da['IR108_BT'].squeeze().values.T
    try:
        orig = u_int.interpolate_data(data, inds, weights, shape)
    except IndexError:
        print('Interpolation problem, continue')
        return
    lon, lat = grid.ll_coordinates

    return {'data' : orig, 'lon' : lon, 'lat' : lat, 'inds': inds, 'weights' : weights, 'shape' : shape, 'grid' : grid}



def run_one(filepath):


    da = xr.open_dataset(filepath, decode_times=False)

    #ipdb.set_trace()

    ff = os.path.basename(filepath)
    year = ff[10:14]
    month = ff[14:16]
    day = ff[16:18]
    hour = ff[19:21]
    minute = ff[21:23]

    date = [datetime.datetime(int(year), int(month), int(day), int(hour), int(minute))]

    interp = interpolate(da)

    wObj = cores.dataset('METEOSAT3K_veraLS')
    wObj.read_img(interp['data'], interp['lon'], interp['lat'], edge_smoothing=False)
    wObj.applyWavelet()

    small_scale = wObj.scaleWeighting(wtype='nflicsv2')
    ss = wObj.to_dataarray(date=date, names=['small_scale', 'tir'])
    dominant = wObj.scaleWeighting(wtype='dominant')
    dom = wObj.to_dataarray(date=date, names=['dom', 'tir'])

    newds = xr.merge([ss,dom])

    comp = dict(zlib=True, complevel=5)
    enc = {var: comp for var in newds.data_vars}

    outfile = filepath.replace('108', 'wavelet')
    outfile = outfile.replace('real_time_data', 'real_time_wavelet')


    if not os.path.exists(os.path.dirname(outfile)):
        os.makedirs(os.path.dirname(outfile))

    if os.path.isfile(outfile):
        print('File exists, continue')
        return

    newds.to_netcdf(path=outfile, mode='w', encoding=enc, format='NETCDF4')
    print('Saved ' + outfile)




def multi(directory):

    files = glob.glob(directory + '*/*/*/*.nc')

    dummy = xr.open_dataset(files[0])
    interp = interpolate(dummy)

    passit = []
    for f in files:
        passit.append((f, interp['inds'], interp['weights'], interp['shape'], interp['grid']))

    pool = multiprocessing.Pool(processes=3)
    res = pool.map(_loop, passit)

    # res=[]
    # for l in passit:
    #
    #     res.append(file_loop(l))

    pool.close()



def _loop(passit):

    filepath, inds, weights, shape, grid = passit

    print('Doing ', filepath)


    da = xr.open_dataset(filepath, decode_times=False)

    #ipdb.set_trace()

    ff = os.path.basename(filepath)
    year = ff[10:14]
    month = ff[14:16]
    day = ff[16:18]
    hour = ff[19:21]
    minute = ff[21:23]

    date = [datetime.datetime(int(year), int(month), int(day), int(hour), int(minute))]

    data = da['IR108_BT'].squeeze().values.T
    try:
        data = u_int.interpolate_data(data, inds, weights, shape)
    except IndexError:
        print('Interpolation problem, continue')
        return
    lon, lat = grid.ll_coordinates


    wObj = cores.dataset('METEOSAT3K_veraLS')
    wObj.read_img(data, lon, lat, edge_smoothing=False)
    wObj.applyWavelet()

    small_scale = wObj.scaleWeighting(wtype='nflicsv2')
    ss = wObj.to_dataarray(date=date, names=['small_scale', 'tir'])
    dominant = wObj.scaleWeighting(wtype='dominant')
    dom = wObj.to_dataarray(date=date, names=['dom', 'tir'])

    newds = xr.merge([ss,dom])

    comp = dict(zlib=True, complevel=5)
    enc = {var: comp for var in newds.data_vars}

    outfile = filepath.replace('108', 'wavelet')
    outfile = outfile.replace('real_time_data', 'real_time_wavelet')

    path = os.path.dirname(outfile)
    if not os.path.exists(path):
        os.makedirs(path)

    if os.path.isfile(outfile):
        print('File exists, continue')
        return

    newds.to_netcdf(path=outfile, mode='w', encoding=enc, format='NETCDF4')
    print('Saved ' + outfile)






