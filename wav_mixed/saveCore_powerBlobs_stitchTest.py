# -*- coding: utf-8 -*-


from eod import msg, mfg
from utils import u_grid, u_interpolate as u_int
from utils import constants as cnst
import matplotlib.pyplot as plt
import xarray as xr
import pickle


m = mfg.ReadMfg(cnst.network_data +'/data/OBS/MFG_JJAS/', y1=1992, y2=1993, months=[6,7])
files  = m.fpath

gridll = pkl.load( open (cnst.network_data + 'data/OBS/saves/VERA_msg_latlon_18W12E_1N17N.p', 'rb'))

mdic = m.read_data(files[0], llbox=[-25, 20, 2, 25])  #[-14, 2.5, 4, 11.5]

# make salem grid
grid = u_grid.make(gridll['lon'].values, gridll['lat'].values, 5000)
inds, weights, shape = u_int.interpolation_weights_grid(mdic['lon'].values, mdic['lat'].values, grid)
gridd = (inds,weights,shape, grid)
ind = 35
array = mdic['t'][ind,:,:].values

outt = u_int.interpolate_data(array, inds, weights, shape)
lon, lat = grid.ll_coordinates

da = xr.DataArray(outt, coords={'lat': lat[:,0], 'lon':lon[0,:]}, dims=['lat', 'lon'])

da.to_netcdf('/users/global/cornkle/test.nc')