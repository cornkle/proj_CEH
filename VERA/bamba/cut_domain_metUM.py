import xarray as xr
import pdb
import matplotlib.pyplot as plt
import numpy as np
from utils import constants as cnst, u_arrays, u_grid
import salem



veg = '/users/global/cornkle/w2018_bamba/qrparm.cci.4km.nc'
topo = cnst.WA_TOPO_1MIN

filepart = '/scratch/ssf/xmhkga.pb20140407.nc'
keep_vars2d = ['STASH_m01s00i024', 'STASH_m01s03i237', 'STASH_m01s03i236', 'STASH_m01s16i222', 'STASH_m01s03i226' , 'STASH_m01s03i225'] # LST, spec hum 1.5m, T 1.5m, p at sl,10m v wind, 10m u wind
keep_vars3d = ['STASH_m01s30i203', 'STASH_m01s30i201', 'STASH_m01s30i202', 'STASH_m01s30i204', 'STASH_m01s30i205', 'STASH_m01s30i206', 'STASH_m01s30i208' ] # u on P, v on P, T on P, spec H on P, rel H on P, omega on P

vars2d_names = ['lst', 'q2', 'T2', 'slp', 'v10', 'u10']
vars3d_names = ['w_pl', 'u_pl', 'v_pl', 'T_pl', 'q_pl', 'rh_pl', 'omega_pl']
joined_names = vars2d_names+vars3d_names
name_dict = {}

for id, n in enumerate(keep_vars2d+keep_vars3d):

    name_dict[n] = joined_names[id]

box = [325, 475, 93, 235 ]  # x1, x2, y1, y2


vegdat = xr.open_dataarray(veg)
varsdat = xr.open_mfdataset(filepart, concat_dim='time')


#varsdat = varsdat.sel(time=slice('2014-04-06', '2014-04-07'))

vegdat = vegdat.isel(longitude=slice(box[0], box[1]), latitude=slice(box[2], box[3]))
varsdat = varsdat.isel(grid_longitude_t=slice(box[0], box[1]), grid_latitude_t=slice(box[2], box[3]))
varsdat = varsdat.isel(grid_longitude_uv=slice(box[0], box[1]), grid_latitude_uv=slice(box[2], box[3]))
varsdat = varsdat[keep_vars2d + keep_vars3d]
varsdat = varsdat.sel(P_ECMWF=slice(950, 550))


varsdat.rename(name_dict, inplace=True)

# vegdat.latitude = varsdat.latitude_t.values.mean(axis=1)
# vegdat.longitude = varsdat.longitude_t.values.mean(axis=0)
varsdat.coords['grid_latitude_t'] = varsdat.latitude_t.values.mean(axis=1)
varsdat.coords['grid_longitude_t'] = varsdat.longitude_t.values.mean(axis=0)
varsdat.coords['grid_latitude_uv'] = varsdat.latitude_uv.values.mean(axis=1)
varsdat.coords['grid_longitude_uv'] = varsdat.longitude_uv.values.mean(axis=0)

varsdat['forest_frac'] = (('pseudo_dim', 'grid_latitude_t', 'grid_longitude_t'), vegdat.values[0,:,:][None,...])

factor = (1000/varsdat.P_ECMWF.values)**0.286
test = np.repeat(factor.T, 142*150, axis=0)
test.shape = (10,142,150)

varsdat['theta_pl'] = (('TH1', 'P_ECMWF', 'grid_latitude_t', 'grid_longitude_t'), varsdat['T_pl']*test)

# plt.figure()
# #varsdat.q2[0,:,:].squeeze().plot.contourf()
# #pdb.set_trace()
# varsdat.theta_pl[0,:,40, :].squeeze().plot.contourf(cmap='jet')
#
# plt.figure()
# #varsdat.q2[0,:,:].squeeze().plot.contourf()
# #pdb.set_trace()
# varsdat.T_pl[0,:,40, :].squeeze().plot.contourf(cmap='jet')

varsdat.to_netcdf('/users/global/cornkle/w2018_bamba/mini_forest/pb20140407_forest.nc')

