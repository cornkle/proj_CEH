import xarray as xr
import pdb




vpast = '/prj/vera/ancils/afforested_lc/final_veg_v8_guassian2p0_cci.nc'
vfuture = '/users/global/cornkle/w2018_bamba/qrparm.cci.4km.nc'
topo = '/users/global/cornkle/w2018_bamba/topo.nc'
dummy = '/users/global/cornkle/w2018_bamba/xmhkga.pc20140405_00.nc'
out = '/users/global/cornkle/w2018_bamba/clean_files_raw/'


box = [150, 900, 93, 400]

vpast_ds = xr.open_dataset(vpast)
vcurrent_ds = xr.open_dataset(vfuture)
topo_ds = xr.open_dataset(topo, decode_times=False)
dummy_ds = xr.open_dataset(dummy)

dummy_ds = dummy_ds[['longitude_t', 'latitude_t']]
pdb.set_trace()
dummy_ds['vfraction_past'] = (['grid_latitude_t', 'grid_longitude_t'], vpast_ds['land_cover_fraction'][0,:,:].squeeze().values)
dummy_ds['vfraction_current'] = (['grid_latitude_t', 'grid_longitude_t'], vcurrent_ds['land_cover_fraction'][0,:,:].squeeze().values)
dummy_ds['topo'] = (['grid_latitude_t', 'grid_longitude_t'],topo_ds['ht'].squeeze().values)
dummy_ds['topo_stddev'] = (['grid_latitude_t', 'grid_longitude_t'], topo_ds['field150'].squeeze().values)


dummy_ds = dummy_ds.isel(grid_longitude_t=slice(box[0], box[1]), grid_latitude_t=slice(box[2], box[3]))

dummy_ds.to_netcdf(out+'ancils_vera.nc')



