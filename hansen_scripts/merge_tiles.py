from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
import os
from datetime import datetime
from tqdm import tqdm


def mosaic_cover2000_tiles(tile_directory, mask_water=False):

    global_lon = np.arange(-180, 180, 0.05) + 0.5*0.05
    global_lat = np.arange(-90, 90, 0.05)[::-1] + 0.5*0.05

    tile_files = os.listdir(tile_directory)
    cover_tile_files = [f for f in tile_files if 'Hansen_GFC_percent_cover_2000' in f]

    f = Dataset('../data/hansen_0pt05/global_percent_cover_2000_0pt05deg.nc', 'w', format='NETCDF4')
    f.description = f'Percentage tree canopy cover in 2000 on 0.05deg CMG.'
    today = datetime.today()
    f.history = "Created " + today.strftime("%d/%m/%y")

    f.createDimension('lat', global_lat.size)
    f.createDimension('lon', global_lon.size)
    latitude = f.createVariable('lat', 'f4', 'lat')
    latitude[:] = global_lat
    latitude.standard_name = 'lat'
    latitude.long_name = 'latitude'
    latitude.units = 'degrees_north'
    longitude = f.createVariable('lon', 'f4', 'lon')
    longitude[:] = global_lon
    longitude.standard_name = 'lon'
    longitude.long_name = 'longitude'
    longitude.units = 'degrees_east'
    cover_2000 = f.createVariable(f'tree_cover_2000', 'f4', ('lat', 'lon'), fill_value=-999)
    cover_2000.standard_name = 'tree_cover_2000'
    cover_2000.units = 'Percent'
    cover_2000[:] = np.zeros((global_lat.size, global_lon.size), dtype=float)
    if mask_water:
        cover_2000 -= 999

    for tile in tqdm(cover_tile_files, desc='cover 2000'):
        with Dataset(f'{tile_directory}/{tile}', 'r') as tile_data:
            tile_lats = tile_data.variables['lat'][:]
            tile_lons = tile_data.variables['lon'][:]
            tile_loss = tile_data.variables['tree_cover_2000'][:]
            if not mask_water:
                tile_loss = ma.filled(tile_loss, 0.)
            lat_idx = np.argmin(np.abs(global_lat - tile_lats.min()))
            lon_idx = np.argmin(np.abs(global_lon - tile_lons.min()))
            cover_2000[lat_idx-199:lat_idx+1, lon_idx:lon_idx+200] = tile_loss

    cover_2000[:] = np.flipud(cover_2000[:])
    latitude[:] = latitude[::-1]

    f.close()


def mosaic_loss_tiles(tile_directory, mask_water=False):

    global_lon = np.arange(-180, 180, 0.05) + 0.5*0.05
    global_lat = np.arange(-90, 90, 0.05)[::-1] + 0.5*0.05

    tile_files = os.listdir(tile_directory)
    loss_tile_files = [f for f in tile_files if 'Hansen_GFC_percent_net_loss' in f]

    f = Dataset('../data/hansen_0pt05/global_percent_forest_loss_0pt05deg.nc', 'w', format='NETCDF4')
    f.description = f'Percentage net loss of global forest coverage from 2000-2019 on 0.05deg CMG.'
    today = datetime.today()
    f.history = "Created " + today.strftime("%d/%m/%y")

    f.createDimension('lat', global_lat.size)
    f.createDimension('lon', global_lon.size)
    latitude = f.createVariable('lat', 'f4', 'lat')
    latitude[:] = global_lat
    latitude.standard_name = 'lat'
    latitude.long_name = 'latitude'
    latitude.units = 'degrees_north'
    longitude = f.createVariable('lon', 'f4', 'lon')
    longitude[:] = global_lon
    longitude.standard_name = 'lon'
    longitude.long_name = 'longitude'
    longitude.units = 'degrees_east'
    forest_change = f.createVariable(f'forest_cover_loss', 'f4', ('lat', 'lon'), fill_value=-999)
    forest_change.standard_name = 'forest_cover_loss'
    forest_change.units = 'Percent'
    forest_change[:] = np.zeros((global_lat.size, global_lon.size), dtype=float)
    if mask_water:
        forest_change -= 999

    for tile in tqdm(loss_tile_files, desc='percent loss'):
        with Dataset(f'{tile_directory}/{tile}', 'r') as tile_data:
            tile_lats = tile_data.variables['lat'][:]
            tile_lons = tile_data.variables['lon'][:]
            tile_loss = tile_data.variables['forest_cover_loss'][:]
            if not mask_water:
                tile_loss = ma.filled(tile_loss, 0.)
            lat_idx = np.argmin(np.abs(global_lat - tile_lats.min()))
            lon_idx = np.argmin(np.abs(global_lon - tile_lons.min()))
            forest_change[lat_idx-199:lat_idx+1, lon_idx:lon_idx+200] = tile_loss

    forest_change[:] = np.flipud(forest_change[:])
    latitude[:] = latitude[::-1]

    f.close()


if __name__ == '__main__':
    mosaic_cover2000_tiles('../data/hansen_0pt05')
    mosaic_loss_tiles('../data/hansen_0pt05')