"""
Percentage forest cover loss on 0.05 degree grid.

This script works with the Hansen et al. (2013) Global Forest Change dataset
v1.7 (2000-2019). The data must be downloaded before running this script.
This can be done using download_hansen_tiles.py: required fields
are lossyear, gain, cover2000 and datamask. All downloaded files should
be placed in a single directory.

Once the files are downloaded into directory /my/example/directory/,
this script can be run as
>> python aggregate_forest_loss.py /my/example/directory/

This will use the Hansen et al. 30m-resolution data and compute the
loss in percent forest cover from 2000-2019 at 0.05 degree resolution.
For example, if a 0.05 degree pixel has 80% tree cover in 2000
and 60% tree cover in 2019, its computed value will be 20.
Negative values represent net gain. If a 0.05 degree pixel contains
no land in the 30m data, it will be masked.

It is assumed that a gain in forest at the 30m scale results
in 100% forest cover (for that 30m pixel). Therefore the net loss value
may be treated as a lower bound.

For each downloaded 10x10 degree tile, a corresponding NetCDF file
will be created in the directory ../data/hansen0pt05 containing
the percent loss data.

Bethan L. Harris, UK Centre for Ecology & Hydrology, 15th April 2021.
"""

import os, rasterio, sys
import numpy as np
from tqdm import tqdm
from netCDF4 import Dataset
from datetime import datetime
from multiprocessing import Pool


def get_data(filename):
    """
    Read data and coordinates from .tif file. Assumes required data is in band 1.

    Parameters:
    filename (str): path to .tif file to read.
    Returns:
    lon_west (float): longitude of western edge of data (degrees east).
    lat_south (float): latitude of southern edge of data (degrees north).
    """

    with rasterio.open(filename) as dataset:
        data = dataset.read(1)
        width = dataset.width
        height = dataset.height
        lon_west, _ = dataset.transform * (0, 0)
        _, lat_south = dataset.transform * (dataset.width, dataset.height)
    return lon_west, lat_south, data


def pixel_areas(lons, lats):
    """
    Compute areas of grid boxes comprising a latitude-longitude grid using spherical geometry.

    Assumes resolution of 1/4000 deg (as in Hansen et al. tiles).
    Parameters:
    lons (numpy 1D array of floats): longitudes of grid points (degrees east).
    lats (numpy 1D array of floats): latitudes of grid points (degrees north).
    Returns:
    grid_box_area (numpy 2D array of floats): area of each lon/lat grid box (m^2).
    """

    deg_spacing = 0.00025
    radius_earth = 6371000.
    lat_grid, _ = np.meshgrid(lats, lons, indexing='ij')
    grid_box_area = (np.pi / 180.) * (radius_earth ** 2) * deg_spacing * np.abs(
    np.sin(np.deg2rad(lat_grid+0.5*deg_spacing)) - np.sin(np.deg2rad(lat_grid-0.5*deg_spacing)))
    return grid_box_area


def aggregate_tile(loss_filename):
    """
    For a 10x10 degree 30m-resolution tile from Hansen et al. dataset,
    read the forest loss, gain and 2000 tree cover data
    and aggregate this to percentage forest cover loss from 2000-2019 on a 0.05 degree grid.

    Percentage cover loss is defined as 
    (percent forest cover of 0.05 deg box in 2019) - (percent forest cover of 0.05 deg box in 2000).
    Negative values of loss are possible and correspond to net forest cover gain.
    Cover loss is filled with -999 if there is no land in the grid box.
    Parameters:
    loss_filename (str): path to the file containing the Hansen et al. forest loss year data for tile.
    Returns:
    lon_west (float): longitude of western edge of tile (degrees east).
    lat_south (float): latitude of southern edge of tile (degrees north).
    cover2000_array (numpy array of floats, 200x200): percent forest cover in 2000 at 0.05 deg.
    loss_array (numpy array of floats, 200x200): percent forest cover loss at 0.05 deg.
    """

    # Read data for tile: forest loss year, forest gain, tree cover in 2000, water mask.
    lon_west, lat_south, loss_data = get_data(loss_filename)
    gain_filename = loss_filename.replace('lossyear', 'gain')
    cover2000_filename = loss_filename.replace('lossyear', 'treecover2000')
    datamask_filename = loss_filename.replace('lossyear', 'datamask')
    _, _, gain_data = get_data(gain_filename)
    _, _, cover2000_data = get_data(cover2000_filename)
    _, _, datamask = get_data(datamask_filename)
    # Create arrays of 30m pixel lats/lons for computing pixel areas
    deg_spacing = 0.00025
    lats = np.arange(lat_south, lat_south + 10, deg_spacing)[::-1] + 0.5 * deg_spacing 
    lons = np.arange(lon_west, lon_west + 10, deg_spacing) + 0.5 * deg_spacing
    # Find which 30m pixels experience net loss/net gain of forest from 2000-2013 
    cover2000_array = np.zeros((200, 200), dtype=float)
    loss_array = np.zeros((200, 200), dtype=float)
    loss_pre_2012 = np.logical_and(loss_data > 0,  loss_data < 13)
    loss_in_2013 = (loss_data == 13)
    # Net loss if forest lost in 2013, or lost in/before 2012 and no gain from 2000-2012
    overall_loss_bool = np.logical_or(loss_in_2013, np.logical_and(loss_pre_2012, gain_data == 0))
    # Net gain if forest gained from 2000-2012 and no loss after 2012
    overall_gain_bool = np.logical_and(gain_data == 1, loss_data < 13)
    # Iterate through each 0.05 deg grid box covering the tile to aggregate loss/gain to net loss percentage
    block_size = 200
    number_blocks = 200
    for i in range(number_blocks):
        for j in range(number_blocks):
            loss_data_pixel = overall_loss_bool[i*block_size:(i+1)*block_size, j*block_size:(j+1)*block_size]
            gain_data_pixel = overall_gain_bool[i*block_size:(i+1)*block_size, j*block_size:(j+1)*block_size]
            cover2000_pixel = cover2000_data[i*block_size:(i+1)*block_size, j*block_size:(j+1)*block_size]/100.
            datamask_pixel = datamask[i*block_size:(i+1)*block_size, j*block_size:(j+1)*block_size]
            pixel_lats = lats[i*block_size:(i+1)*block_size]
            pixel_lons = lons[j*block_size:(j+1)*block_size]
            area = pixel_areas(pixel_lons, pixel_lats)
            total_land_area = area[datamask_pixel==1].sum()
            total_pixel_area = area.sum()
            if total_land_area > 0.:
                # Create percentage tree cover map for 2013
                final_cover = np.copy(cover2000_pixel)
                final_cover[loss_data_pixel] = 0. # Assume net forest loss results in 0% cover
                final_cover[gain_data_pixel] = 1. # Assume net forest gain results in 100% cover
                # Compute change in percentage cover 2000-2013.
                original_forest_cover = (cover2000_pixel * area).sum()/total_pixel_area
                final_forest_cover = (final_cover * area).sum()/total_pixel_area
                loss_forest_cover = original_forest_cover - final_forest_cover
                cover2000_array[i, j] = np.round(100*original_forest_cover, 2)
                loss_array[i, j] = np.round(100*loss_forest_cover, 2)
            else: # 0.05 deg box is entirely over water
                cover2000_array[i, j] = -999
                loss_array[i, j] = -999
    return lon_west, lat_south, cover2000_array, loss_array


def loss_tiff_list(data_directory):
    """
    Get a list of paths to all the files in the download directory that contain data for forest loss year.

    Parameters:
    data_directory (str): directory to which Hansen et al. tiles have been downloaded.
    Returns:
    all_loss_files (list of str): list of file paths for forest loss tiles.
    """

    all_files = os.listdir(data_directory)
    all_loss_files = ['{}/{}'.format(data_directory, f) for f in all_files if f.endswith('.tif') and 'lossyear' in f]
    return all_loss_files


def write_aggregated_tile(loss_tile_filename):
    """
    For a 10x10 degree 30m-resolution tile from Hansen et al. dataset,
    compute percentage forest cover loss on a 0.05 deg grid and save as NetCDF file.

    NetCDF files are saved in directory '../data/hansen0pt05' with filename format
    "Hansen_GFC_percent_net_loss_{LON}_{LAT}.nc", where LON is the longitude and LAT is the latitude of
    the north-west corner of the tile (as in the original Hansen et al. download tiles).
    Parameters:
    loss_tile_filename (str): path to the file containing the Hansen et al. forest loss year data for tile.
    Returns:
    None
    """

    # Check if save directory exists (create if not) and if tile has already been aggregated
    save_directory = '../data/hansen_0pt05/2013'
    if not os.path.isdir(save_directory):
        os.makedirs(save_directory)
    save_filename_cover2000 = f'{save_directory}/Hansen_GFC_percent_cover_2000_{loss_tile_filename[-12:-4]}.nc'
    save_filename_loss = f'{save_directory}/Hansen_GFC_percent_net_loss_2013_{loss_tile_filename[-12:-4]}.nc'
    if os.path.exists(save_filename_loss):
        print(f'skipping {save_filename_loss}: already exists')
    else: 
        # Compute percentage cover loss aggregated to 0.05 deg
        lon_west, lat_south, cover2000_data, loss_data = aggregate_tile(loss_tile_filename)
        # Save tiles of 0.05 deg data
        lon_east = lon_west + 10
        lat_north = lat_south + 10
        today = datetime.today()
        # Create 0.05 deg lat/lon grids for tile
        tile_lon = np.arange(lon_west, lon_east, 0.05) + 0.5*0.05
        tile_lat = np.arange(lat_south, lat_north, 0.05)[::-1] + 0.5*0.05
        # Save aggregated data for 2000 tree  cover
        # cover2000_file = Dataset(save_filename_cover2000, 'w', format='NETCDF4')
        # cover2000_file.history = "Created " + today.strftime("%d/%m/%y")
        # cover2000_file.description = f'Percentage tree canopy cover in 2000 on 0.05 deg CMG.'
        # cover2000_file.createDimension('lat', tile_lat.size)
        # cover2000_file.createDimension('lon', tile_lon.size)
        # latitude = cover2000_file.createVariable('lat', 'f4', 'lat')
        # latitude[:] = tile_lat
        # latitude.standard_name = 'lat'
        # latitude.long_name = 'latitude'
        # latitude.units = 'degrees_north'
        # longitude = cover2000_file.createVariable('lon', 'f4', 'lon')
        # longitude[:] = tile_lon
        # longitude.standard_name = 'lon'
        # longitude.long_name = 'longitude'
        # longitude.units = 'degrees_east'
        # forest_change = cover2000_file.createVariable(f'tree_cover_2000', 'f4', ('lat', 'lon'), fill_value=-999)
        # forest_change.standard_name = 'tree_cover_2000'
        # forest_change.units = 'Percent'
        # forest_change[:] = cover2000_data
        # cover2000_file.close()
        # Save aggregated data for percent forest loss 2000-2013
        loss_file = Dataset(save_filename_loss, 'w', format='NETCDF4')
        loss_file.history = "Created " + today.strftime("%d/%m/%y")
        loss_file.description = f'Percentage loss of global forest coverage from 2000-2013 on 0.05 deg CMG.'
        loss_file.createDimension('lat', tile_lat.size)
        loss_file.createDimension('lon', tile_lon.size)
        latitude = loss_file.createVariable('lat', 'f4', 'lat')
        latitude[:] = tile_lat
        latitude.standard_name = 'lat'
        latitude.long_name = 'latitude'
        latitude.units = 'degrees_north'
        longitude = loss_file.createVariable('lon', 'f4', 'lon')
        longitude[:] = tile_lon
        longitude.standard_name = 'lon'
        longitude.long_name = 'longitude'
        longitude.units = 'degrees_east'
        forest_change = loss_file.createVariable(f'forest_cover_loss', 'f4', ('lat', 'lon'), fill_value=-999)
        forest_change.standard_name = 'forest_cover_loss'
        forest_change.units = 'Percent'
        forest_change[:] = loss_data
        loss_file.close()


def write_all_tiles(data_directory):
    """
    For every 10x10 degree 30m-resolution tile from Hansen et al. dataset,
    compute percentage forest cover loss on a 0.05 deg grid and save as NetCDF file.

    NetCDF files are saved in directory '../data/hansen0pt05' with filename format
    "Hansen_GFC_percent_net_loss_{LON}_{LAT}.nc", where LON is the longitude and LAT is the latitude of
    the north-west corner of the tile (as in the original Hansen et al. download tiles).
    Parameters:
    data_directory (str): directory to which Hansen et al. tiles have been downloaded.
    Returns:
    None
    """

    all_loss_tiles = loss_tiff_list(data_directory)
    # Run 2 tiles at once. May need to be adjusted dependent on machine memory
    # If it's not working for unclear reasons, try processes=1
    with Pool(processes=2) as pool:
        pool.map(write_aggregated_tile, all_loss_tiles)


if __name__ == '__main__':
    # Read directory containing the downloaded Hansen et al. data from command line arguments
    # and compute 0.05 deg percent forest cover loss for all tiles
    hansen_download_directory = sys.argv[1]
    write_all_tiles(hansen_download_directory)
