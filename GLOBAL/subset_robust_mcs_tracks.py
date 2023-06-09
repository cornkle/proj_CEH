import sys
import os
import numpy as np
import xarray as xr
import glob
from GLOBAL import glob_util
from utils import constants as cnst


def subset_file(infile, outdir, lon_box, lat_box):

    status = 0

    # Make output file name

    # Strip the path to get input file name
    infilename = os.path.basename(infile)
    outfile = f'{outdir}/{infilename}'

    if os.path.isfile(outfile):
        print('File exists')
        return status

    print('Starting', outfile)

    # Read input MCS track file
    ds = xr.open_dataset(infile)

    # Get track initiation lat/lon location
    lon0 = ds.meanlon.isel(times=0)
    lat0 = ds.meanlat.isel(times=0)

    # Get track index for MCS initiation location within a region
    trackid = np.where((lon0 >= np.min(lon_box)) & (lon0 <= np.max(lon_box)) & \
                        (lat0 >= np.min(lat_box)) & (lat0 <= np.max(lat_box)))[0]


    # Subset the tracks from the dataset
    dsout = ds.sel(tracks=trackid, drop=True)
    # Add subset lat/lon box to global attribute
    dsout.attrs['lon_box'] = lon_box
    dsout.attrs['lat_box'] = lat_box



    # Write output file
    comp = dict(zlib=True, complevel=5)
    encoding = {var: comp for var in dsout.data_vars}
    dsout.to_netcdf(path=outfile, mode='w', encoding=encoding, format='NETCDF4', unlimited_dims='tracks')
    print(f'Output saved: {outfile}')

    status = 1
    return status


def main(reg):
    """

    :param reg: sub-region string as defined in glob_util
    """
    # Read from input

    consts = glob_util.MREGIONS[reg]
    box = consts[0]

    inpath =  '/media/ck/LStorage/MCS_Feng/global_v2/tracks/global/' #sys.argv[5]
    outdir = '/media/ck/LStorage/MCS_Feng/global_v2/tracks/custom/'+reg+'/' #sys.argv[6]

    # Subset region boundary
    lon_box = [box[0], box[1]]
    lat_box = [box[2], box[3]]
    os.makedirs(outdir, exist_ok=True)

    # Call function
    for infile in glob.glob(inpath+'*.nc'):
        status = subset_file(infile, outdir, lon_box, lat_box)



def filter(path):
    """
    Used to remove variables from existing track files.
    :param path: path to track files
    :return: slimmed-down track files
    """
    for f in glob.glob(path+'/*1.nc'):
        ds = xr.open_dataset(f)
        drops = ['datetimestring', 'movement_r', 'movement_theta', 'movement_r_meters_per_second',
                              'movement_time_lag', 'movement_storm_x', 'movement_storm_y', 'pf_nuniqpix', 'location_idx', 'pixel_duration',
                             'pixel_pcp', 'pf_skewness', 'direction', 'eccentricity', 'mcs_status', 'uspeed', 'vspeed']
        ds = ds.drop(drops)

        comp = dict(zlib=True, complevel=5)
        encoding = {var: comp for var in ds.data_vars}

        fout = f.replace('.nc', '_cut.nc')
        ds.to_netcdf(path=fout, mode='w', encoding=encoding, format='NETCDF4', unlimited_dims='tracks')



for reg in glob_util.MREGIONS.keys():
    main(reg)