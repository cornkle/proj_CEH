"""
This script subsets robust MCS tracks initiated within a region and saves to a new netCDF file.
"""
__author__ = "Zhe.Feng@pnnl.gov"
__date__ = "08-Dec-2020"

import sys
import os
import numpy as np
import xarray as xr
import glob

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


def main(reg, inds,lon1, lon2, lat1, lat2):
    # Read from input
    lon_min = lon1 #float(sys.argv[1])
    lon_max = lon2 #float(sys.argv[2])
    lat_min = lat1 #float(sys.argv[3])
    lat_max = lat2 #float(sys.argv[4])
    inpath =  '/media/ck/Elements/global/MCS_Feng/tracks/'+inds[1]+'/'#sys.argv[5]
    outdir = '/media/ck/Elements/global/MCS_Feng/tracks/custom/'+reg+'/' #sys.argv[6]

    # Subset region boundary
    lon_box = [lon_min, lon_max]
    lat_box = [lat_min, lat_max]
    # lon_box = [-15, 45]
    # lat_box = [-20, 30]
    # lon_box = [-15, 40]
    # lat_box = [-15, 25]

    # Define output directory
    # outdir = f'/global/cscratch1/sd/feng045/waccem/mcs_region/spac/subset_africa_stats/'
    os.makedirs(outdir, exist_ok=True)

    # Call function
    for infile in glob.glob(inpath+'*.nc'):
        status = subset_file(infile, outdir, lon_box, lat_box)

mregions = {'WAf' : [[-18,25,4,25], 'spac', 0], # last is hourly offset to UCT # 12
 'SAf' : [[20,35, -35,-15], 'spac', 2], # 10
 'india' : [[70,90, 5,30], 'asia', 5], # 7
 'china' : [[105,115,25,40], 'asia', 8 ], # 4
 'australia' : [[120,140,-23, -11], 'asia', 9], # 3
 'sub_SA' : [[-68,-47, -40, -20.5], 'spac', -4] , # 16
 'trop_SA' : [[-75, -50, -20, -5], 'spac', -5], # 17
 'GPlains' : [[-100,-90,32,47], 'nam', -6] # # 18

}

for reg in ['WAf', 'SAf', 'india', 'china', 'australia', 'sub_SA', 'trop_SA', 'GPlains']:
    inds = mregions[reg]
    lon1, lon2, lat1, lat2 = inds[0]
    main(reg, inds, lon1, lon2, lat1, lat2)