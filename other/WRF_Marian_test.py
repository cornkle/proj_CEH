import glob
import ipdb
import os
import salem
import sys

varlist = ['ALBEDO', 'LH',  'HFX', 'OLR', 'PBLH',
           'PSFC' , 'Q2','QFX','PRCP','SNOWC',
           'LWDNB', 'SWDNB','SWUPB', 'LWUPB', 'T2C', 'TH2', 'TSK', 'U10', 'V10', 'RH2', 'col_int_QVAPOR', 'col_int_QICE']

hydro = ['SFROFF', 'UDROFF', 'TR', 'RUNSB', 'RUNSF']

#LANDMASK, SST, 'LU_INDEX',

static = ["LANDMASK", "LU_INDEX", "IVGTYP", "HGT"]

def rewrite_wrf(disk_path):

    #for f in glob.glob(disk_path):
    f = disk_path
    print('Doing', f)

    outp = f.replace('wrfout_', 'wrfout_small_')

    fbase = os.path.basename(f)
    #     if os.path.isfile(outp + fbase):
    #         print('File exists, continue')
    #         continue
    ds = salem.open_wrf_dataset(f, decode_times=False)
    ds = ds[varlist]
    new_lat = ds['lat'].data[:, 0]
    new_lon = ds['lon'].data[0, :]
    south_north = ds['south_north']
    west_east = ds['west_east']

    ds = ds.rename({'south_north': 'latitude', 'west_east': 'longitude'})
    ds = ds.assign_coords({'latitude': new_lat, 'longitude': new_lon})
    ds = ds.assign_coords({'south_north': ('latitude', south_north), 'west_east': ('longitude', west_east)})

    comp = dict(zlib=True, complevel=5)
    encoding = {var: comp for var in ds.data_vars}
    ds.to_netcdf(outp , mode='w', encoding=encoding, format='NETCDF4')
    print(outp, 'written')
    # ds_list.append(dsout)
    del ds



rewrite_wrf(sys.argv[1])
