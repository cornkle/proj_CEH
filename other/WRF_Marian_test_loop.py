import glob
import ipdb
import os
import salem
import sys
import multiprocessing

varlist_2d = ['ALBEDO', 'LH',  'HFX', 'OLR', 'PBLH',
           'PSFC' , 'Q2','QFX','PRCP','SNOWC',
           'LWDNB', 'SWDNB','SWUPB', 'LWUPB', 'T2C', 'TH2', 'TSK', 'U10', 'V10']

static = ["LANDMASK", "LU_INDEX", "IVGTYP", "HGT"]
varlist_3d = ['U', 'V']

boxes = {'western' : [-18,10] , 'central' : [-10,0], 'eastern' : [0,10]
######################  INOUT CHANGES HERE
WRF_IN = '/home/cornkle/WRF_testfiles/'
WRF_OUT = '/home/cornkle/WRF_extractfiles/
VARLIST = '3d'
##################

def multiproc():

    pool = multiprocessing.Pool(processes=5)
    files = glob.glob(WRF_IN+'wrfout*')
    res = pool.map(rewrite_wrf, files)

    for f in files:
        ds = rewrite_wrf(f)



def rewrite_wrf(files):

    print('Doing', f)
    if varlist == '2d':
       varl = varlist_2d
    elif varlist == '3d':
       varl = varlist_3d
    elif varlist == 'static':
       varl = static
    elif varlist is None:
       print('No varlist set, please choose "2d", "3d" or "static"')
       return

    outp = f.replace('wrfout_', 'wrfout_small_'+varlist+'_')

    if os.path.isfile(outp):
	     print('File exists, continue')
	     return
    ds = salem.open_wrf_dataset(f, decode_times=False)
    ds = ds[varl]
    new_lat = ds['lat'].data[:, 0]
    new_lon = ds['lon'].data[0, :]
    south_north = ds['south_north']
    west_east = ds['west_east']

    ds = ds.rename({'south_north': 'latitude', 'west_east': 'longitude'})
    ds = ds.assign_coords({'latitude': new_lat.data, 'longitude': new_lon.data})
    ds = ds.assign_coords({'south_north': ('latitude', south_north.data), 'west_east': ('longitude', west_east.data)})
    
    ds = ds.drop('xtime')
    
    if box is not None:
       ds.sel(longitude=slice(box[0], box[1]), latitude=slice(box[2],box[3]))

    comp = dict(zlib=True, complevel=5)
    encoding = {var: comp for var in ds.data_vars}
    ds.to_netcdf(outp , mode='w', encoding=encoding, format='NETCDF4')
    print(outp, 'written')
    del ds


