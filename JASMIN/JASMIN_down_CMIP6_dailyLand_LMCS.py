from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
import zarr
import gcsfs
import os
import glob
import pdb


storage = '/gws/nopw/j04/lmcs/CMIP6/daily_surface/'
df = pd.read_csv('https://storage.googleapis.com/cmip6/cmip6-zarr-consolidated-stores.csv')
df.head()

var = 'hfss'    # [prw, ua ]
sce = 'ssp585' # ssp585 , historical - run scenario first as fewer models did that!
tab = 'day'
variable = "'"+var+"'"
scenario = "'"+sce+"'"
table_id = "'"+tab+"'" 

df_prw = df.query("table_id == "+table_id+" & variable_id == "+variable+" & experiment_id == "+scenario +" & member_id == 'r1i1p1f1'")

def read_simulation(zstore):

    # this only needs to be created once
    gcs = gcsfs.GCSFileSystem(token='anon')

    # create a mutable-mapping-style interface to the store
    mapper = gcs.get_mapper(zstore)

    # open it using xarray and zarr
    ds = xr.open_zarr(mapper, consolidated=True)
    return ds


############## for all 2d JAS mean variables
save_path = storage+sce+'/'+var+'/'
model_list = []
model_names = []
file_name = []

for ids in range(len(df_prw)):  #len(df_prw)
    
    zstore = df_prw.zstore.values[ids]
    
    if 'r1i1p1' not in zstore:
        print('No r1i1p1 experiment, continue')
        continue
    #ipdb.set_trace()
    if sce == 'historical':
        if not zstore.split('/')[-8] in tllist:
            #ipdb.set_trace()
            print('Model does not exist for prw future, continue')
            continue
            
    if sce=='historical':
        y1 = 1980
        y2 = 2010
    else:
        y1 = 2070
        y2 = 2100
    
    ds = read_simulation(zstore)
    fn = ds.attrs['parent_mip_era']+'_'+ds.attrs['source_id']+'_'+ds.attrs['experiment_id']+'_'+ds.attrs['frequency']+'_'+ds.attrs['nominal_resolution']+'_'+ds.attrs['variant_label']+'_'+str(ds['time.year'].values[0])+'-'+str(ds['time.year'].values[-1])+'_'+ds.attrs['variable_id']+'_cutY'+str(y1)+'-'+str(y2)+'.nc'
    
    
    model_test = glob.glob(save_path+ds.attrs['parent_mip_era']+'_'+ds.attrs['source_id']+'*_r1i1p1*.nc')
    if model_test != []:
        print('Model already exists, continue')
    
    if os.path.isfile(save_path+fn):
        print(fn+' exists, continue.')
        continue
        
    if ds.attrs['source_id'] not in tnames.values:  # only download where ssp858 prw is available
        continue
        

#     if ('EC-Earth3-Veg' in ds.attrs['source_id']) | ('HadGEM3-GC31-MM' in ds.attrs['source_id']):     #CNRM-CM6-1-HR
#         continue
        
    print('Doing ', fn)
    
        
    ismask = (ds['time.year']>=y1) & (ds['time.year']<=y2) #& (ds['time.month']>=7) & (ds['time.month']<=9)

    ds = ds.isel(time=ismask)
    
    ds_jas_shift = shift_lons_data(ds[var])
    ds_jas_WA = ds_jas_shift#.sel(lat=slice(0,33), lon=slice(-19,25))# , lon=slice(-19,25))

    #ds_jas_WA = ds_jas_WA.where(ds_jas_WA > (0.1/3600)).quantile(0.95, dim='time', skipna=True) #groupby('time.year').mean('time')
    
    ds_jas_WA.attrs = ds.attrs
    
    file_name.append(fn)
    enc = {var: {'complevel': 5, 'shuffle': True, 'zlib': True}}
    ds_jas_WA.to_netcdf(save_path + fn, encoding=enc)
