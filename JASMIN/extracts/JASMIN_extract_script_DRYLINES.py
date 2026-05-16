import xarray as xr
import glob
import os
import itertools
from JASMIN import constants as cnst, MetUM_variables as mv
import numpy as np
import pdb

def cut_box():

    fpath = '/home/users/cornkle/CP4home/'
    tag = 'fut'
    full_box = [10,38,-35,-16]
    months = [10,11,12,1,2,3]
    vars_read = ['q2', 'lsRain', 'lw_out_PBLtop','p_srfc', 'u10', 'v10', 't2', 'SM', 'sh'] # ['rh_pl' , 'sh', 'colDryMass', 'colWetMass', 
    for vv in vars_read:
        
        inpath = fpath + 'CP4'+tag+'_EDW/' +  vv + '/'
        outpath = fpath + 'DRYLINE_SA_'+tag + '/' + vv + '/'
    
        
        if not os.path.exists(outpath):
           os.makedirs(outpath)
        for ff in glob.glob(inpath+'*.nc'):
            fmonth = ff[-11:-9]

            fout = ff.replace('CP4'+tag+'_EDW', 'DRYLINE_SA_'+tag)
            print('Doing', fout)
            if os.path.isfile(fout):
                print('File exists, continue')
                continue

            if int(fmonth) not in months:
               print('Wrong month')
               continue
            da = xr.open_dataarray(ff)
            da = da.sel(longitude=slice(full_box[0],full_box[1]), latitude=slice(full_box[2], full_box[3]), time=da.time.dt.month.isin(months)).load()
            
            da.name = vv
            encoding = {vv: {'complevel': 5, 'zlib': True}}
            try:
            	da.to_netcdf(fout , format='NETCDF4', encoding=encoding)
            except:
                pdb.set_trace()
            da.close()
            del da

            print('Wrote '+ fout)



