import xarray as xr
import glob
import os
import itertools
from JASMIN import constants as cnst, MetUM_variables as mv
import numpy as np
import pdb

def run(orig_names=False):

    fpath = '/home/users/cornkle/runscript/in'
    outpath = '/home/users/cornkle/runscript/out_central'

    cloud_box = [6.5+360,18+360, -6.5, 5]
    temp_box = [0+360,30+360, -15, 15]
    months = [5,10] # March-May
    years = [1997,1998,1999,2000,2001,2002,2003,2004,2005,2006] # 2003
    pickhours = [0,3,6,9,12,15,18,21]
    picklevels = [600,700,750,800,850,900,925,950]
    dic = {
       #'lsRain' : ([temp_box], ['keep'], [], []),
       #'totRain' : ([temp_box], ['keep'], [], []),

       #'lowCloudFrac' : ([cloud_box], ['keep'], [], []),
       'lw_out_PBLtop' : ([cloud_box], ['keep'], [], []),
       # 't2' : ([cloud_box], ['keep'], [], []),
       # 'q2' : ([cloud_box], ['keep'], [], []),
       # 'slp' : ([temp_box], ['keep'], [], pickhours),
       # 'u10': ([cloud_box], ['keep'], [], []),
       # 'v10': ([cloud_box], ['keep'], [], []),
       # 'v_pl' : ([temp_box], ['keep'], picklevels, pickhours),
       # 'u_pl' : ([temp_box], ['keep'], picklevels, pickhours),
       # 't_pl' : ([temp_box], ['keep'], picklevels, pickhours),
       #'geoH_pl' : ([temp_box], ['keep'], picklevels, pickhours),
       # 'q_pl' : ([temp_box], ['keep'], picklevels, pickhours),
       # 'rh_pl' : ([temp_box], ['keep'], picklevels, pickhours)
    }
    keys = dic.keys()

    for k in keys:

        info = cnst.VARDIC[k]
        dinfo = dic[k]
        var = mv.create_CP4_filename(k)

        if not orig_names:
            pathvar = k
        else:
            pathvar = var
        
        infolder = fpath+os.sep + pathvar
        outfolder = outpath +os.sep +  k 
        files = glob.glob(infolder + os.sep + var+'*.nc' )
        #pdb.set_trace()
        if not os.path.exists(outfolder):
           os.makedirs(outfolder)
        for f in files:
            
            fname = os.path.basename(f)
            outname = fname.replace(var, k+'_CentralA_')
            if ("_pl" in var) & ("1997" in fname):
               print('Skipping 1997 for _pl variable')
               continue
            outfile = outfolder + os.sep + outname
            print('Trying '+outfile)
            if os.path.isfile(outfile):
                print('File already exists, continue.')
                continue
            #if '199809220100-199809230000' in outfile:
            # continue
            if years is not None:
                   print('Isyear', outfile[-15:-11])
                   #pdb.set_trace()
                   if int(outfile[-15:-11]) not in years:
                       print('Wrong year, continue')
                       continue
            
            try:
             ds = xr.open_dataset(f)
            except OSError:
             print('Netcdf OSError')
             pdb.set_trace()
             continue

            if (ds['time.month'][0]<months[0]) | (ds['time.month'][0]>months[1]):
               print('Wrong month, continue') 
               continue
            

            if dinfo[3] != []:
                ds = ds.isel(time=(([np.in1d(ds['time.hour'].values, dinfo[3])][0]) & (ds['time.minute']==0))) 
            box = dinfo[0]

            for id, b in enumerate(box):

                agg = dinfo[1][id]
                pres = dinfo[2]
                cut = ds.sel(longitude=slice(b[0], b[1]), latitude=slice(b[2], b[3]))
                  
                try:
                  da = cut[var]
                except KeyError:
                   try:
                      da = cut['c03238'] # stupid t2 problem
                   except KeyError:
                      try:
                         da = cut['a04203'] # stupid lsRain_hFreq problem
                      except:
                         print('STUPID KEY ERROR') 
                         pdb.set_trace()

                if pres != []:
                   # da = da.sel(pressure=pres)
                     da = da.sel(pressure=slice(pres[0],pres[-1]))

                if agg != 'keep':
                    da.values[da.values==0] = np.nan
                    #da = da.mean('time')
                    da = da.resample(time='D').mean('time').isel(time=slice(0,-1))
                    #pdb.set_trace()

                comp = dict(zlib=True, complevel=5)

                da.name = k
                #pdb.set_trace()
                #da.longitude.values = da.longitude.values-360
                da = da.assign_coords(longitude=(da.longitude-360))
                #encoding = {var: comp for var in da.data_vars}
                encoding = {k: {'complevel': 5, 'zlib': True}}
                if not os.path.exists(outfolder):
                    os.makedirs(outfolder)
                try:  
                   da.to_netcdf(outfolder + os.sep + outname , format='NETCDF4', encoding=encoding)
                except:
                   pdb.set_trace()
                da.close()

                print('Wrote '+ outfolder + os.sep + outname)



