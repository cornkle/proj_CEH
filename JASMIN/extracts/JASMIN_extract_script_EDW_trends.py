import xarray as xr
import glob
import os
import itertools
from JASMIN import constants as cnst, MetUM_variables as mv
import numpy as np
import pdb

def run(orig_names=False):

    fpath = '/home/users/cornkle/runscript/in'
    outpath = '/home/users/cornkle/runscript/out_EDW'

    local_box = [-18+360,14+360, 3.5, 14]
    temp_box = [-18+360,35+360, 3.5, 30]
    #hq_box = [-18+360,-10+360, 12, 17]
    hq_box = [-18+360,25+360, 3.5, 20]
    months = [1,12] # March-May
    full_box = [-18+360, 52+360, -36,39]

    dic = {

        't2' : ([full_box], ['keep'], [], []),
        'q2' : ([full_box], ['keep'], [], []),
        #'t2_daily' : ([full_box], ['keep'], [], []),
        'u10': ([temp_box], ['keep'], [], [12,3,6]),
        'v10': ([temp_box], ['keep'], [], [12,3,6]),
        'lw_out_PBLtop' : ([full_box], ['keep'], [], []),
        #'v_pl' : ([temp_box], ['keep'], [200, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 925], [12,6]),
        #'u_pl' : ([temp_box], ['keep'], [200, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 925], [12,6]),
        #'t_pl' : ([temp_box], ['keep'], [200, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 925], [12,6]),
        #'geoH_pl' : ([temp_box], ['keep'], [200, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 925], [12,6]),
        #'omega_pl' : ([temp_box], ['keep'], [650,300], [12,3,6]),
        'lsRain' : ([full_box], ['keep'], [], []),
        #'totRain' : ([temp_box], ['keep'], [], []),
        #'q_pl' : ([temp_box], ['keep'], [200, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 925], [12,6]),
        #'colDryMass' : ([temp_box], ['keep'], [], [12,3,6]),
        #'colWetMass' : ([temp_box], ['keep'], [], [12,3,6]),
        #'lsRain_hFreq' : ([hq_box], ['keep'], [], []),
        #'colDryMass' : ([temp_box], ['daily'], [], []),
        #'colWetMass' : ([temp_box], ['daily'], [], []),
        #'u_pl' : ([temp_box], ['daily'], [600, 650,925], []),
        #'rh_pl' : ([temp_box], ['keep'], [400, 500, 600, 700, 850, 925], [12]),
        #'sh' : ([full_box], ['keep'], [], [10,12,15,18]),
        #'lh' : ([full_box], ['keep'], [], [10,12,15,18]),
        #'SM' : ([full_box], ['keep'], [], [10,12,15,18])
        'p_srfc' : ([full_box], ['keep'], [], []),
        'u10' : ([full_box], ['keep'], [], []),
        'v10' : ([full_box], ['keep'], [], [])
    }
    keys = dic.keys()

    for k in keys:
        print('Doing ', k)
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
        for f in files:
            
            fname = os.path.basename(f)
            outname = fname.replace(var, k+'_fullPL_')
            outfile = outfolder + os.sep + outname
            print('Trying '+outfile)
            if os.path.isfile(outfile):
                print('File already exists, continue.')
                continue
            #if '199809220100-199809230000' in outfile:
            # continue
            try:
             ds = xr.open_dataset(f).load()
            except OSError:
             print('Netcdf OSError')
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
                da = da.assign_coords(longitude=da.longitude.values-360)
                #encoding = {var: comp for var in da.data_vars}
                encoding = {k: {'complevel': 5, 'zlib': True}}
                if not os.path.exists(outfolder):
                    os.makedirs(outfolder)

                da.to_netcdf(outfolder + os.sep + outname , format='NETCDF4', encoding=encoding)
                da.close()
                del da

                print('Wrote '+ outfolder + os.sep + outname)



