import xarray as xr
import glob
import os
import itertools
from JASMIN import constants as cnst, MetUM_variables as mv
import numpy as np
import pdb

def run(orig_names=False):

    fpath = '/home/users/cornkle/runscript/in25'
    outpath = '/home/users/cornkle/runscript/out25_day'

    temp_box = [-18+360,35+360, 3.5, 30]
    months = [6,9] # March-May

    dic = {

        #'v_pl' : ([temp_box], ['daily'], [ 650, 850, 925], []),
        #'lsRain' : ([temp_box], ['daily'], [], []),
        #'totRain' : ([temp_box], ['daily'], [], []),
        #'colDryMass' : ([temp_box], ['daily'], [], []),
        #'colWetMass' : ([temp_box], ['daily'], [], []),
        #'u_pl' : ([temp_box], ['daily'], [ 650, 850, 925], []),
        'sh' : ([temp_box], ['daily'], [], []),
        'lh' : ([temp_box], ['daily'], [], []),
        'SM' : ([temp_box], ['daily'], [], []),
        'q2' : ([temp_box], ['daily'], [], []),
        't2' : ([temp_box], ['daily'], [], []),
        't_pl' : ([temp_box], ['daily'], [650,850,925], []),
        'q_pl' : ([temp_box], ['daily'], [650,850,925], []),
        'v_pl' : ([temp_box], ['daily'], [650,850,925], []),

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
        outfolder = outpath +os.sep +  k + '_daily' 
        files = glob.glob(infolder + os.sep + var+'*.nc' )
        for f in files:
            
            fname = os.path.basename(f)
            outname = fname.replace(var, k+'_fullPL_')
            outfile = outfolder + os.sep + outname
            fmonth = outfile[-11:-9]
            #pdb.set_trace()
            if (int(fmonth)<months[0]) | (int(fmonth)>months[1]):
                print('Wrong month, continue')
                continue

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
            try:
                mins = np.unique(ds[var]['time.minute'])[0]  # 0 in most cases, but may be 30 for fluxes
            except:
                if var == 'e08223':
                    var = 'SM'
                    mins = np.unique(ds[var]['time.minute'])[0]
             
            if dinfo[3] != []:
                ds = ds.isel(time=(([np.in1d(ds['time.hour'].values, dinfo[3])][0]) & (ds['time.minute']==mins))) 
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
                           da = cut['a04203'] # stupid lsRain_hFreq proble
                        except KeyError:
                           print('KEY ERROR, name missing')
                           pdb.set_trace()

                if pres != []:
                   # da = da.sel(pressure=pres)
                     da = da.sel(pressure=slice(pres[0],pres[-1]))
                if agg != 'keep':
                    da.values[da.values==0] = np.nan
                    if "Rain" not in k:

                    #da = da.mean('time')
                        da = da.resample(time='D').mean('time')#.isel(time=slice(0,-1))
                    else:
                        da = da.resample(time='D').sum('time')#.isel(time=slice(0,-1)) # mm/day
                    
                    if ("25" not in fpath) | ("_pl" in k):                      
                        if len(da.time)>1:
                           da = da.isel(time=slice(0,-1))
                #pdb.set_trace()
                comp = dict(zlib=True, complevel=5)

                da.name = k
                da = da.assign_coords(longitude=da.longitude.values-360)
                #encoding = {var: comp for var in da.data_vars}
                encoding = {k: {'complevel': 5, 'zlib': True}}
                if not os.path.exists(outfolder):
                    os.makedirs(outfolder)

                da.to_netcdf(outfolder + os.sep + outname , format='NETCDF4', encoding=encoding)
                da.close()
                del da

                print('Wrote '+ outfolder + os.sep + outname)



