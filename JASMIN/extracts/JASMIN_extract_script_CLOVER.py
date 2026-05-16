import xarray as xr
import glob
import os
import itertools
from JASMIN import constants as cnst, MetUM_variables as mv
import numpy as np
import pdb

def run(orig_names=False):

    fpath = '/home/users/cornkle/runscript/fut_in'
    outpath = '/home/users/cornkle/runscript/fut_out'

    local_box = [-18+360,14+360, 3.5, 14]
    temp_box = [-18+360,35+360, 3.5, 30]
    #hq_box = [-18+360,-10+360, 12, 17]
    hq_box = [-18+360,25+360, 3.5, 20]
    months = [8,8] # March-May 6,9

    dic = {

        #'t2' : ([temp_box], ['keep'], [], [12,3,6]),
        #'u10': ([temp_box], ['keep'], [], [12,3,6]),
        #'v10': ([temp_box], ['keep'], [], [12,3,6]),
        #'lw_out_PBLtop' : ([temp_box], ['keep'], [], []),
        #'v_pl' : ([temp_box], ['keep'], [200, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 925], [12,6]),
        #'u_pl' : ([temp_box], ['keep'], [200, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 925], [12,6]),
        #'t_pl' : ([temp_box], ['keep'], [200, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 925], [12,6]),
        #'geoH_pl' : ([temp_box], ['keep'], [200, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 925], [12,6]),
        'omega_pl' : ([temp_box], ['keep'], [300,650,850], [6,9,12,15,18]),
        #'lsRain' : ([temp_box], ['keep'], [], []),
        #'totRain' : ([temp_box], ['keep'], [], []),
        #'q_pl' : ([temp_box], ['keep'], [200, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 925], [12,6]),
        #'colDryMass' : ([temp_box], ['keep'], [], [12,3,6]),
        #'colWetMass' : ([temp_box], ['keep'], [], [12,3,6]),
        #'lsRain_hFreq' : ([hq_box], ['keep'], [], []),
        #'colDryMass' : ([temp_box], ['daily'], [], []),
        #'colWetMass' : ([temp_box], ['daily'], [], []),
        #'u_pl' : ([temp_box], ['daily'], [600, 650,925], []),
        #'rh_pl' : ([temp_box], ['keep'], [400, 500, 600, 700, 850, 925], [12]),
        #'sh' : ([temp_box], ['keep'], [], [12,15,18,21]),
        #'lh' : ([temp_box], ['keep'], [], [12,15,18,21]),
        #'SM' : ([temp_box], ['keep'], [], [12,15,18,21]),
        #'q2' : ([temp_box], ['keep'], [], [12,3,6]),
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
            if (var == 'e08223') & ('fut_in25' in infolder):
                ds = ds.rename({'sm' : 'e08223', 't':'time', 'level6' : 'depth'})

            mins = np.unique(ds[var]['time.minute'])[0]  # 0 in most cases, but may be 30 for fluxes
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
                    #da = da.mean('time')
                    da = da.resample(time='D').mean('time').isel(time=slice(0,-1))
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



