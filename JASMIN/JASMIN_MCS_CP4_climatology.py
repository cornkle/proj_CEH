import numpy as np
import xarray as xr
import os
import ipdb
import glob
import itertools
import datetime
import sys
import pandas as pd
from JASMIN import MetUM_variables as mu

def load_file(ffile,var):
    try:
        ds=xr.open_dataset(ffile)[var]
    except:
        #Deal with funny landseamask file
        ds=xr.open_dataset(ffile,decode_times=False)[var]
        ds=ds.rename({"rlat":"latitude","rlon":"longitude"})
    ds=ds.assign_coords({"longitude":(ds.longitude-360)})

    try:
        ds=ds.isel(pressure=slice(None,None,-1))
    except:
        pass
    return ds


def filtering(dar, v, outv):


    if (v == 'SM'):
        dar = dar.sel(depth=0.05)
        dar = dar.where(dar < 500, other=np.nan)

    if 'pressure' in dar.coords:

        try:
            dar.values[dar.values == 0] = np.nan  # potential missing value maskout
        except ValueError:
            ipdb.set_trace()
            print('Pressure value error!')
            return

    return dar



def file_save(cp_dir, out_dir, vars, datestring, box):

    #loop through every var
    for idx, outv in enumerate(vars):

        print('Variable ', outv)

        if ('_pl' in outv) | (outv == 'pblH'):
            if HOUR not in [0,3,6,9,12,15,18,21]:
                print('3-hourly data only for pressure level, continue with next var')

        orig_v = mu.create_CP4_filename(outv)

        savefile = out_dir + os.sep +  outv+'_'+datestring[4:6] + '-' + datestring[6:8] + '_'+ str(HOUR) + '_'+FTAG+'.nc'

        if len(glob.glob(savefile)) > 0:
            print(savefile, ' already exists, continue!')
            continue
        
        try:
            filepath = glob.glob(cp_dir+os.sep+str(outv)+os.sep+'*_'+datestring+'*.nc')[0]
        except IndexError:
            #ipdb.set_trace()
            print('No file found for var, continue')
            if (outv!='pblH'):
                ipdb.set_trace()
            continue

        #ddate = dar.time.values.item()
        dt_main = datetime.datetime(int(datestring[0:4]), int(datestring[4:6]), int(datestring[6:8]),
                               HOUR)
        mean_arr = []
        for yy in range(1998,2007):
                print('Doing year', yy)
                per_year = []
                dt = datetime.datetime(yy, int(datestring[4:6]), int(datestring[6:8]), HOUR)
                             
                for dd in [7,6,5,4,3,2,1,0]:
                    ndate = dt - pd.Timedelta(str(dd)+'days')
                    if ndate.day >30:
                        dday = 30
                    else:
                        dday = ndate.day
                    ndatestring = str(ndate.year)+str(ndate.month).zfill(2)+str(dday).zfill(2)
                    try:
                        nfile = glob.glob(cp_dir+os.sep+str(outv)+os.sep+'*_'+ndatestring+'*.nc')[0]
                    except:
                        continue
                    try:
                        meanarr = load_file(nfile, orig_v)
                        mmeanarr = meanarr.sel(longitude=slice(box[0],box[1]), latitude=slice(box[2],box[3]))
                    except OSError:
                        print('Couldnt find clim file, continue')
                        continue
                    mmeanarr = mmeanarr[mmeanarr['time.hour'] == HOUR].squeeze()
    
                    mmeanarr = filtering(mmeanarr, outv, outv)
                    per_year.append(mmeanarr)
    
                for dd in [7,6,5,4,3,2,1]:
                    ndate = dt + pd.Timedelta(str(dd)+'days')
                    if ndate.day >30:
                        dday = 30
                    else:
                        dday = ndate.day
                    ndatestring = str(ndate.year)+str(ndate.month).zfill(2)+str(dday).zfill(2)
                    try:
                        nfile = glob.glob(cp_dir+os.sep+str(outv)+os.sep+'*_'+ndatestring+'*.nc')[0]
                    except:
                        continue
                    try:
                        meanarr = load_file(nfile, orig_v)
                        mmeanarr = meanarr.sel(longitude=slice(box[0],box[1]), latitude=slice(box[2],box[3]))
                    except OSError:
                        print('Couldnt find clim file, continue')
                        continue
                    mmeanarr = mmeanarr[mmeanarr['time.hour'] == HOUR].squeeze()
    
                    mmeanarr = filtering(mmeanarr, outv, outv)
                    per_year.append(mmeanarr)
                print('Saved year', yy)
                print('Cases in clim mean:', len(per_year))
                mean_arr.append(xr.concat(per_year, dim='time').mean('time'))
                del meanarr
                del per_year
        try:
            fullmean = xr.concat(mean_arr, dim='years').mean('years')
        except:
            print(dt_main, 'Concatenation failed, continue')
        fullmean = fullmean.to_dataset()
        comp = dict(zlib=True, complevel=5)
        encoding = {var: comp for var in fullmean.data_vars}
        
        fullmean.to_netcdf(path=savefile, mode='w', encoding=encoding)
        print('Saved ' + savefile)
        del fullmean
        del mean_arr


### Inputs:
##### Provide hour when you run script!!! as sys.argv[1]
fdir = str(sys.argv[2])
if fdir == 'hist':
    FTAG = 'historical'
else:
    FTAG = 'future'

main = '/home/users/cornkle/linked_CP4/'
main_lmcs = '/home/users/cornkle/lmcs/cklein/CP_models/MCS_files/'
data_path = main + '/'+fdir
out_path = main_lmcs + 'climatology/'+FTAG+'/'
box = [-19, 30, 4, 26]  # W- E , S - N geographical coordinates box

years = np.array(np.arange(2003,2004), dtype=str)
months = (['07', '08', '09'])
days = np.array(np.arange(1,31), dtype=str)

HOUR= int(sys.argv[1]) # hour to extract


vars = ['pblH','sh', 'lh', 't2', 'q2', 'SM', 'u10', 'v10', 't_pl','colWetMass', 'colDryMass', 'lw_out_PBLtop','lsRain','t_pl', 'u_pl', 'v_pl', 'q_pl', 'omega_pl', 'geoH_pl']
datelist = []
for y,m,d in itertools.product(years, months, days):
    datelist.append(y+m+str(d).zfill(2))

for d in datelist:

    print('Doing ', d)
    file_save(data_path, out_path, vars, d, box)

