import xarray as xr
import glob
import os
import numpy as np
import pdb
import datetime
import sys
import pickle as pkl
from utils import constants as cnst

var = sys.argv[1]
hour = '17'
hist = sorted(glob.glob(cnst.lmcs_drive + '/CP_models/MCS_files/MODELS/CP4_box/CP4_allHours_historical_5000km2_-50_WAf_box_v2/*_'+hour+'*.nc'))
fut = sorted(glob.glob(cnst.lmcs_drive + '/CP_models/MCS_files/MODELS/CP4_box/CP4_allHours_future_5000km2_-50_WAf_box_v2/*_'+hour+'*.nc'))

tags = ['CP4hist', 'CP4fut']

strct = {'CP4hist' : [], 'CP4fut' : []}
strct_anom = {'CP4hist' : [], 'CP4fut' : []}

for idx, dats in enumerate([hist, fut]):
    u_dates = []
    cp4_full_path = cnst.other_drive + 'CP4/CP4_WestAfrica/'+tags[idx]+'/'
    cp4_all_files = sorted(glob.glob(cp4_full_path + var + '/' + var +'_*.nc'))
    
    for file in dats[0:100]:  ########restricted files for testing
        fname = os.path.basename(file)
        u_dates.append(fname[0:4]+fname[5:7]+fname[8:10])
    loop_dates = np.unique(u_dates)
    for ld in loop_dates:
        storm_files = []
        for file in dats:
            if ld[0:4]+'-'+ld[4:6]+'-'+ld[6:8] in file:
                storm_files.append(file)

        
        try:
            cp4_file = glob.glob(cp4_full_path + var + '/' + var +'_*_'+ld+'*.nc')[0]
        except:
            print('Couldnt open file')
            pdb.set_trace()
            continue

        cp4_pos = cp4_all_files.index(cp4_file)

        files_toread = np.array(cp4_all_files)[cp4_pos-10:cp4_pos+10]

        cp4_agg = xr.open_mfdataset(files_toread)[var]
            
        for sf in storm_files:

            sda = xr.open_dataset(sf)
            lon = sda.minlon
            lat = sda.minlat

            cp4_box = cp4_agg.sel(longitude=slice(lon-0.35, lon+0.34), latitude=slice(lat-0.35, lat+0.35)).mean(['latitude','longitude'])
            cp4_mean = cp4_box.groupby('time.hour').mean('time')
            cp4_anom = cp4_box.groupby('time.hour')-cp4_mean

            strct[tags[idx]].append(cp4_box.values)
            strct_anom[tags[idx]].append(cp4_anom.values)

pkl.dump(strct, open(cnst.network_data+ 'data/LMCS/CP4_study_saves/'+var+'_timeseries_-10to+10days_'+hour+'.p','wb'))
pkl.dump(strct_anom, open(cnst.network_data+ 'data/LMCS/CP4_study_saves/'+var+'_anom_timeseries_-10to+10days_'+hour+'.p','wb'))