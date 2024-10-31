import xarray as xr
import glob
import os
import numpy as np
import pdb
import datetime
import sys
import pickle as pkl
import sys
sys.path.append('/home/users/cornkle/pythonWorkspace/')
from LMCS.land_wavelet import wclass

var = sys.argv[1]
hour = '17'
hist = sorted(glob.glob('/home/users/cornkle/lmcs/cklein/CP_models/MCS_files/moved_fromCEH/CP4_box/CP4_allHours_historical_5000km2_-50_WAf_box_v2/*_'+hour+'*.nc'))
fut = sorted(glob.glob('/home/users/cornkle/lmcs/cklein/CP_models/MCS_files/moved_fromCEH/CP4_box/CP4_allHours_future_5000km2_-50_WAf_box_v2/*_'+hour+'*.nc'))

tags = ['CP4hist', 'CP4fut']

strct = {'CP4hist' : [], 'CP4fut' : []}
strct_anom = {'CP4hist' : [], 'CP4fut' : []}

for idx, dats in enumerate([hist, fut]):
    u_dates = []
    cp4_full_path = '/home/users/cornkle/CP4home/'+tags[idx]+'/'
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

        cp4_agg = xr.open_mfdataset(files_toread)
        cp4_agg = cp4_agg.sel(time=cp4_agg['time.hour']==12)
        
            
        for sf in storm_files:

            sda = xr.open_dataset(sf)
            lon = sda.minlon
            lat = sda.minlat

            cp4_latlon = cp4_agg.sel(longitude=lon, latitude=lat, method='nearest')
            idx = np.where(cp4_agg.longitude == cp4_latlon.longitude)
            idy = np.where(cp4_agg.latitude == cp4_latlon.latitude)

            distx = 57  # 57 = 250 km at 4.4 res, 500km across
            disty = 57
            try:
                cp4_box = cp4_agg.isel(latitude=slice(ypos - disty, ypos + disty + 1),
                             longitude=slice(xpos - distx, xpos + distx + 1))
            except IndexError:
            continue
            if (len(cp4_box.latitude) != disty * 2 + 1) | (len(cp4_box.longitude) != distx * 2 + 1):
                print(cp4_box)
            continue
            cp4_box = cp4_box.where((cp4_box['lsRain']*3600 < 0.005) & cp4_box['lw_out_PBLtop']/100 >=-30)

            cp4_mcs_wav = sda[var] - sda[var].mean()

            cp4_box_wav = cp4_box[var] - cp4_box[var].mean()

            

            

            
            
            cp4_mean = cp4_box.groupby('time.hour').mean('time')
            cp4_anom = cp4_box.groupby('time.hour')-cp4_mean

            strct[tags[idx]].append(cp4_box.values)
            strct_anom[tags[idx]].append(cp4_anom.values)

pkl.dump(strct, open('/home/users/cornkle/lmcs/cklein/CP_models/MCS_files/'+var+'_timeseries_-10to+10days_'+hour+'.p','wb'))
pkl.dump(strct_anom, open('/home/users/cornkle/lmcs/cklein/CP_models/MCS_files/'+var+'_anom_timeseries_-10to+10days_'+hour+'.p','wb'))