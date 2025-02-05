import xarray as xr
import glob
import os
import numpy as np
import pdb
import datetime
import sys
import pickle as pkl
import ipdb
import pandas as pd

sys.path.append('/home/users/cornkle/pythonWorkspace/')

from JASMIN import MetUM_variables as mu

ivar = str(sys.argv[1])
var = mu.create_CP4_filename(ivar)


hour = '17'
hist = sorted(glob.glob('/gws/nopw/j04/lmcs/cklein/CP_models/MCS_files/CP4_box_mean_JASMIN/CP4_historical_5000km2_-50_box_mean_v3/*-08-*_'+hour+':*.nc'))
fut = sorted(glob.glob('/gws/nopw/j04/lmcs/cklein/CP_models/MCS_files/CP4_box_mean_JASMIN/CP4_future_5000km2_-50_box_mean_v3/*-08-*_'+hour+':*.nc'))

tags = ['hist', 'future']

strct = {'hist' : [], 'future' : []}
strct_anom = {'hist' : [], 'future' : []}

for idx, dats in enumerate([hist, fut]):
    u_dates = []
    cp4_full_path = '/home/users/cornkle/linked_CP4/'+tags[idx]+'/'
    #cp4_all_files = sorted(glob.glob(cp4_full_path + ivar + '/' + var +'_*.nc'))
    
    for file in dats[400:500]:  ########restricted files for testing
        fname = os.path.basename(file)
        u_dates.append(fname[0:4]+fname[5:7]+fname[8:10])
    loop_dates = np.unique(u_dates)
    
    for ld in loop_dates:
        print('Doing', ld)
        storm_files = []
        for file in dats:
            if ld[0:4]+'-'+ld[4:6]+'-'+ld[6:8] in file:
                storm_files.append(file)

        try:
            cp4_file = glob.glob(cp4_full_path + ivar + '/' + var +'_*_'+ld+'*.nc')[0]
        except:
            print('Couldnt open file')
            ipdb.set_trace()
            continue

        # cp4_pos = cp4_all_files.index(cp4_file)

        # files_toread = np.array(cp4_all_files)[cp4_pos-10:cp4_pos+10]
    
        files_toread = []
        dt = datetime.datetime(int(ld[0:4]), int(ld[4:6]), int(ld[6:8]))
        for dd in np.arange(-7,8):
            ndate = dt + pd.Timedelta(str(dd)+'days')
            if ndate.day >30:
                dday = 30
            else:
                dday = ndate.day
           
            ndatestring = str(ndate.year)+str(ndate.month).zfill(2)+str(dday).zfill(2)
          
            getf = glob.glob(cp4_full_path + ivar + '/' + var +'_*_'+ndatestring+'*.nc')[0]
            files_toread.append(getf)

        print(len(files_toread))

        if len(files_toread) < 15:
            continue

        cp4_agg = xr.open_mfdataset(files_toread)
      
        cp4_agg=cp4_agg.assign_coords({"longitude":(cp4_agg.longitude-360)})
        #cp4_agg = cp4_agg.sel(time=cp4_agg['time.hour']==12)
        
            
        for sf in storm_files:

            sda = xr.open_dataset(sf)
            
            if 'pressure' in sda.coords:
                sda = sda.sel(pressure=950)
            lon = sda.minlon
            lat = sda.minlat

            cp4_latlon = cp4_agg.sel(longitude=lon, latitude=lat, method='nearest')
            xpos = np.where(cp4_agg.longitude == cp4_latlon.longitude)[0][0]
            ypos = np.where(cp4_agg.latitude == cp4_latlon.latitude)[0][0]

            #ipdb.set_trace()

            distx = 10  # 57 = 250 km at 4.4 res, 500km across
            disty = 10
            try:
                cp4_box = cp4_agg.isel(latitude=slice(ypos - disty, ypos + disty + 1),
                             longitude=slice(xpos - distx, xpos + distx + 1))
            except IndexError:
                continue
            if (len(cp4_box.latitude) != disty * 2 + 1) | (len(cp4_box.longitude) != distx * 2 + 1):
                print(cp4_box)
                continue
            #cp4_box = cp4_box.where((cp4_box[mu.create_CP4_filename('lsRain')]*3600 < 0.005))
            
            cp4_box = cp4_box[var].mean(['longitude','latitude'])
            cp4_mean = cp4_box.groupby('time.hour').mean('time')

            cp4_anom = cp4_box.groupby('time.hour')-cp4_mean
            #ipdb.set_trace()

            strct[tags[idx]].append(cp4_box.values)
            strct_anom[tags[idx]].append(cp4_anom.values)
            
pkl.dump(strct, open('/gws/nopw/j04/lmcs/cklein/CP_models/CP4_timeseries/'+ivar+'_timeseries_-10to+10days_'+hour+'.p','wb'))
pkl.dump(strct_anom, open('/gws/nopw/j04/lmcs/cklein/CP_models/CP4_timeseries/'+ivar+'_anom_timeseries_-10to+10days_'+hour+'.p','wb'))
print('Saved timeseries files!')