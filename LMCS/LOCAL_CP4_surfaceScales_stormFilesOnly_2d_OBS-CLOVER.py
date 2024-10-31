import xarray as xr
import glob
import os
import numpy as np
import ipdb
import datetime
import sys
import pickle as pkl
import sys
#sys.path.append('/home/users/cornkle/pythonWorkspace/')
from land_wavelet import wclass
from utils import constants as cnst, u_arrays
from wavelet import util, wav
import pandas as pd

var = 'LSTA'
hour = '17'
box = [-9, 9, 12, 17]

hist = sorted(glob.glob(cnst.lmcs_drive + '/MCS_5000km2_tables/ERA5_added_12LT/WAf/*.csv'))
hist = sorted(glob.glob('/home/ck/DIR/cornkle/data/CLOVER/CLOVER/saves/bulk_-50_5000km2_GPM_5-26N_16W17E_p15_ERA0.7_TCWV_hourly_wDate*.p'))

strct_values = {'OBS' : []}
strct_scale = {'OBS' : []}
strct_power = {'OBS' : []}

strct_values_y = {'OBS' : []}
strct_scale_y = {'OBS' : []}
strct_power_y = {'OBS' : []}

strct_power_2d = {'power' : [], 'coeffs' : [], 'lsta':[]}

tag = 'posOnly'

flsta = '/home/ck/DIR/cornkle/data/OBS/MSG_LSTA/lsta_netcdf_1330/'



for idx, dats in enumerate([hist]):
    u_dates = []

    for sf in dats:  ########restricted files for testing  [0:5000]
        #ipdb.set_trace()
        fname = os.path.basename(sf)
        year = fname[-6:-2]
        if int(year) not in list(np.arange(2006,2011)):
            continue
        print('Doing', year)
        try:
            sda = pkl.load(open(glob.glob(sf)[0], 'rb'))
        except:
            ipdb.set_trace()

        sda = pd.DataFrame(sda)



        #sda = sda[(sda['lt_hour'] == 17) & ((sda['month']>=6) & (sda['month']<=9)) & ((sda['tminlon']>=box[0]) & (sda['tminlon']<=box[1])) & ((sda['tminlat']>=box[2]) & (sda['tminlat']<=box[3]))]

        mask = (np.array(sda['hour']) == 17) & ((sda['month']>=6) & (sda['month']<=9)) & ((np.array(sda['clon']) >= box[0]) & (np.array(sda['clon']) <= box[1])) & ((np.array(sda['clat']) >= box[2]) & (np.array(sda['clat']) <= box[3]))
        sda = sda[mask]
        utdate = []
        for yy, mm, dd in zip(sda['year'], sda['month'], sda['day']):
            utdate.append(str(yy)+'-'+str(mm).zfill(2)+'-'+str(dd).zfill(2))
        dunique = np.unique(utdate)

        for date in dunique:

            tdate = date[0:4]+date[5:7]+date[8:10]
            try:
                lsta = xr.open_dataset(flsta+'lsta_daily_'+tdate+'.nc')
            except:
                print('LSTA date missing, continue', tdate)
                continue
            lsta = (lsta['LSTA'] - lsta['LSTA'].mean()).squeeze()

            #lsta = lsta.where(np.isfinite(lsta), other=0)
            mask = np.array(utdate)==date
           # ipdb.set_trace()
            mcs_dates = sda[mask]

            for lon, lat in zip(mcs_dates['clon'], mcs_dates['clat']):


                pos = lsta.sel(lat=lat, lon=lon, method='nearest', tolerance=0.02)
                lonpos = np.where(lsta.lon.values==float(pos.lon.values))
                latpos = np.where(lsta.lat.values==float(pos.lat.values))

                # distx = 57  # 57 = 250 km at 4.4 res, 500km across
                # disty = 57
                #ipdb.set_trace()
                filt = u_arrays.cut_kernel(lsta.values, int(lonpos[0]), int(latpos[0]), 81)
                filt = filt - np.nanmean(filt)
                # ipdb.set_trace()
                # try:
                #
                #     filt = lsta.sel(lat=slice(pos.lat - 2.25, pos.lat + 2.25),
                #                      lon=slice(pos.lon - 2.25, pos.lon + 2.25)).squeeze()
                # except IndexError:
                #     continue

                dist = 163
                midpoint = 82
                dy = 16

                if (filt.shape[0] != dist) | (filt.shape[0] != dist):
                    print(filt, 'wrong dimension box')
                    ipdb.set_trace()
                    continue

                if np.sum(np.isnan(filt)) >= filt.size*0.5:
                    print('Too many missing values, continue')
                    continue

                filt = xr.DataArray(data=filt, coords={'lon': np.arange(dist), 'lat': np.arange(dist)}, dims=['lat', 'lon'])

                lines = filt.isel(lat=slice(midpoint - dy, midpoint + dy)).mean('lat')
                linesy = filt.isel(lon=slice(midpoint - dy, midpoint + dy)).mean('lon')



                cp4_mcs_wav = filt.where(np.isfinite(filt), other=0)

                wObj = wclass.landwav('CP4_OBS_VARS')

                wObj.read_img(cp4_mcs_wav.values, cp4_mcs_wav.lon.values, cp4_mcs_wav.lat.values)
                coeffs, power, scales, period = wObj.applyWavelet(normed='scale_squared', le_thresh=0, fill=0)

                # pos = cp4_mcs_wav==0
                # coeffs[:,pos] = np.nan
                # power[:, pos] = np.nan

                coeffline = np.nanmean(coeffs[:,midpoint-dy:midpoint+dy,:], axis=1)
                powerline = np.nanmean(power[:,midpoint-dy:midpoint+dy,:], axis=1)

                coeffliney = np.nanmean(coeffs[:,:,midpoint-dy:midpoint+dy], axis=2)
                powerliney = np.nanmean(power[:,:, midpoint-dy:midpoint+dy], axis=2)
                #ipdb.set_trace()
                # CONTINUE DEVELOPING CODE HERE!!
                strct_values['OBS'].append(lines)
                strct_scale['OBS'].append(coeffline)
                strct_power['OBS'].append(powerline)

                strct_values_y['OBS'].append(linesy)
                strct_scale_y['OBS'].append(coeffliney)
                strct_power_y['OBS'].append(powerliney)

                strct_power_2d['power'].append(power)
                strct_power_2d['coeffs'].append(coeffs)
                strct_power_2d['lsta'].append(filt.values)

            del lsta

#
strct_scale['scales'] = scales
strct_power['scales'] = scales
pkl.dump(strct_values, open(cnst.network_data+ 'data/LMCS/CP4_study_saves/'+var+'_xcross_-11to11_'+hour+'_2d_x'+tag+'.p','wb'))
pkl.dump(strct_scale, open(cnst.network_data+ 'data/LMCS/CP4_study_saves/'+var+'_wcoeffs_-11to11_'+hour+'_2d_x'+tag+'.p','wb'))
pkl.dump(strct_power, open(cnst.network_data+ 'data/LMCS/CP4_study_saves/'+var+'_power_-11to11_'+hour+'_2d_x'+tag+'.p','wb'))

strct_scale_y['scales'] = scales
strct_power_y['scales'] = scales
pkl.dump(strct_values_y, open(cnst.network_data+ 'data/LMCS/CP4_study_saves/'+var+'_xcross_-11to11_'+hour+'_2d_y'+tag+'.p','wb'))
pkl.dump(strct_scale_y, open(cnst.network_data+ 'data/LMCS/CP4_study_saves/'+var+'_wcoeffs_-11to11_'+hour+'_2d_y'+tag+'.p','wb'))
pkl.dump(strct_power_y, open(cnst.network_data+ 'data/LMCS/CP4_study_saves/'+var+'_power_-11to11_'+hour+'_2d_y'+tag+'.p','wb'))

pkl.dump(strct_power_2d, open(cnst.network_data+ 'data/LMCS/CP4_study_saves/'+var+'_power_-11to11_'+hour+'_2d_xy'+tag+'.p','wb'))