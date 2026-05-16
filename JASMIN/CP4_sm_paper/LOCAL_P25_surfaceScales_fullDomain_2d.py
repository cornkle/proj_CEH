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
from utils import constants as cnst
import matplotlib.pyplot as plt

var = sys.argv[1]
sign = sys.argv[2]
norming = sys.argv[3] # 'stddev' or 'scale_squared'
hour = '17'

#for yy in range(1997,2006):
hist = sorted(glob.glob(cnst.other_drive + '/CP4/CP4_WestAfrica/CP25hist/'+var+'/'+var+'*.nc')) #+'_'+str(yy)+'*.nc'))
fut = sorted(glob.glob(cnst.other_drive + '/CP4/CP4_WestAfrica/CP25fut/'+var+'/'+var+'*.nc')) #+'_'+str(yy)+'*.nc'))

tags = ['CP4hist', 'CP4fut']

strct_coeff = {'CP4hist' : [], 'CP4fut' : []}
strct_power = {'CP4hist' : [], 'CP4fut' : []}

if sign == 'pos':
    tag = 'posOnly'
if sign == 'neg':
    tag = 'negOnly'
if sign == 'both':
    tag = 'both'

for idx, dats in enumerate([hist, fut]):
    print('DOING', tags[idx])
    u_dates = []
    for idp, sf in enumerate(dats):  ########restricted files for testing
        fname = os.path.basename(sf).split('_')[-1]
        u_date = fname[0:4]+fname[4:6]+fname[6:8]

        if int(fname[4:6]) not in [8]:
            continue
        # if int(fname[6:8]) not in [1,5,10,15,20,25]:
        #     continue

        print('Doing', u_date)

        sda_all = xr.open_dataset(sf)[var]

        for sda in sda_all:



            if sda['time.day'] not in [1,5,10,15,20,25]:
                continue

            print('Doing', sda.time)

            if 'depth' in sda.coords:
                sda = sda.isel(depth=0).squeeze()
                sda= sda.where(sda < 500, other=np.nan)

            try:
                lines = sda.sel(time=(sda['time.hour']==12), latitude=slice(5,25), longitude=slice(-15,25)).squeeze()
            except:
                lines = sda.sel(latitude=slice(5, 25), longitude=slice(5, 25)).squeeze()
            lines = lines#-lines.mean('longitude')
            cp4_mcs_wav = (lines )
            cp4_mcs_wav = cp4_mcs_wav.where(np.isfinite(cp4_mcs_wav), other=0)

            wObj = wclass.landwav('P25_VARS')

            wObj.read_img(cp4_mcs_wav.values, cp4_mcs_wav.longitude.values, cp4_mcs_wav.latitude.values)
            if sign == 'pos':
                coeffs, power, scales, period = wObj.applyWavelet(normed=norming, le_thresh=0, fill=0)
            if sign == 'neg':
                coeffs, power, scales, period = wObj.applyWavelet(normed=norming,ge_thresh=0, fill=0)
            if sign == 'both':
                coeffs, power, scales, period = wObj.applyWavelet(normed=norming, ge_thresh=None, fill=0)


            ds = xr.Dataset()
            ds['power'] =xr.DataArray(data=power, coords={'scales':scales, 'longitude':cp4_mcs_wav.longitude, 'latitude':cp4_mcs_wav.latitude}, dims=['scales','latitude','longitude'])
            ds['coeffs'] = xr.DataArray(data=coeffs, coords={'scales':scales, 'longitude':cp4_mcs_wav.longitude, 'latitude':cp4_mcs_wav.latitude}, dims=['scales','latitude','longitude'])



            ds = ds.sel(latitude=slice(8,16), longitude=slice(-12,20))#latitude=slice(9,17), longitude=slice(-12,12))

            # f = plt.figure()
            # ax = f.add_subplot(121)
            # plt.pcolormesh(cp4_mcs_wav)
            # plt.colorbar()
            # ax = f.add_subplot(122)
            # plt.pcolormesh(ds['coeffs'].isel(scales=slice(0, 3)).mean('scales'))
            # plt.colorbar()
            # plt.show()

            ds = ds.groupby('scales').mean(['latitude','longitude'])


            # CONTINUE DEVELOPING CODE HERE!!

            strct_coeff[tags[idx]].append(ds['coeffs'].values)
            strct_power[tags[idx]].append(ds['power'].values)

        del sda_all
        del wObj
        del ds


strct_coeff['scales'] = scales
strct_power['scales'] = scales

pkl.dump(strct_coeff, open(cnst.network_data+ 'data/LMCS/CP4_study_saves/'+var+'_wcoeffs_fullDomain_'+hour+'_'+tag+'_'+norming+'_17N_P25_south.p','wb')) #'_'+str(yy)+'.p','wb'))
pkl.dump(strct_power, open(cnst.network_data+ 'data/LMCS/CP4_study_saves/'+var+'_power_fullDomain_'+hour+'_'+tag+'_'+norming+'_17N_P25_south.p','wb')) #'_'+str(yy)+'.p','wb'))