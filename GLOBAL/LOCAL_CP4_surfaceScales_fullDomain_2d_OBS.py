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
hour = '17'

#for yy in range(1997,2006):
hist = sorted(glob.glob('/home/ck/DIR/cornkle/data/OBS/MSG_LSTA/lsta_netcdf_1330/*.nc'))
tags = ['OBS']

strct_coeff = {'OBS' : []}
strct_power = {'OBS' : []}

tag = 'negOnly'

for idx, dats in enumerate([hist]):
    u_dates = []
    for idp, sf in enumerate(dats):  ########restricted files for testing
        fname = os.path.basename(sf).split('_')[-1]
        #ipdb.set_trace()
        u_date = fname[0:4]+fname[4:6]+fname[6:8]

        if int(fname[4:6]) not in [8]:
            continue
        if int(fname[6:8]) not in np.arange(1,31,3):
            continue

        print('Doing', u_date)

        sda = xr.open_dataset(sf).squeeze()
        lines = sda[var]
        #sdc = xr.open_mfdataset(dats[idp-10:idp+10])
        #lines_clim = sdc[var].sel(time=(sdc['time.hour'] == 12), lat=slice(7, 25), lon=slice(-15, 25)).mean('time').load()

        cp4_mcs_wav = (lines )#- lines_clim)

        cp4_mcs_wav = cp4_mcs_wav.where(np.isfinite(cp4_mcs_wav), other=0)#-#np.nanmean(cp4_mcs_wav)nter
        # inter = np.where(np.isnan(cp4_mcs_wav))
        # cp4_mcs_wav[inter] = 0
        #ipdb.set_trace()
        wObj = wclass.landwav('CP4_OBS_VARS')

        wObj.read_img(cp4_mcs_wav.values, cp4_mcs_wav.lon.values, cp4_mcs_wav.lat.values)
        coeffs, power, scales, period = wObj.applyWavelet(normed='scale_squared', ge_thresh=0, fill=0)


        ds = xr.Dataset()
        ds['power'] =xr.DataArray(data=power, coords={'scales':scales, 'lon':cp4_mcs_wav.lon, 'lat':cp4_mcs_wav.lat}, dims=['scales','lat','lon'])
        ds['coeffs'] = xr.DataArray(data=coeffs, coords={'scales':scales, 'lon':cp4_mcs_wav.lon, 'lat':cp4_mcs_wav.lat}, dims=['scales','lat','lon'])



        ds = ds.sel(lat=slice(12,17), lon=slice(-8,14))#latitude=slice(9,17), longitude=slice(-12,12))

        # f = plt.figure()
        # ax = f.add_subplot(121)
        # plt.pcolormesh(cp4_mcs_wav)
        # plt.colorbar()
        # ax = f.add_subplot(122)
        # plt.pcolormesh(ds['coeffs'].isel(scales=slice(0, 3)).mean('scales'))
        # plt.colorbar()
        # plt.show()

        ds = ds.groupby('scales').mean(['lat','lon'])


        # CONTINUE DEVELOPING CODE HERE!!

        strct_coeff[tags[idx]].append(ds['coeffs'].values)
        strct_power[tags[idx]].append(ds['power'].values)

        del sda
        del wObj
        del ds

strct_coeff['scales'] = scales
strct_power['scales'] = scales

pkl.dump(strct_coeff, open(cnst.network_data+ 'data/LMCS/CP4_study_saves/OBS_'+var+'_wcoeffs_fullDomain_'+hour+'_'+tag+'.p','wb')) #'_'+str(yy)+'.p','wb'))
pkl.dump(strct_power, open(cnst.network_data+ 'data/LMCS/CP4_study_saves/OBS_'+var+'_power_fullDomain_'+hour+'_'+tag+'.p','wb')) #'_'+str(yy)+'.p','wb'))