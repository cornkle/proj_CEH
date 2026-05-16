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
from wavelet import util, wav

var = sys.argv[1]
hour = '17'
hist = sorted(glob.glob(cnst.lmcs_drive + '/CP_models/MCS_files/MODELS/CP4_box_anom/CP4_allHours_historical_5000km2_-50_WAf_box_anom_v2/*_'+hour+'*.nc'))
fut = sorted(glob.glob(cnst.lmcs_drive + '/CP_models/MCS_files/MODELS/CP4_box_anom/CP4_allHours_future_5000km2_-50_WAf_box_anom_v2/*_'+hour+'*.nc'))

tags = ['CP4hist', 'CP4fut']

strct_values = {'CP4hist' : [], 'CP4fut' : []}
strct_scale = {'CP4hist' : [], 'CP4fut' : []}
strct_power = {'CP4hist' : [], 'CP4fut' : []}

for idx, dats in enumerate([hist, fut]):
    u_dates = []

    for sf in dats[0:5000]:  ########restricted files for testing  [0:5000]
        #ipdb.set_trace()
        fname = os.path.basename(sf)
        u_date = fname[0:4]+fname[5:7]+fname[8:10]

        sda = xr.open_dataset(sf)
        midpoint = 58
        dy = 11
        lines = sda[var].isel(latitude=slice(midpoint-dy, midpoint+dy)).mean('latitude')

        cp4_mcs_wav = sda[var].where(np.isfinite(lines), other=0)# - np.nanmean(sda[var])
        # inter = np.where(np.isnan(cp4_mcs_wav))
        # cp4_mcs_wav[inter] = 0.01

        wObj = wclass.landwav('CP4_VARS')

        wObj.read_img(cp4_mcs_wav.values, cp4_mcs_wav.longitude.values, cp4_mcs_wav.latitude.values)
        coeffs, power, scales, period = wObj.applyWavelet_1d(normed='scale', direction='x', le_thresh=0, fill=0)

        coeffline = np.nanmean(coeffs[:,midpoint-dy:midpoint+dy,:], axis=1)
        powerline = np.nanmean(power[:,midpoint-dy:midpoint+dy,:], axis=1)
        #ipdb.set_trace()
        # CONTINUE DEVELOPING CODE HERE!!
        strct_values[tags[idx]].append(lines)
        strct_scale[tags[idx]].append(coeffline)
        strct_power[tags[idx]].append(powerline)

        del sda

    strct_scale['scales'] = scales
    strct_power['scales'] = scales
pkl.dump(strct_values, open(cnst.network_data+ 'data/LMCS/CP4_study_saves/'+var+'_xcross_-11to11_'+hour+'_1d_x.p','wb'))
pkl.dump(strct_scale, open(cnst.network_data+ 'data/LMCS/CP4_study_saves/'+var+'_wcoeffs_-11to11_'+hour+'_1d_x.p','wb'))
pkl.dump(strct_power, open(cnst.network_data+ 'data/LMCS/CP4_study_saves/'+var+'_power_-11to11_'+hour+'_1d_x.p','wb'))