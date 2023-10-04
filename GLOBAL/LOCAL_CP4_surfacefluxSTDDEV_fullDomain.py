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
hour = '12'

#for yy in range(1997,2006):
hist = sorted(glob.glob(cnst.other_drive + '/CP4/CP4_WestAfrica/CP4hist/'+var+'/'+var+'*.nc')) #+'_'+str(yy)+'*.nc'))
fut = sorted(glob.glob(cnst.other_drive + '/CP4/CP4_WestAfrica/CP4fut/'+var+'/'+var+'*.nc')) #+'_'+str(yy)+'*.nc'))

tags = ['CP4hist', 'CP4fut']

strct_std = {'CP4hist' : [], 'CP4fut' : []}

for idx, dats in enumerate([hist, fut]):
    u_dates = []
    for idp, sf in enumerate(dats):  ########restricted files for testing
        fname = os.path.basename(sf).split('_')[-1]
        u_date = fname[0:4]+fname[4:6]+fname[6:8]

        if int(fname[4:6]) not in [8]:
            continue
        if int(fname[6:8]) not in [15]:
            continue

        print('Doing', u_date)

        # sda = xr.open_dataset(sf)
        # lines = sda[var].sel(time=(sda['time.hour']==12), latitude=slice(8,15), longitude=slice(-15,25)).squeeze()
        sdc = xr.open_mfdataset(dats[idp-15:idp+15])
        lines_clim = sdc[var].sel(time=(sdc['time.hour'] == 12), latitude=slice(8, 15), longitude=slice(-15, 25)).std('time').load()

        cp4_mcs_wav = lines_clim.values.flatten()#- lines_clim)

        strct_std[tags[idx]].extend(cp4_mcs_wav)


pkl.dump(strct_std, open(cnst.network_data+ 'data/LMCS/CP4_study_saves/'+var+'_stddev_'+hour+'_Aug.p','wb'))