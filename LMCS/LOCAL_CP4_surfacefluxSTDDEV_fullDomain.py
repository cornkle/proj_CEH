import xarray as xr
import glob
import os
import numpy as np
import ipdb
import datetime
import sys
import pickle as pkl
import sys
import salem
#sys.path.append('/home/users/cornkle/pythonWorkspace/')
from land_wavelet import wclass
from utils import constants as cnst
import matplotlib.pyplot as plt

var = sys.argv[1]
hour = '12'

#for yy in range(1997,2006):
hist = sorted(glob.glob(cnst.other_drive + '/CP4/CP4_WestAfrica/CP4hist/'+var+'/'+var+'*.nc')) #+'_'+str(yy)+'*.nc'))
fut = sorted(glob.glob(cnst.other_drive + '/CP4/CP4_WestAfrica/CP4fut/'+var+'/'+var+'*.nc')) #+'_'+str(yy)+'*.nc'))

dummy = xr.open_dataarray(glob.glob(cnst.other_drive + '/CP4/CP4_WestAfrica/CP25fut/t2/t2*.nc')[0]).isel(time=0)
dummy2 = xr.open_dataarray(glob.glob(cnst.other_drive + '/CP4/CP4_WestAfrica/CP4fut/t2/t2*.nc')[0]).isel(time=0).sel(latitude=slice(8, 15), longitude=slice(-12, 20))

regrid, lut = dummy.salem.lookup_transform(dummy2, return_lut=True)

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
        try:
            sdc = xr.open_mfdataset(dats[idp-15:idp+15])
            if 'depth' in sdc.coords:
                sdc = sdc.isel(depth=0)
                sdc = sdc['SM'].where(sdc['SM'] < 500, other=np.nan)
            else:
                sdc = sdc[var]
        except:
            continue

        lines_clim = sdc.sel(time=(sdc['time.hour'] == 12), latitude=slice(8, 15), longitude=slice(-12, 20)).load()
        #ipdb.set_trace()
        lines_clim_25 = dummy.salem.lookup_transform(lines_clim, lut=lut)
        lines_clim_25 = lines_clim_25.std('time')

        cp4_mcs_wav = lines_clim_25.values.flatten()#- lines_clim)

        strct_std[tags[idx]].extend(cp4_mcs_wav[np.isfinite(cp4_mcs_wav)])


pkl.dump(strct_std, open(cnst.network_data+ 'data/LMCS/CP4_study_saves/'+var+'_stddev_'+hour+'_Aug_CP4_on25.p','wb'))