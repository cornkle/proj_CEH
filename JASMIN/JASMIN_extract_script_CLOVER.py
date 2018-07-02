import xarray as xr
import glob
import os
import itertools
from JASMIN import constants as cnst, MetUM_variables as mv
import numpy as np
import pdb
fpath = '/users/global/cornkle/figs/CLOVER/CP4testfiles'

def run(fpath, orig_names=False):


    local_box = [-18+360,14+360, 3.5, 14]
    temp_box = [-18+360,35+360, 3.5, 30]

    months = [3,5] # March-May
    month_tag = 'MAM'

    dic = {

        't2' : ([temp_box], ['keep'], [], [12,3]),
        'lw_out_PBLtop' : ([local_box], ['keep'], [], []),
        'u_pl' : ([temp_box], ['keep', 'keep'], [650, 850], [12,3]),
        't_pl' : ([temp_box], ['keep'], [925], [12,3]),
        'omega_pl' : ([local_box], ['keep'], [650,300], []),
        'lsRain' : ([local_box], ['keep'], [], []),

    }
    keys = dic.keys()

    for k in keys:

        info = cnst.VARDIC[k]
        dinfo = dic[k]
        if orig_names:
            var = mv.create_CP4_filename(k)
        else:
            var = k
        infolder = fpath+os.sep + var
        outfolder = fpath +os.sep +  k
        files = glob.glob(infolder + os.sep + var+'*.nc' )

        for f in files:

            ds = xr.open_dataset(f)

            if (ds['time.month'][0]<months[0]) | (ds['time.month'][0]>months[1]):
                continue

            if dinfo[3] != []:
                ds = ds.isel(time=[np.in1d(ds['time.hour'].values, dinfo[3])][0])
                pdb.set_trace()

            fname = os.path.basename(f)
            box = dinfo[0]

            for id, b in enumerate(box):

                agg = dinfo[1][id]
                pres = dinfo[2]
                cut = ds.sel(longitude=slice(b[0], b[1]), latitude=slice(b[2], b[3]))
                da = cut[var]

                if pres != []:
                    da = da.sel(pressure=pres)

                if agg != 'keep':
                    da = da.resample('24H', base=16, dim='time', skipna=True, how='mean')

                outname = fname.replace(var, k+'_'+month_tag)

                comp = dict(zlib=True, complevel=5)

                da.name = k
                da.longitude.values = da.longitude.values-360
                #encoding = {var: comp for var in da.data_vars}
                encoding = {k: {'complevel': 5, 'zlib': True}}
                if not os.path.exists(outfolder):
                    os.makedirs(outfolder)

                da.to_netcdf(outfolder + os.sep + outname , format='NETCDF4', encoding=encoding)
                da.close()

                print('Wrote '+ outfolder + os.sep + outname)


