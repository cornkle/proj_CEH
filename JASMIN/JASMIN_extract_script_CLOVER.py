import xarray as xr
import glob
import os
import itertools
from JASMIN import constants as cnst, MetUM_variables as mv
import pdb


#veg = '/users/global/cornkle/w2018_bamba/qrparm.cci.4km.nc'

### 2d vars , xmh*.pc*.nc files
folder = '/users/global/cornkle/figs/CLOVER/CP4testfiles'

local_box = [-18+360,14+360, 3.5, 14]
temp_box = [-18+360,48+360, -36, 32]

months = [3,5] # March-May

dic = {

    't2_daily' : ([temp_box], ['keep'], []),
    'lw_out_PBLtop' : ([local_box], ['keep'], []),
    'u_pl' : ([temp_box, local_box], ['agg', 'keep'], [650, 850]),
    't_pl' : ([local_box], ['agg'], [925]),
    'omega_pl' : ([local_box], ['keep'], []),
    'lsRain' : ([local_box], ['keep'], []),

}
keys = dic.keys()

for k in keys:

    info = cnst.VARDIC[k]
    dinfo = dic[k]
    var = mv.create_CP4_filename(k)
    infolder = folder+os.sep + var
    outfolder = folder +os.sep +  k
    files = glob.glob(infolder + os.sep + var+'*.nc' )

    for f in files:

        ds = xr.open_dataset(f)

        if (ds['time.month'][0]<months[0]) | (ds['time.month'][0]>months[1]):
            continue

        fname = os.path.basename(f)
        box = dinfo[0]

        for id, b in enumerate(box):

            agg = dinfo[1][id]
            pres = dinfo[2]
            cut = ds.sel(longitude=slice(b[0], b[1]), latitude=slice(b[2], b[3]))
            da = cut[var]

            if pres != []:
                da = da.sel(pressure=pres)

            # if agg != 'keep':
            #     da = da.resample('d', dim='time', how='mean')

            outname = fname.replace(var, k)

            comp = dict(zlib=True, complevel=5)

            da.name = k
            #encoding = {var: comp for var in da.data_vars}
            encoding = {k: {'complevel': 5, 'zlib': True}}
            if not os.path.exists(outfolder):
                os.makedirs(outfolder)

            da.to_netcdf(outfolder + os.sep + outname , format='NETCDF4', encoding=encoding)
            da.close()



