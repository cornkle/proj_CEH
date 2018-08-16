import xarray as xr
import glob
import os
import itertools
from JASMIN import constants as cnst, MetUM_variables as mv
import numpy as np
import pdb
from JASMIN import constants

def run(orig_names=False):

    fpath = '/home/users/cornkle/runscript/in'
    outpath = '/home/users/cornkle/runscript/out_full_3d'

    temp_box = [-18+360,35+360, 3.5, 30]

    dic = constants.VARDIC
    keys = dic.keys()

    for k in keys:

        info = cnst.VARDIC[k]
        dinfo = dic[k]
        var = mv.create_CP4_filename(k)

        if not orig_names:
            pathvar = k
        else:
            pathvar = var

        infolder = fpath+os.sep + pathvar
        outfolder = outpath +os.sep +  k
        files = glob.glob(infolder + os.sep + var+'*.nc' )
        for f in files:

            fname = os.path.basename(f)
            outname = fname.replace(var, k)
            outfile = outfolder + os.sep + outname
            if os.path.isfile(outfile):
                print('File already exists, continue.')
                continue
            datestr = fname[-15:-8]

            if (datestr!='20020510') & (datestr!='20020511') & (datestr!='20020512'):
                continue

            ds = xr.open_dataset(f)
            b = temp_box

            cut = ds.sel(longitude=slice(b[0], b[1]), latitude=slice(b[2], b[3]))
            try:
              da = cut[var]
            except KeyError:
               try:
                  da = cut['c03238'] # stupid t2 problem
               except KeyError:
                  pdb.set_trace()


            da.name = k
            da.longitude.values = da.longitude.values-360
            #encoding = {var: comp for var in da.data_vars}
            encoding = {k: {'complevel': 5, 'zlib': True}}
            if not os.path.exists(outfolder):
                os.makedirs(outfolder)

            da.to_netcdf(outfolder + os.sep + outname , format='NETCDF4', encoding=encoding)
            da.close()

            print('Wrote '+ outfolder + os.sep + outname)



