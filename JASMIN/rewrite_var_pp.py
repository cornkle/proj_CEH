import numpy as np
import os
import glob
import ipdb
import iris

years = np.arange(1997,2007)
months = np.arange(1,13)
days = np.arange(1,31)

data = 'P25_hist/'
dir = '/media/ck/Elements/'#'/gws/nopw/j04/impala/users/cklein/CLOVER/'
ppfiles = dir + 'pp_01/'#'pp_raw/'+data
ncfiles = dir #+ 'pp_to_nc/'+data

var = 'a09103'
stream = 'ay488'

for y in years:
    for m in months:

        infiles = ppfiles + stream + 'a.pa' + str(y) + str(m).zfill(2) + '*.pp'
        outfiles = ncfiles + var +'_A1hr_mean_'+stream+'_25km_'+str(y)+str(m).zfill(2)+'010030-'+str(y)+str(m).zfill(2)+'302330.nc'

        print('Reading ', glob.glob(infiles))

        print('Doing ', outfiles)

        if os.path.isfile(outfiles):
            print(outfiles, ' exists, continue!')
            continue

        files = glob.glob(infiles)

        field = iris.load_cube(files)
        field.var_name = var

        # if os.path.isfile(outfiles):
        #     os.remove(outfiles)

        iris.save(field, '/media/ck/Elements/iris_test.nc', complevel=4, zlib=False)
        print('Written: ', var, outfiles)