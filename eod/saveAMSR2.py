# -*- coding: utf-8 -*-


import multiprocessing
import glob
from eod import rewrite_data
from shared.utils import constants as cnst

def saveNetcdf(day=True):

    if day:
        daystring = '_A_'
        dstring = 'day'
    else:
        daystring = '_D_'
        dstring = 'night'


    sm_folder = cnst.network_data + 'data/OBS/AMSRE/aqua/raw_'+dstring+'_AMSR2'
    pool = multiprocessing.Pool(processes=7)

    files = glob.glob(sm_folder+'/LPRM-AMSR*'+daystring+'*.nc4')
    print('start loop')

    for f in files:
        rewrite_data.rewrite_AMSR2(f)

    #res = pool.map(rewrite_data.rewrite_AMSR2, files)
