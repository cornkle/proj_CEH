# -*- coding: utf-8 -*-


import multiprocessing
import glob
from eod import rewrite_data
import ipdb
import os
import xarray as xr
import numpy as np
import pandas as pd
from utils import constants as cnst, u_met
from scipy.interpolate import griddata
import salem


def saveMonthlyClim():
    tag = 'hist'
    data_path = '/media/ck/Elements/Africa/WestAfrica/CP4/CP4'+tag+'/'
    vars = [ 'tcwv', 't_pl', 'u_pl', 'v_pl', 'q_pl','lw_out_PBLtop', 'q_pl', 'lsRain'] #

    out = cnst.network_data + 'data/CP4/CLOVER/CP4_monthly_'+tag+'/'

    for v in vars:
        mf = xr.open_mfdataset(data_path+v+os.sep+"*.nc", combine='nested', concat_dim='time') #

        clim = mf.groupby('time.month').mean('time')

        for m in clim.month.values:
            #ipdb.set_trace()
            month = clim.sel(month=m).load()

            comp = dict(zlib=True, complevel=5)
            encoding = {var: comp for var in month.data_vars}

            outname = out + v + os.sep
            if not os.path.isdir(outname):
                os.mkdir(outname)
            outfile = outname + v + '_monthlyClim_'+tag+'_4km_'+str(m).zfill(2)+'.nc'
            month.to_netcdf(outfile, mode='w', encoding=encoding, format='NETCDF4')


def saveMonthly():
    tag = 'fut'
    data_path = '/media/ck/Elements/Africa/WestAfrica/CP4/CP4'+tag+'/'
    vars = [ 'tcwv', 'q_pl' ] # 'tcwv','t_pl', 'u_pl', 'v_pl','lw_out_PBLtop','lsRain'

    out = cnst.network_data + 'data/CP4/CLOVER/CP4_monthly_'+tag+'/'

    ff = glob.glob(data_path + 'lsRain' + os.sep + "*_" + str(2000) + str(8).zfill(2) + "*.nc")

    dummy = xr.open_dataset(ff[0])

    era_srfc = xr.open_dataset(
        '/home/ck/DIR/mymachine/ERA5/monthly/synoptic/srfc_1979-2019_monthly_synop_07x07.nc')

    dummy_on_ERA, lut = era_srfc.salem.lookup_transform(dummy['lsRain'], return_lut=True)

    for v in vars:
        for y in range(1998,2007):
            for m in range(1,13):


                try:
                    mf = xr.open_mfdataset(data_path+v+os.sep+"*_"+str(y)+str(m).zfill(2)+"*.nc", combine='by_coords', concat_dim='time') #
                except OSError:
                    continue
                #ipdb.set_trace()
                print('Doing ',v,' and did mfaggregation')
                try:
                    month = mf.sel(time=mf['time.hour']==12, pressure=[925,850,650]).load()
                except:
                    month = mf.sel(time=mf['time.hour'] == 12).load()


                month = month.mean('time')

                month[v].values[month[v].values<-10000] = np.nan

                month_on_e = era_srfc.salem.lookup_transform(month, lut=lut)

                comp = dict(zlib=True, complevel=5)
                encoding = {var: comp for var in month_on_e.data_vars}

                outname = out + v + os.sep
                if not os.path.isdir(outname):
                    os.mkdir(outname)
                outfile = outname + v + '_monthly_12UTC_0.7deg_'+tag+'_4km_'+str(y)+str(m).zfill(2)+'.nc'
                print('Saving', outfile)
                month_on_e.to_netcdf(outfile, mode='w', encoding=encoding, format='NETCDF4')

                del mf
                del month



