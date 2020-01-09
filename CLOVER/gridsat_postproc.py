import pandas as pd
import numpy as np
import xarray as xr
import ipdb
import matplotlib.pyplot as plt
import cartopy
from utils import u_plot as up
import cartopy.crs as ccrs
import os
import matplotlib as mpl
import pickle as pkl
from utils import constants as cnst
from utils import u_arrays as ua
from scipy.ndimage.measurements import label


def month():
    y1 = 1982
    y2 =2017#2017
    years = list(range(1983,1985)) #+ list(range(2004,2014))

    msg_folder = cnst.GRIDSAT
    fname='aggs/gridsat_WA_-70_-50-5000km2_monthly.nc'


    da = None
    for y in years:
        y = str(y)
        da1 = xr.open_dataset(cnst.GRIDSAT+'gridsat_WA_-50_' + y + '.nc')
        print('Doing '+y)
        da1['tir'] = da1['tir'].where(da1['tir'] <= -70)

        da1 = da1.resample(time='m').mean('time')
        try:
            da = xr.concat([da, da1], 'time')
        except TypeError:
            da = da1.copy()


    enc = {'tir': {'complevel': 5,  'zlib': True}}
    da.to_netcdf(msg_folder + fname, encoding=enc)


def month_boxmean():

    years = list(range(1983,2018))

    msg_folder = cnst.GRIDSAT

    #ipdb.set_trace()
    #if not os.path.isfile(msg_folder + fname):
    da = None
    da_box = None
    hov_box = None
    for y in years:
        y = str(y)
        da1 = xr.open_dataset(cnst.GRIDSAT + 'gridsat_WA_-40_1000km2_15-21UTC' + y + '.nc')
        print('Doing ' + y)
        da1['tir'] = da1['tir'].where((da1['tir'] <= -40) & (da1['tir'] >= -108) )
        #da1['tir'].values[da1['tir'].values == 0] = np.nan

        da_res = da1.resample(time='m').mean('time')
        WA_box = [-13,13,4.5,8]
        SAW_box = [15,25,-26,-18] #[20,30,-30,-10]
        SAE_box = [25,33,-28,-10]
        box = SAW_box
        boxed = da1['tir'].sel(lat=slice(box[2],box[3]), lon=slice(box[0],box[1])).resample(time='m').mean()

        try:
            da = xr.concat([da, da_res], 'time')
        except TypeError:
            da = da_res.copy()

        try:
            da_box = xr.concat([da_box, boxed], 'time')
        except TypeError:
            da_box = boxed.copy()
        da_box.attrs['box'] = box
        da_box.to_netcdf(msg_folder + 'aggs/SAboxWest_meanT-40_1000km2.nc')
        #da.to_netcdf(msg_folder + 'aggs/SAb_meanT-40_1000km2.nc')



def month_mean_hov():

        years = list(range(1983,2018))

        msg_folder = cnst.GRIDSAT
        fname = 'aggs/gridsat_WA_-50_monthly_mean.nc'

        hov_box = None
        for y in years:
            y = str(y)
            da1 = xr.open_dataset(cnst.GRIDSAT + 'gridsat_WA_-50_' + y + '.nc')
            print('Doing ' + y)
            da1['tir'] = da1['tir'].where((da1['tir'] <= -60) & (da1['tir'] >= -108) )

            WA_box = [-12,12,4.5,20]
            SA_box = [25,33,-28,-10]
            input = WA_box
            hov_boxed = da1['tir'].sel(lat=slice(input[2],input[3]), lon=slice(input[0],input[1])).resample(time='m').mean(['lon','time'])

            out = xr.DataArray(hov_boxed.values, coords={'month':hov_boxed['time.month'].values, 'lat':hov_boxed.lat}, dims=['month', 'lat'])


            try:
                hov_box = xr.concat([hov_box, out], 'year')
            except TypeError:
                hov_box = out.copy()

        hov_box.year.values = hov_box.year.values+years[0]
        hov_box.to_netcdf(msg_folder + 'aggs/WAbox_meanT-60_hov_5000km2.nc')



def month_mean_climatology():

    years = list(range(1983,2018))

    msg_folder = cnst.GRIDSAT
    fname = 'aggs/gridsat_WA_-50_monthly_mean.nc'

    #if not os.path.isfile(msg_folder + fname):
    da = None
    da_box = None
    for y in years:
        y = str(y)
        da1 = xr.open_dataset(cnst.GRIDSAT + 'gridsat_WA_-50_' + y + '.nc')
        print('Doing ' + y)
        da1['tir'] = da1['tir'].where((da1['tir'] <= -50) & (da1['tir'] >= -108) )
        #da1['tir'].values[da1['tir'].values < -70] = 1


        da_res = da1.resample(time='m').mean('time')

        #boxed = da1['tir'].sel(lat=slice(4.5,8), lon=slice(-10,10)).resample(time='m').mean()
        boxed = da1['tir'].sel(lat=slice(5, 8), lon=slice(-12, 10)).resample(time='m').mean()

        try:
            da = xr.concat([da, da_res], 'time')
        except TypeError:
            da = da_res.copy()

        try:
            da_box = xr.concat([da_box, boxed], 'time')
        except TypeError:
            da_box = boxed.copy()

        # pkl.dump(np.array(boxed),
        #          open('/users/global/cornkle/data/CLOVER/saves/box_13W-13E-4-8N_meanT-50_from5000km2.p',
        #               'wb'))

        enc = {'tir': {'complevel': 5, 'zlib': True}}
        #da_box.to_netcdf(msg_folder + 'box_13W-13E-4-8N_meanT-50_from5000km2.nc')
        da_box.to_netcdf(msg_folder + 'box_12W-10E-5-8N_meanT-50_from5000km2.nc')
        da.to_netcdf(msg_folder + fname, encoding=enc)


def month_count():

    years = list(range(1983, 2018))

    msg_folder = cnst.GRIDSAT
    fname = 'aggs/gridsat_WA_-70_monthly_mean_5000km2.nc' #65_monthly_count_-40base_15-21UTC_1000km2.nc'

    if not os.path.isfile(msg_folder + fname):
        da = None
        for y in years:
            y = str(y)
            da1 = xr.open_dataset(cnst.GRIDSAT + 'gridsat_WA_' + y + '.nc')  # _-40_1000km2_15-21UTC
            print('Doing ' + y)
            da1['tir'].values = da1['tir'].values#/100  ONLY FOR NEWER FILES
            da1['tir'] = da1['tir'].where((da1['tir'] <= -70) & (da1['tir'] >= -108), other=0) #-65
            da1['tir'].values[da1['tir'].values < 0] = 1
            #ipdb.set_trace()
            da1 = da1.resample(time='m').mean('time')
            try:
                da = xr.concat([da, da1], 'time')
            except TypeError:
                da = da1.copy()

        enc = {'tir': {'complevel': 5, 'zlib': True}}
        da.to_netcdf(msg_folder + fname, encoding=enc)


def month_count_sum():

    years = list(range(1983, 2018))

    msg_folder = cnst.GRIDSAT

    for y in years:
        y = str(y)
        da1 = xr.open_dataset(cnst.GRIDSAT + 'gridsat_WA_-70_5000km2_15-21UTC' + y + '.nc')
        print('Doing ' + y)
        da1['tir'].values = da1['tir'].values/100
        da1['tir'] = da1['tir'].where((da1['tir'] <= -70) & (da1['tir'] >= -108))
        da1['tir'].values[da1['tir'].values<=-70] = 1

        da1 = da1.resample(time='d').max('time').resample(time='m').mean('day') # percentage of days per month

        fname = '/gridsat_WA_-70_5000km2_15-21UTC' + y + '_monthSum.nc'
        enc = {'tir': {'complevel': 5, 'zlib': True}}
        da1.to_netcdf(msg_folder + fname, encoding=enc)


def month_count_concat():
    msg_folder = cnst.GRIDSAT
    fname = 'aggs/gridsat_WA_-70_monthly_count_15-21UTC_5000km2.nc'
    da = xr.open_mfdataset(cnst.GRIDSAT + 'monthSum/gridsat_WA_-70_*_monthSum.nc')

    enc = {'tir': {'complevel': 5, 'zlib': True}}
    da.to_netcdf(msg_folder + fname, encoding=enc)



    # else:
    #     ds = xr.open_dataset(msg_folder + fname)
    #     da = ds['tir']
    # pdb.set_trace()
    # da.values[da.values == 0] = np.nan
    # da.sel(lat=slice(11, 20))
    # mean = da['t'].mean(dim=['lat', 'lon'])
    #
    # mean.plot()
    #
    # plt.savefig('/users/global/cornkle/figs/CLOVER/GRIDSAT_cold_clouds/tests/trend_mcs-50.png', dpi=300)


def storm_count(area=False):
    msg_folder = cnst.GRIDSAT
    fname = msg_folder + 'gridsat_WA_-40_1000km2_15-21UTC'

    def makedic():
        mdic = {}
        for m in range(1,13):
            mdic[m] = []
        return mdic

    dic75 = makedic()
    dic70 = makedic()
    dic60 = makedic()
    dic50 = makedic()
    dic40 = makedic()

    area75 = makedic()
    area70 = makedic()
    area60 = makedic()
    area50 = makedic()
    area40 = makedic()

    for y in range(1983,2018): #2018
        ds = xr.open_dataset(fname + str(y) + '.nc')
        for m in range(1,13):


            da = ds['tir'][(ds['time.month'] == m) & (ds['time.hour']==18)]#(ds['time.hour']>=15) & (ds['time.hour']<=21)]
            da.values = da.values / 100

            val = 0
            storm = 0
            ar = []
            pixel = 78  # 78 # 78 = 5000km2 # 15000 = 253
            for d in da:

                cut = d.sel(lat=slice(5.2, 8), lon=slice(-10, 12))  # 4.5,8.5
                labels, goodinds = ua.blob_define(cut.values, -40, minmax_area=[pixel, 25000],
                                                  max_area=None)  # 7.7x7.7km = 64km2 per pix in gridsat?
                if area:
                    for gi in goodinds:
                        ar.append(np.sum(labels == gi))

                storm += np.float(goodinds.size)
                val += 1
            dic40[m].append(storm/val)
            if area:
                try:
                    area40[m].append(np.percentile(ar, 90))
                except IndexError:
                    area40[m].append(np.nan)


            val = 0
            storm = 0
            ar = []
            pixel = 78 #78 # 78 = 5000km2 # 15000 = 253
            for d in da:

                cut = d.sel(lat=slice(5.2, 8), lon=slice(-10, 12)) # 4.5,8.5
                labels, goodinds = ua.blob_define(cut.values, -50, minmax_area=[pixel, 25000],
                                              max_area=None)  # 7.7x7.7km = 64km2 per pix in gridsat?
                if area:
                    for gi in goodinds:
                        ar.append(np.sum(labels==gi))

                storm += np.float(goodinds.size)
                val += 1
            dic50[m].append(storm/val)
            if area:
                try:
                    area50[m].append(np.percentile(ar, 90))
                except IndexError:
                    area50[m].append(np.nan)

            val = 0
            storm = 0
            ar = []
            for d in da:
                cut = d.sel(lat=slice(5.2, 8), lon=slice(-10, 12))
                labels, goodinds = ua.blob_define(cut.values, -60, minmax_area=[pixel, 25000],
                                              max_area=None)  # 7.7x7.7km = 64km2 per pix in gridsat?
                if area:
                    for gi in goodinds:
                        ar.append(np.sum(labels==gi))

                storm += np.float(goodinds.size)
                print(m, goodinds, val)
                val += 1
            dic60[m].append(storm/val)
            if area:
                try:
                    area60[m].append(np.percentile(ar, 90))
                except IndexError:
                    area60[m].append(np.nan)

            val = 0
            storm = 0
            ar = []
            for d in da:
                cut = d.sel(lat=slice(5.2, 8), lon=slice(-10, 12))
                labels, goodinds = ua.blob_define(cut.values, -70, minmax_area=[pixel, 25000],
                                              max_area=None)  # 7.7x7.7km = 64km2 per pix in gridsat?
                if area:
                    for gi in goodinds:
                        ar.append(np.sum(labels==gi))

                storm += np.float(goodinds.size)
                val += 1
            dic70[m].append(storm/val)
            if area:
                try:
                    area70[m].append(np.percentile(ar, 90))
                except IndexError:
                    area70[m].append(np.nan)

            val = 0
            storm = 0
            ar = []
            for d in da:
                cut = d.sel(lat=slice(5.2, 8), lon=slice(-10, 12))
                labels, goodinds = ua.blob_define(cut.values, -75, minmax_area=[pixel, 25000],
                                              max_area=None)  # 7.7x7.7km = 64km2 per pix in gridsat?
                if area:
                    for gi in goodinds:
                        ar.append(np.sum(labels==gi))

                storm += np.float(goodinds.size)
                val += 1
            dic75[m].append(storm/val)
            if area:
                try:
                    area75[m].append(np.percentile(ar, 90))
                except IndexError:
                    area75[m].append(np.nan)
    print(40, dic40[3])
    print(50, dic50[3])
    print(60, dic60[3])
    print(70, dic70[3])
    print(75, dic75[3])


    pkl.dump(dic40, open(cnst.network_data + 'data/CLOVER/saves/storm_count_10W-12E_5-8N_-40C_5000km2_18.p', #4f5-8f5N
                        'wb'))

    pkl.dump(dic50, open(cnst.network_data + 'data/CLOVER/saves/storm_count_10W-12E_5-8N_-50C_5000km2_18.p', #4f5-8f5N
                        'wb'))

    pkl.dump(dic60, open(cnst.network_data + 'data/CLOVER/saves/storm_count_10W-12E_5-8N_-60C_5000km2_18.p',
                         'wb'))

    pkl.dump(dic70, open(cnst.network_data + 'data/CLOVER/saves/storm_count_10W-12E_5-8N_-70C_5000km2_18.p',
                         'wb'))

    pkl.dump(dic75, open(cnst.network_data + 'data/CLOVER/saves/storm_count_10W-12E_5-8N_-75C_5000km2_18.p',
                        'wb'))



    if area:

        pkl.dump(area40, open(cnst.network_data + 'data/CLOVER/saves/storm_90centArea_12W-10E_5-8N_-40C_5000km2_1800.p', #4f5-8f5N
                            'wb'))

        pkl.dump(area50, open(cnst.network_data + 'data/CLOVER/saves/storm_90centArea_12W-10E_5-8N_-50C_5000km2_1800.p', #4f5-8f5N
                            'wb'))

        pkl.dump(area60, open(cnst.network_data + 'data/CLOVER/saves/storm_90centArea_12W-10E_5-8N_-60C_5000km2_1800.p',
                             'wb'))

        pkl.dump(area70, open(cnst.network_data + 'data/CLOVER/saves/storm_90centArea_12W-10E_5-8N_-70C_5000km2_1800.p',
                             'wb'))

        pkl.dump(area75, open(cnst.network_data + 'data/CLOVER/saves/storm_90centArea_12W-10E_5-8N_-75C_5000km2_1800.p',
                            'wb'))



def storm_count_hov():
    msg_folder = cnst.GRIDSAT
    fname = msg_folder + 'gridsat_WA_-40_1000km2_15-21UTC'

    def makedic():
        mdic = {}
        for m in range(1,13):
            mdic[m] = []
        return mdic

    dic75 = makedic()
    dic70 = makedic()
    dic60 = makedic()
    dic50 = makedic()
    dic40 = makedic()


    for y in range(1983,2018): #2018
        ds = xr.open_dataset(fname + str(y) + '.nc')
        for m in range(1,13):


            da = ds['tir'][(ds['time.month'] == m) & (ds['time.hour']==18)]#(ds['time.hour']>=15) & (ds['time.hour']<=21)]
            da.values = da.values / 100
            da = da.sel(lat=slice(5.2, 8), lon=slice(-10, 12))

            val = 0
            storm = np.array([0]*da.shape[1])
            ar = []
            pixel = 78  # 78 # 78 = 5000km2 # 15000 = 253
            for d in da:

                cut = d.sel(lat=slice(5.2, 8), lon=slice(-10, 12))
                labels, goodinds = ua.blob_define(cut.values, -40, minmax_area=[pixel, 25000],
                                                  max_area=None)  # 7.7x7.7km = 64km2 per pix in gridsat?
                unarr = []
                for ll in labels:
                    isun = np.unique(ll)
                    num = isun.size
                    unarr.append(num)


                unarr = np.array(unarr)
                # if np.sum(unarr-1) > 0:
                #     ipdb.set_trace()

                storm += unarr-1
                val += 1

                print('-40 storm', storm)

            dic40[m].append(storm/val)



            val = 0
            storm = np.array([0]*da.shape[1])
            ar = []
            pixel = 78 #78 # 78 = 5000km2 # 15000 = 253
            for d in da:

                cut = d.sel(lat=slice(5.2, 8), lon=slice(-10, 12)) # 4.5,8.5
                labels, goodinds = ua.blob_define(cut.values, -50, minmax_area=[pixel, 25000],
                                              max_area=None)  # 7.7x7.7km = 64km2 per pix in gridsat?

                unarr = []
                for ll in labels:
                    isun = np.unique(ll)
                    num = isun.size
                    unarr.append(num)


                unarr = np.array(unarr)

                storm += unarr-1
                val += 1
            dic50[m].append(storm/val)


            val = 0
            storm = np.array([0]*da.shape[1])
            ar = []
            for d in da:
                cut = d.sel(lat=slice(5.2, 8), lon=slice(-10, 12))
                labels, goodinds = ua.blob_define(cut.values, -60, minmax_area=[pixel, 25000],
                                              max_area=None)  # 7.7x7.7km = 64km2 per pix in gridsat?
                unarr = []
                for ll in labels:
                    isun = np.unique(ll)
                    num = isun.size
                    unarr.append(num)


                unarr = np.array(unarr)

                storm += unarr-1
                val += 1

            dic60[m].append(storm/val)


            val = 0
            storm = np.array([0]*da.shape[1])
            ar = []
            for d in da:
                cut = d.sel(lat=slice(5.2, 8), lon=slice(-10, 12))
                labels, goodinds = ua.blob_define(cut.values, -70, minmax_area=[pixel, 25000],
                                              max_area=None)  # 7.7x7.7km = 64km2 per pix in gridsat?
                unarr = []
                for ll in labels:
                    isun = np.unique(ll)
                    num = isun.size
                    unarr.append(num)


                unarr = np.array(unarr)

                storm += unarr-1
                val += 1

            dic70[m].append(storm/val)


            val = 0
            storm = np.array([0]*da.shape[1])
            ar = []
            for d in da:
                cut = d.sel(lat=slice(5.2, 8), lon=slice(-10, 12))
                labels, goodinds = ua.blob_define(cut.values, -75, minmax_area=[pixel, 25000],
                                              max_area=None)  # 7.7x7.7km = 64km2 per pix in gridsat?
                unarr = []
                for ll in labels:
                    isun = np.unique(ll)
                    num = isun.size
                    unarr.append(num)


                unarr = np.array(unarr)

                storm += unarr-1
                val += 1

            dic75[m].append(storm/val)

    print(40, dic40[3])
    print(50, dic50[3])
    print(60, dic60[3])
    print(70, dic70[3])
    print(75, dic75[3])


    pkl.dump(dic40, open(cnst.network_data + 'data/CLOVER/saves/storm_HOVcount_10W-12E_5-8N_-40C_5000km2_18.p', #4f5-8f5N
                        'wb'))

    pkl.dump(dic50, open(cnst.network_data + 'data/CLOVER/saves/storm_HOVcount_10W-12E_5-8N_-50C_5000km2_18.p', #4f5-8f5N
                        'wb'))

    pkl.dump(dic60, open(cnst.network_data + 'data/CLOVER/saves/storm_HOVcount_10W-12E_5-8N_-60C_5000km2_18.p',
                         'wb'))

    pkl.dump(dic70, open(cnst.network_data + 'data/CLOVER/saves/storm_HOVcount_10W-12E_5-8N_-70C_5000km2_18.p',
                         'wb'))

    pkl.dump(dic75, open(cnst.network_data + 'data/CLOVER/saves/storm_HOVcount_10W-12E_5-8N_-75C_5000km2_18.p',
                        'wb'))

