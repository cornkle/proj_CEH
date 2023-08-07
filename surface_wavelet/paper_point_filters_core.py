# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 10:15:40 2016

@author: cornkle
"""

import numpy as np
import xarray as xr

import matplotlib
import multiprocessing
import ipdb
import pandas as pd
from utils import u_arrays, constants as cnst, u_met
import pickle as pkl
from utils import u_arrays as ua, u_darrays as uda
import salem
import matplotlib.pyplot as plt


matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)
matplotlib.rcParams['hatch.linewidth'] = 0.1

def run_hours():

    l = [15,16,17,18,19,20,21,22,23, 0,1,2,3,4,5,6,7]  #
    for ll in l:
        composite(ll)


key = '2hOverlap'
box = [9,19.5,-11.5,11.5]


def rewrite_list(hour):
        path = '/home/ck/DIR/cornkle/figs/LSTA/corrected_LSTA/new/wavelet_coefficients/core_txt/'
        dic = pkl.load(
            open(path + "cores_gt15000km2_table_"+str(hour)+"_" + key + ".p", "rb"))
        new = dic.copy()

        for k in new.keys():
            new[k] = []

        for k in dic.keys():
            lists = dic[k]
            for l in lists:
                new[k].extend(l)

        pkl.dump(new, open(path + "cores_gt15000km2_table_"+str(hour)+"_" + key + ".p", "wb"))


        df = pd.DataFrame.from_dict(new)
        df = df.reindex(columns=['year', 'month', 'day', 'hour', 'lon', 'lat', 'xloc', 'yloc', 'area', 'csize', 't', 'storm_id', 'topo', 'dtime'])
        df.to_csv(path + "cores_gt15000km2_table_1640_580_"+str(hour)+"_" + key + ".csv", na_rep=-999, index_label='id')


def composite(hour):
    pool = multiprocessing.Pool(processes=5)

    file = cnst.MCS_POINTS_DOM #MCS_TMIN #
    path = cnst.network_data + 'figs/LSTA/corrected_LSTA/new/wavelet_coefficients/' #corrected_LSTA/wavelet/large_scale

    msg = xr.open_dataarray(file)
    msg = msg[(msg['time.hour'] == hour) & (msg['time.minute'] == 0) & (
        msg['time.year'] >= 2006) & (msg['time.year'] <= 2010) & (msg['time.month'] >= 6) & (msg['time.month'] <= 9) ]

    msg = msg.sel(lat=slice(box[0],box[1]), lon=slice(box[2], box[3]))

    res = pool.map(file_loop, msg)
    pool.close()
    print('return parallel')
    # res = []
    # for m in msg[0:30]:
    #     out = file_loop(m)
    #     res.append(out)

    res = [x for x in res if x is not None]
    dic = {}

    res = np.array(res)

    dic_names = [ 'year', 'month', 'day', 'hour','lon', 'lat', 'area' , 'csize', 'xloc', 'yloc', 't','storm_id', 'topo', 'dtime']

    for id, l in enumerate(dic_names):

            dic[l] = np.squeeze(res[:,id,...])

    pkl.dump(dic, open(path+"core_txt/cores_gt15000km2_table_"+str(hour)+"_" + key + ".p", "wb"))  #"+str(hour)+"
    print('Save file written!')

    rewrite_list(hour)




def file_loop(fi):


    print('Doing day: ', fi.time)
    msg_latlon = np.load(cnst.network_data + 'data/OBS/MSG_WA30/MSG_1640_580_lat_lon.npz')
    mlon = msg_latlon['lon'].flatten()
    mlat = msg_latlon['lat'].flatten()

    msg_coords = np.array(list(zip(mlon,mlat)))

    date = pd.to_datetime(
        str(fi['time.year'].values) + str(fi['time.month'].values).zfill(2) + str(fi['time.day'].values).zfill(2))

    dayd = pd.Timedelta('1 days')
    if fi['time.hour'].values.size != 1:
        'hour array too big, problem!'

    if (fi['time.hour'].values) <= 13:
        print('Nighttime')
        daybefore = date - dayd
    else:
        print('Daytime')
        daybefore = date

    fdate = str(daybefore.year) + str(daybefore.month).zfill(2) + str(daybefore.day).zfill(2)

    outdate = pd.to_datetime(
        str(fi['time.year'].values) + str(fi['time.month'].values).zfill(2) + str(fi['time.day'].values).zfill(2))

    pos = np.where((fi.values >= 5) & (fi.values <= 65)) # (fi.values >= 5) & (fi.values < 65) #(fi.values >= 5) & (fi.values < 65)

    if (np.sum(pos[0]) == 0):
        print('No blobs found')
        return None


    ###############################Blob loop
    ref = int(fi['time.hour'].values)
    #ipdb.set_trace()
    if ref >= 14:
        eh= 12-ref
    else:
        eh = 12-(ref+24)


    mcs_hour = xr.open_dataarray(cnst.MCS_HOUR_DAILY)  ###5km grid
    mcsimage = xr.open_dataarray(cnst.MCS_15K)
    mcsimage = mcsimage.sel(time=fi.time, lat=slice(box[0],box[1]), lon=slice(box[2], box[3]))
    mcs_hour = mcs_hour.sel(lat=slice(box[0],box[1]), lon=slice(box[2], box[3]))

    topo = xr.open_dataset(cnst.WA_TOPO_3KM)
    ttopo = topo['h']

    #size filter
    labels, goodinds = ua.blob_define(mcsimage.values, -50, minmax_area=[600,100000], max_area=None)


    area_list = []
    csize_list = []
    temperature_list = []
    yearlist = []
    daylist = []
    monthlist = []
    hourlist = []
    xpos_3k = []
    ypos_3k = []
    storm_id = []
    topo = []
    core_size = []
    dtime = []


    xloc = []
    yloc = []

    for y, x in zip(pos[0], pos[1]):

        #size filter
        if (labels[y,x] not in goodinds) | (labels[y,x] == 0):
            print('MCS too small!!')
            continue
        # temperature filter
        if (mcsimage.values[y,x] > -60):
            print('Core too warm!!')
            continue

        si = labels[y,x]
        lat = fi['lat'][y]
        lon = fi['lon'][x]

        h = ttopo.sel(lat=lat, lon=lon, method='nearest', tolerance=0.04)

        #overlap filter
        try:
            mhour = mcs_hour.sel(time=pd.datetime(int(fdate[0:4]), int(fdate[4:6]), int(fdate[6:8]), 14), lon=lon,
                                lat=lat).values
        except:
            continue

        if mhour <= 13:
            mhour += 24

        chour = fi['time.hour'].values


        if (chour >= 0) & (chour <= 13):
            chour += 24

        deltatime = chour - mhour

        # if (mhour+2 < chour): #| (np.isnan(mhour)
        #     print('Core overlaps: earliest:', mhour, ' core: ', chour)
        #     continue

        if np.isnan(mhour):
            plt.figure()
            plt.contourf(mcsimage)
            plt.contour(mhour, levels=[14,chour])
            plt.plot(y,x,'ro')


        isnear = ua.closest_point(np.array([lon,lat]), msg_coords)
        ispoint = msg_coords[isnear]
        loc = np.unravel_index(isnear,msg_latlon['lon'].shape)


        area_list.append(np.sum(labels==labels[y,x])*25)
        temperature_list.append(mcsimage.values[y,x])
        csize_list.append(fi.values[y,x])
        xpos_3k.append(float(lon.values))
        ypos_3k.append(float(lat.values))

        hourlist.append(int(fi['time.hour'].values))
        monthlist.append(int(fi['time.month'].values))
        daylist.append(int(fi['time.day'].values))
        yearlist.append(int(fi['time.year'].values))
        xloc.append(loc[1])
        yloc.append(loc[0])
        storm_id.append(si)
        topo.append(h.values)
        dtime.append(deltatime)

        print('llstart', (float(lon.values), float(lat.values)), 'isMSG', ispoint)


    return (yearlist, monthlist, daylist , hourlist, xpos_3k, ypos_3k,area_list,csize_list, xloc,yloc,temperature_list, storm_id, topo, dtime)
# rcns_sum, rcwe_sum, cns_sum, cwe_sum,


if __name__ == "__main__":
    run_hours()
