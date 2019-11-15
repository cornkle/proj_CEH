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


matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)
matplotlib.rcParams['hatch.linewidth'] = 0.1

def run_hours():

    l = [14,15,16,17,18,19,20,21,22,23,0,1,2,3,4,5,6,7]
    for ll in l:
        composite(ll)


def rewrite_list():
        dic = pkl.load(
            open('/home/ck/DIR/cornkle/figs/LSTA/corrected_LSTA/new/cores_-60_filteredPoints_gt25000km2_table.p', "rb"))
        new = dic.copy()
        for k in new.keys():
            new[k] = []

        for k in dic.keys():
            lists = dic[k]
            for l in lists:
                new[k].extend(l)

        pkl.dump(new, open('/home/ck/DIR/cornkle/figs/LSTA/corrected_LSTA/new/cores_-60_filteredPoints_gt25000km2_table_new.p', "wb"))


def composite():
    pool = multiprocessing.Pool(processes=5)

    file = cnst.MCS_POINTS_DOM #MCS_TMIN #
    path = cnst.network_data + 'figs/LSTA/corrected_LSTA/new/' #corrected_LSTA/wavelet/large_scale

    #hour = hour

    msg = xr.open_dataarray(file)
    msg = msg[ (msg['time.minute'] == 0) & (
        msg['time.year'] >= 2006) & (msg['time.year'] <= 2010) & (msg['time.month'] >= 6) & (msg['time.month'] <= 9) ]  #(msg['time.hour'] == hour) &

    msg = msg.sel(lat=slice(10.2,19.3), lon=slice(-9.8, 9.8))

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

    dic_names = [ 'date','x_3k', 'y_3k', 'x_5k', 'y_5k', 'rlist' , 'elist', 'mlist',  'area' , 'csize', 'hour']

    for id, l in enumerate(dic_names):

            dic[l] = np.squeeze(res[:,id,...])

    pkl.dump(dic, open(path+"/cores_-60_filteredPoints_gt25000km2_table.p", "wb"))  #"+str(hour)+"
    print('Save file written!')



def cut_kernel(xpos, ypos, dist, msg, era ,cmorph):

    msg_kernel = u_arrays.cut_kernel(msg,xpos, ypos,dist)

    if msg_kernel.shape != (dist*2+1, dist*2+1):
        print('Kernels shape wrong!')
        ipdb.set_trace()

    era_kernel = u_arrays.cut_kernel(era['tciw'].values,xpos, ypos,dist)

    cmorph_kernel = u_arrays.cut_kernel(cmorph,xpos, ypos,dist)

    return msg_kernel, era_kernel, cmorph_kernel



def get_era_mcs(date, ehour, refhour):


    date = date.replace(hour=refhour)
    cm = xr.Dataset()

    if ehour > 0:
        edate = date + pd.Timedelta(str(ehour) + ' hours')
    else:
        edate = date - pd.Timedelta(str(np.abs(ehour)) + ' hours')

    t1 = edate
    file = cnst.ERA5

    try:
        css = xr.open_dataset(file + 'hourly/surface/ERA5_'+str(edate.year)+'_' + str(edate.month).zfill(2) + '_srfc.nc')
        css = uda.flip_lat(css)
    except:
        return None


    srfc_clim = xr.open_dataset(file + 'monthly/synop_selfmade/CLIM_2006-2010_new/ERA5_2006-2010_CLIM_'+str(edate.month).zfill(2)+'-'+str(edate.hour).zfill(2)+'_srfc_rw.nc').load()
    ## latitude in surface is already flipped, not for pressure levels though... ?!


    css = css.sel(longitude=slice(-13, 13), latitude=slice(8, 21))
    srfc_clim = srfc_clim.sel(longitude=slice(-13, 13), latitude=slice(8, 21))

    css = css.sel(time=t1)
    srfc_clim = srfc_clim.squeeze()

    cm['tciw'] = css['tciw'].squeeze() - srfc_clim['tciw'].squeeze()

    del srfc_clim
    del css

    return cm



def get_msg_mcs(date, ehour, refhour):


    date = date.replace(hour=refhour)

    if ehour > 0:
        edate = date + pd.Timedelta(str(ehour) + ' hours')
    else:
        edate = date - pd.Timedelta(str(np.abs(ehour)) + ' hours')
        #edate = edate.replace(hour=ehour)

    t1 = edate - pd.Timedelta('1 hours')
    t2 = edate + pd.Timedelta('1 hours')

    file = cnst.MCS_15K# MCS_15K #_POINTS_DOM
    msg = xr.open_dataarray(file)
    try:
        msg = msg.sel(time=slice(t1.strftime("%Y-%m-%dT%H"), t2.strftime("%Y-%m-%dT%H")))
    except OverflowError:
        return None

    pos = np.where((msg.values <= -40) )

    out = np.zeros_like(msg)
    out[pos] = 1
    out = np.sum(out, axis=0)
    out[out>0]=1

    msg = msg.sum(axis=0)*0

    xout = msg.copy()
    del msg
    xout.name = 'probs'
    xout.values = out

    return xout


def get_CMORPH(date):


    tdic = {14 : ('32 hours', '4 hours'),
            15: ('33 hours', '5 hours'),
            16 : ('34 hours', '6 hours'),  # 6am prev - 10am storm day
            17: ('35 hours', '7 hours'),
            18 : ('36 hours', '8 hours'),
            19 : ('37 hours', '9 hours'),
            20: ('38 hours', '10 hours'),
            21: ('39 hours', '11 hours'),
            22: ('40 hours', '12 hours'),
            23: ('41 hours', '13 hours'),
            0: ('42 hours', '14 hours'),
            1: ('43 hours', '15 hours'),
            2: ('44 hours', '16 hours'),
            3: ('45 hours', '17 hours'),
            4: ('46 hours', '17 hours'),
            5: ('47 hours', '18 hours'),
            6: ('48 hours', '19 hours'),
            7: ('49 hours', '20 hours')}

    before = pd.Timedelta(tdic[date.hour][0])
    before2 = pd.Timedelta(tdic[date.hour][1])

    t1 = date - before
    t2 = date - before2


    file = cnst.CMORPH
    try:
        cmm = xr.open_dataarray(file + 'CMORPH_WA_' + str(date.year) + '.nc')
    except:
        return None
    cmm = cmm.sel( time=slice(t1, t2)).sum(dim='time')
    # cmm = cmm.sel(lat=slice(10.9, 19), lon=slice(-9.8, 9.8))

    cm = cmm
    pos = np.where(cm.values>=10)

    out = np.zeros_like(cm)
    out[pos] = 1

    xout = cm.copy()
    xout.name = 'probs'
    xout.values = out

    del cm
    return xout



def file_loop(fi):


    print('Doing day: ', fi.time)

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

    try:
        lsta = xr.open_dataset(cnst.LSTA_NEW + 'lsta_daily_' + fdate + '.nc')  # 3km grid
    except OSError:
        return None
    print('Doing ' + 'lsta_daily_' + fdate + '.nc')

    topo = xr.open_dataset(cnst.LSTA_TOPO)
    ttopo = topo['h']
    grad = np.gradient(ttopo.values)
    gradsum = abs(grad[0]) + abs(grad[1])

    lsta_da = lsta['LSTA'].squeeze()


    if (np.sum(np.isfinite(lsta_da)) / lsta_da.size) < 0.05:
        print('Not enough valid')
        return None

    lsta_da.values[np.isnan(lsta_da.values)] = 0

    lsta_da.values[ttopo.values >= 450] = np.nan
    lsta_da.values[gradsum > 30] = np.nan
    pos = np.where((fi.values >= 5) & (fi.values <= 65)) # (fi.values >= 5) & (fi.values < 65) #(fi.values >= 5) & (fi.values < 65)

    if (np.sum(pos) == 0):
        print('No blobs found')
        return None


    ###############################Blob loop
    ref = int(fi['time.hour'].values)
    #ipdb.set_trace()
    if ref >= 14:
        eh= 12-ref
    else:
        eh = 12-(ref+24)

    era_mcs = get_era_mcs(daybefore, eh, ref)
    print('Era5 collect')

    try:
        era_on_lsta = lsta_da.salem.transform(era_mcs)
    except RuntimeError:
        print('Era5 on LSTA interpolation problem')
        return None
    del era_mcs

    msg_mcs = get_msg_mcs(daybefore, eh, ref)

    msg_on_lsta = lsta_da.salem.transform(msg_mcs, interp='nearest')
    del msg_mcs

    cmorph = get_CMORPH(daybefore)
    try:
        cmorph_on_lsta = lsta_da.salem.transform(cmorph)
    except RuntimeError:
        return None
    del cmorph
    del lsta


    mcs_hour = xr.open_dataarray(cnst.MCS_HOUR_DAILY)  ###5km grid
    mcsimage = xr.open_dataarray(cnst.MCS_15K)
    mcsimage = mcsimage.sel(time=fi.time, lat=slice(10.2,19.3), lon=slice(-9.8, 9.8))

    #size filter
    labels, goodinds = ua.blob_define(mcsimage.values, -50, minmax_area=[1000,25000], max_area=None)

    mlist = []
    rlist = []
    elist = []
    area_list = []
    csize_list = []
    temperature_list = []
    date_list = []
    xpos_3k = []
    ypos_3k = []
    xpos_5k = []
    ypos_5k = []
    hourlist = []

    for y, x in zip(pos[0], pos[1]):

        #size filter
        if (labels[y,x] not in goodinds) | (labels[y,x] == 0):
            print('MCS too small!!')
            continue
        # temperature filter
        if (mcsimage.values[y,x] > -60):
            print('Core too warm!!')
            continue


        lat = fi['lat'][y]
        lon = fi['lon'][x]

        #overlap filter
        mhour = mcs_hour.sel(time=pd.datetime(int(fdate[0:4]), int(fdate[4:6]), int(fdate[6:8]), 14), lon=lon,
                             lat=lat).values

        if mhour <= 13:
            mhour += 24

        if mhour == 0:
            mhour += 24

        chour = fi['time.hour'].values


        if (chour >= 0) & (chour <= 13):
            chour += 24
        if (mhour < chour) | (np.isnan(mhour)):
            print('Core overlaps: earliest:', mhour, ' core: ', chour)
            continue

        point = lsta_da.sel(lat=lat, lon=lon, method='nearest')
        plat = point['lat'].values
        plon = point['lon'].values

        xpos = np.where(lsta_da['lon'].values == plon)
        xpos = int(xpos[0])
        ypos = np.where(lsta_da['lat'].values == plat)
        ypos = int(ypos[0])

        dist = 200
        try:
            msg_kernel, era_kernel, cmorph_kernel = cut_kernel(xpos, ypos, dist, msg_on_lsta, era_on_lsta, cmorph_on_lsta)  #3km res kernels
        except TypeError:
            print('Kernel error')
            continue

        mflag = 0
        eflag = 0
        rflag = 0

        if np.nansum(msg_kernel[:,dist+5::])>=2:   # filter out cases with msg MCSs at 12
            mflag = 1

        if np.nanmax(era_kernel[:,dist+5::])>=0.05:   # filter out cases with era MCSs at 12
            eflag = 1

        rainarea = cmorph_kernel[110:150,181:221]  # [-270-150kmSouth] x120km wide (x120km high)
        if np.nansum(rainarea)/rainarea.size >= 0.25 :   # filter out cases with rainfall to the south
            rflag = 1

        date_list.append(daybefore)
        mlist.append(mflag)
        elist.append(eflag)
        rlist.append(rflag)
        area_list.append(np.sum(labels==labels[y,x])*25)
        temperature_list.append(mcsimage.values[y,x])
        csize_list.append(fi.values[y,x])
        xpos_3k.append(xpos)
        ypos_3k.append(ypos)
        xpos_5k.append(x)
        ypos_5k.append(y)
        hourlist.append(ref)

    return (date_list, xpos_3k, ypos_3k, xpos_5k, ypos_5k,rlist,elist,mlist,area_list,csize_list, hourlist)
# rcns_sum, rcwe_sum, cns_sum, cwe_sum,


if __name__ == "__main__":
    run_hours()
