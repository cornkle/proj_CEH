# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 10:15:40 2016

@author: cornkle
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib
import pdb
import pandas as pd
from utils import u_met, u_parallelise, u_gis, u_arrays, constants as cnst
import ipdb
import pickle as pkl


matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)


def diurnal_loop():

    for l in range(15,16): #[1,2,3,4,5]
        print('Doing '+str(l))
        composite(l)

def monthly_loop():

    for mm in range(6,10):
        for h,shift in zip(range(0,12), (12+np.arange(0,12))*(-1)):
            composite(h,shift,mm)



def composite(h, eh,mm):
    #pool = multiprocessing.Pool(processes=8)

    path = cnst.network_data + "figs/LSTA/corrected_LSTA/new/"
    file = cnst.MCS_POINTS_DOM

    hour = h

    msg = xr.open_dataarray(file)
    msg = msg[(msg['time.hour'] == hour ) & (msg['time.minute'] == 0) & (
        msg['time.year'] >= 2006) & (msg['time.year'] <= 2010) & (msg['time.month'] == mm) ]

    msg = msg.sel(lat=slice(10,20), lon=slice(-10,10))

    msg.attrs['refhour'] = h
    msg.attrs['eh'] = eh

    dic = u_parallelise.run_arrays(5,file_loop,msg,['ano', 'regional', 'cnt'])#, 'rano', 'rregional', 'rcnt'])
    #
    # res = []
    # for m in msg:
    #     out = file_loop(m)
    #     res.append(out)
    # return

    for k in dic.keys():
       dic[k] = np.nansum(dic[k], axis=0)
    print('File written')
    pkl.dump(dic, open(path + "/LSTA_permonth/composite_noMCSfilter_"+str(mm).zfill(2)+'_'+str(hour).zfill(2)+".p", "wb"))



def cut_kernel(xpos, ypos, arr, date, lon, lat, t, dist, parallax=False, rotate=False, probs=False):

    if parallax:
        km, coords = u_gis.call_parallax_era(date.month, t, lon, lat, 0, 0)
        lx, ly = km

        lx = int(np.round(lx / 3.))
        ly = int(np.round(ly / 3.))  # km into pixels
        xpos = xpos - lx
        ypos = ypos - ly

    kernel = u_arrays.cut_kernel(arr,xpos, ypos,dist)

    if rotate:
        kernel = u_met.era_wind_rotate(kernel,date,lat,lon,level=700, ref_angle=90)
    # plt.figure()
    # plt.imshow(kernel, origin='lower')

    #plt.pause(100000)
    if (np.sum(np.isfinite(kernel)) < 0.01 * kernel.size):
        return

    kernel3 = kernel - np.nanmean(kernel)

    cnt = np.zeros_like(kernel)
    cnt[np.isfinite(kernel)] = 1

    # slot_kernel[np.isnan(kernel)] = 0

    if kernel.shape != (201, 201):
        pdb.set_trace()


    if np.nansum(probs) > 0:
        prob = u_arrays.cut_kernel(probs,xpos, ypos,dist)
        cnt3 = np.zeros_like(kernel)
        cnt3[np.isfinite(prob)] = 1

    else:
        prob = np.zeros_like(kernel)
        cnt3 = np.zeros_like(kernel)


    return kernel, kernel3, cnt, prob, cnt3




def get_previous_hours_msg(date, ehour, refhour):

    # tdic = {18 : ('36 hours', '15 hours'),
    #         19 : ('37 hours', '16 hours'),
    #         20: ('38 hours', '17 hours'),
    #         21: ('39 hours', '18 hours'),
    #         22: ('40 hours', '19 hours'),
    #         23: ('41 hours', '20 hours'),
    #         0: ('42 hours', '21 hours'),
    #         3: ('45 hours', '24 hours'),
    #         6: ('48 hours', '27 hours')}
    # before = pd.Timedelta(tdic[date.hour][0])
    # before2 = pd.Timedelta(tdic[date.hour][1])
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

    #print(prev_time.strftime("%Y-%m-%dT%H"), date.strftime("%Y-%m-%dT%H"))
    pos = np.where((msg.values <= -40) ) #(msg.values >= 5) & (msg.values < 65)) # #

    out = np.zeros_like(msg)
    out[pos] = 1
    out = np.sum(out, axis=0)
    out[out>0]=1
    # if np.sum(out>1) != 0:
    #     'Stop!!!'
    #     pdb.set_trace()

    msg = msg.sum(axis=0)*0

    xout = msg.copy()
    del msg
    xout.name = 'probs'
    xout.values = out

    return xout


def file_loop(fi):

    timestr = str(fi['time.year'].values) + str(fi['time.month'].values).zfill(2) + str(fi['time.day'].values).zfill(2)
    print('Doing day: ', timestr+'_'+str(fi['time.hour'].values))
    date = pd.to_datetime(timestr)

    dayd = pd.Timedelta('1 days')
    if fi['time.hour'].values.size != 1:
        'hour array too big, problem!'

    if (fi['time.hour'].values) <= 16:
        print('Nighttime')
        daybefore = date - dayd
    else:
        print('Daytime')
        daybefore = date

    fdate = str(daybefore.year) + str(daybefore.month).zfill(2) + str(daybefore.day).zfill(2)

    try:
        lsta = xr.open_dataset(cnst.LSTA_NEW + 'lsta_daily_' + fdate + '.nc')
    except OSError:
        print('LSTA file missing')
        return None
    print('Doing '+ 'lsta_daily_' + fdate + '.nc')

    topo = xr.open_dataset(cnst.LSTA_TOPO)
    ttopo = topo['h']
    grad = np.gradient(ttopo.values)
    gradsum = abs(grad[0])+abs(grad[1])


    lsta_da = lsta['LSTA'].squeeze()  # should be LSTA
    #### remove mean from LSTA_NEW

    slot_da = lsta['NbSlot'].squeeze()
    # lsta_da.values[slot_da.values>15] = np.nan # just to test cases with clouds
    if (np.sum(np.isfinite(lsta_da)) / lsta_da.size) < 0.05:
        print('Not enough valid')
        return None


    lsta_da.values[ttopo.values>=450] = np.nan
    lsta_da.values[gradsum>30] = np.nan
    pos = np.where((fi.values >= 5) & (fi.values <= 50) )  #(fi.values >= 5) & (fi.values < 65), fi.values<-40

    if (np.sum(pos) == 0):  #| (len(pos[0]) < 3)
        print('No blobs found')
        return None

    probs_msg = get_previous_hours_msg(date, fi.attrs['eh'], fi.attrs['refhour'])

    msg_on_lsta = lsta.salem.transform(probs_msg, interp='nearest')
    del probs_msg


    kernel2_list = []
    kernel3_list = []
    cnt_list = []

    dist = 100

    # xfi = fi.shape[1]
    # yfi = fi.shape[0]
    # randx = np.random.randint(0,xfi,100)
    #
    # randy = np.random.randint(0,yfi,100)
    # posr = (randy, randx)

    # rkernel2_list = []
    # rkernel3_list = []
    # rcnt_list = []
    #
    # for y, x in zip(posr[0], posr[1]):
    #
    #
    #     lat = fi['lat'][y]
    #     lon = fi['lon'][x]
    #
    #     point = lsta_da.sel(lat=lat, lon=lon, method='nearest')
    #     plat = point['lat'].values
    #     plon = point['lon'].values
    #
    #     xpos = np.where(lsta_da['lon'].values == plon)
    #     xpos = int(xpos[0])
    #     ypos = np.where(lsta_da['lat'].values == plat)
    #     ypos = int(ypos[0])
    #
    #     try:
    #         rkernel2, rkernel3, rcnt= cut_kernel(xpos, ypos, lsta_da, daybefore, plon, plat, -40, parallax=False, rotate=False)
    #     except TypeError:
    #         continue
    #
    #     rkernel2_list.append(rkernel2)
    #     rkernel3_list.append(rkernel3)
    #     rcnt_list.append(rcnt)
    mcs_hour = xr.open_dataarray(cnst.MCS_HOUR_DAILY)
    for y, x in zip(pos[0], pos[1]):

        lat = fi['lat'][y]
        lon = fi['lon'][x]

        # mhour = mcs_hour.sel(time=pd.datetime(int(fdate[0:4]),int(fdate[4:6]), int(fdate[6:8]), 16), lon=lon, lat=lat).values
        # #ipdb.set_trace()
        #
        # if mhour < 16:
        #     mhour +=24
        #
        # if mhour == 0:
        #     mhour +=24
        #
        # chour = fi['time.hour'].values
        #
        # #ipdb.set_trace()
        #
        # if (chour >= 0) & (chour<=15):
        #     chour +=24
        #
        # if (mhour < chour) | (np.isnan(mhour)):
        #     print('Core overlaps: earliest:', mhour, ' core: ',  chour)
        #     continue

        t = fi.sel(lat=lat, lon=lon)

        point = lsta_da.sel(lat=lat, lon=lon, method='nearest')
        plat = point['lat'].values
        plon = point['lon'].values

        xpos = np.where(lsta_da['lon'].values == plon)
        xpos = int(xpos[0])
        ypos = np.where(lsta_da['lat'].values == plat)
        ypos = int(ypos[0])
        try:
            kernel2, kernel3, cnt, msg_kernel, mcnt = cut_kernel(xpos, ypos, lsta_da, daybefore.month, plon, plat, -40, dist,probs=msg_on_lsta, parallax=False, rotate=False)
        except TypeError:
            continue

        if np.nansum(msg_kernel[:,dist::])>=2:   # filter out cases with MCSs at 12
            print('Meteosat MCS continue')
            continue

        kernel2_list.append(kernel2)
        kernel3_list.append(kernel3)
        cnt_list.append(cnt)

    if kernel2_list == []:
        return None

    if len(kernel2_list) == 1:
      return None
    else:

        kernel2_sum = np.nansum(np.stack(kernel2_list, axis=0), axis=0)
        kernel3_sum = np.nansum(np.stack(kernel3_list, axis=0), axis=0)
        cnt_sum = np.nansum(np.stack(cnt_list, axis=0), axis=0)

        # rkernel2_sum = np.nansum((np.stack(rkernel2_list, axis=0)), axis=0)
        # rkernel3_sum = np.nansum((np.stack(rkernel3_list, axis=0)), axis=0)
        # rcnt_sum = np.nansum((np.stack(rcnt_list, axis=0)), axis=0)

    print('Returning')

    return (kernel2_sum, kernel3_sum, cnt_sum)#,  rkernel2_sum, rkernel3_sum, rcnt_sum)



def plot(h):
    hour=h
    path = cnst.network_data + "figs/LSTA/corrected_LSTA/new/LSTA_permonth/"
    dic = pkl.load(open(path + "/composite_"+str(hour).zfill(2)+".p", "rb"))

    extent = dic['ano'].shape[1]/2-1

    f = plt.figure(figsize=(14, 7))
    ax = f.add_subplot(231)

    plt.contourf(dic['regional'] / dic['cnt'], cmap='RdBu_r',  levels=[-0.5,-0.4,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5], extend='both')
    plt.plot(extent, extent, 'bo')

    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 5) - extent) * 3, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')
    plt.title('Regional anomaly, Nb cores: ' + str(np.max(dic['cnt'])) + '| ' + str(hour).zfill(2) + '00UTC, Jun-Sep',
              fontsize=10)

    ax = f.add_subplot(232)

    plt.contourf((dic['regional'] / dic['cnt'])  , cmap='RdBu_r',  levels=[-0.5,-0.4,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5], extend='both')
    plt.plot(extent, extent, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 5) - extent) * 3, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')
    plt.title('Regional anomaly', fontsize=10)

    ax = f.add_subplot(233)

    plt.contourf((dic['ano'] / dic['cnt']), cmap='RdBu_r',  levels=[-0.5,-0.4,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5], extend='both') #-(rkernel2_sum / rcnt_sum)
    plt.plot(extent, extent, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 5) - extent) * 3, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')
    plt.title('Seasonal anomaly',
              fontsize=10)

    ax = f.add_subplot(234)

    plt.contourf((dic['ano'] / dic['cnt']) , cmap='RdBu_r',  levels=[-0.5,-0.4,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5], extend='both')
    plt.plot(extent,extent, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 5) - extent) * 3, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')
    plt.title('Seasonal anomaly - random',
              fontsize=10)

    ax = f.add_subplot(235)

    plt.contourf(dic['cnt'], cmap='viridis') #-(rkernel2_sum / rcnt_sum)
    plt.plot(extent, extent, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 5) - extent) * 3, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')
    plt.title('Valid count',
              fontsize=10)

    ax = f.add_subplot(236)

    plt.contourf(dic['cnt'], cmap='viridis') #-(rkernel2_sum / rcnt_sum)
    plt.plot(extent,extent, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 5) - extent) * 3, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')
    plt.title('Random valid count',
              fontsize=10)

    plt.tight_layout()
    plt.savefig(cnst.network_data + "/figs/LSTA-bullshit/AGU/" + +str(hour).zfill(2)+'_allplots.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png')
    plt.close()

def plot_gewex(h,mm):
    hour=h
    path = cnst.network_data + "figs/LSTA/corrected_LSTA/new/LSTA_permonth/"+str(mm).zfill(2)+'_noMCSfilter'
    dic = pkl.load(open(path + "/composite_noMCSfilter_"+str(mm).zfill(2)+'_'+str(hour).zfill(2)+".p", "rb"))

    extent = (dic['ano'].shape[1]-1)/2

    f = plt.figure(figsize=(7, 5))
    ax = f.add_subplot(111)
    print(dic['ano'].shape)
    plt.contourf((dic['ano'] / dic['cnt']), cmap='RdBu_r', extend='both',levels=[ -0.5,-0.4,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5]) #-(rkernel2_sum / rcnt_sum)  levels=[ -0.5,-0.4,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5]
    plt.plot(extent, extent, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')

    plt.title(str(hour).zfill(2)+'00 UTC | '+str(np.max(dic['cnt']))+' cores', fontsize=17)


    plt.tight_layout()

    plt.savefig(path +'/lsta_' + str(hour).zfill(2)+'_single.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png')
    plt.close()


def plot_all():

    for m in range(6, 10):
        for h in range(15,24):
            plot_gewex(h,m)


def plot_all_multi():



    for m in range(6, 10):
        path = cnst.network_data + "figs/LSTA/corrected_LSTA/new/LSTA_permonth/" + str(m).zfill(2)
        f = plt.figure()

        for ids, h in enumerate([17,18,19,20,21,22,23,0,1,2,3,4]):

                ax = f.add_subplot(4,3,ids+1)

                hour = h

                dic = pkl.load(open(path + "/composite_" + str(m).zfill(2) + '_' + str(hour).zfill(2) + ".p", "rb"))

                extent = (dic['ano'].shape[1] - 1) / 2

                print(dic['ano'].shape)
                plt.contourf((dic['ano'] / dic['cnt']), cmap='RdBu_r', extend='both',
                             levels=[-0.5, -0.4, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4,
                                     0.5])  # -(rkernel2_sum / rcnt_sum)  levels=[ -0.5,-0.4,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5]
                plt.plot(extent, extent, 'bo')
                ax.set_xticklabels(np.array((np.linspace(0, extent * 2, 9) - extent) * 3, dtype=int))
                ax.set_yticklabels(np.array((np.linspace(0, extent * 2, 9) - extent) * 3, dtype=int))
                ax.set_xlabel('km')
                ax.set_ylabel('km')
                plt.colorbar(label='K')

                plt.title(str(hour).zfill(2) + '00 UTC | ' + str(np.max(dic['cnt'])) + ' cores', fontsize=17)

        plt.tight_layout()

        plt.savefig("figs/LSTA/corrected_LSTA/new/LSTA_permonth/multihour_LSTA_"+str(m).zfill(2)+'.png')
        plt.close()