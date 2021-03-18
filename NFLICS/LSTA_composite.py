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
import glob


matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)


def diurnal_loop():

    for h in [14,16,18,20,22,0,2,4,6]:  #range(0,24)

        if h <= 15:
            shift = (12 + h) * (-1)
        else:
            shift = 12 - h

        composite(h, shift)


def composite(h, eh):
    #pool = multiprocessing.Pool(processes=8)

    path = cnst.network_data + 'figs/NFLICS/LSTA_stats_study/'

    for y in range(2004,2016):
        files = glob.glob(cnst.elements_drive + '/Africa/WestAfrica/cores_bigDomain/*_'+str(y)+'_*.nc')

        hour = h

        msg = xr.open_mfdataset(files)

        msg = (msg['small_scale'])[(msg['time.hour'] == hour ) & (msg['time.minute'] == 0) & (msg['time.month'] >= 6) & (msg['time.month'] <= 9)  ]

        msg = msg.sel(lat=slice(6.5,8.5), lon=slice(-11,18)).load()

        msg.attrs['refhour'] = h
        msg.attrs['eh'] = eh
        #
        dic = u_parallelise.run_arrays(5,file_loop,msg,['ano', 'regional', 'cnt'])
        #
        # res = []
        # for m in msg[0:20]:
        #     out = file_loop(m)
        #     res.append(out)
        # ipdb.set_trace()
        # return

        for k in dic.keys():
           dic[k] = np.nansum(dic[k], axis=0)
        print('File written')
        pkl.dump(dic, open(path + "/composite_new_SC110km2_"+str(y)+'_h'+str(hour).zfill(2)+".p", "wb"))



def cut_kernel(xpos, ypos, arr, dist, probs=False):


    kernel = u_arrays.cut_kernel(arr,xpos, ypos,dist)

    if (np.sum(kernel!=0) < 0.001 * kernel.size):
        print('Not enough values, return!')
        return

    kernel3 = kernel - np.nanmean(kernel)

    cnt = np.zeros_like(kernel)
    cnt[kernel!=0] = 1 #kslots[np.isfinite(kernel)]   #1

    # slot_kernel[np.isnan(kernel)] = 0

    if kernel.shape != (2*dist+1, 2*dist+1):
        print('Kernel shape does not fit! Return')
        return

    if np.nansum(probs) > 0:
        prob = u_arrays.cut_kernel(probs,xpos, ypos,dist)
        cnt3 = np.zeros_like(kernel)
        cnt3[prob!=0] = 1

    else:
        prob = np.zeros_like(kernel)
        cnt3 = np.zeros_like(kernel)


    return kernel, kernel3, cnt, prob, cnt3




def get_previous_hours_msg(date, ehour, refhour):

    date = date.replace(hour=refhour)

    if ehour > 0:
        edate = date + pd.Timedelta(str(ehour) + ' hours')
    else:
        edate = date - pd.Timedelta(str(np.abs(ehour)) + ' hours')
        #edate = edate.replace(hour=ehour)

    t1 = edate - pd.Timedelta('1 hours')
    t2 = edate + pd.Timedelta('1 hours')

    file = glob.glob(cnst.elements_drive + '/Africa/WestAfrica/cores_bigDomain/coresPower_MSG_-40_9-130km_-50points_dominant_'+str(date.year)+'*.nc')
    msg = xr.open_mfdataset(file).load()
    try:
        msg = msg['tir'].sel(time=slice(t1.strftime("%Y-%m-%dT%H"), t2.strftime("%Y-%m-%dT%H")))
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
        lsta = xr.open_dataset(cnst.elements_drive + 'Africa/WestAfrica/NFLICS/LSTA_2004-2015/netcdf_onCores/HDF5_LSASAF_ANOM_MSG_LST_MSG-Disk_' + fdate + '1700.nc') #_NEW
    except OSError:
        print('LSTA file missing')
        return
    print('Doing '+ 'HDF5_LSASAF_' + fdate + '.nc')

    lsta_da = lsta['lsta'].sel(lat=slice(4,20), lon=slice(-16.5,20)).squeeze()  # should be LSTA


    #### remove mean from LSTA_NEW

    #slot_da = lsta['NbSlot'].squeeze()

    # lsta_da.values[slot_da.values>15] = np.nan # just to test cases with clouds

    if (np.sum(np.isfinite(lsta_da)) / lsta_da.size) < 0.05:
        print('Not enough valid')
        return

    pos = np.where((fi.values <= -5)) # (fi.values >= -110) &    #(fi.values >= 5) & (fi.values < 65), fi.values<-40

    if (np.sum(pos) == 0):  #| (len(pos[0]) < 3)
        print('No blobs found')
        return

    probs_msg = get_previous_hours_msg(date, fi.attrs['eh'], fi.attrs['refhour'])

    dist = 100

    kernel2_list = np.zeros((dist*2+1, dist*2+1))
    kernel3_list = np.zeros((dist*2+1, dist*2+1))
    cnt_list = np.zeros((dist*2+1, dist*2+1))

    for y, x in zip(pos[0], pos[1]):

        lat = fi['lat'][y]
        lon = fi['lon'][x]

        point = lsta_da.sel(lat=lat, lon=lon, method='nearest')
        plat = point['lat'].values
        plon = point['lon'].values

        xpos = np.where(lsta_da['lon'].values == plon)
        xpos = int(xpos[0])
        ypos = np.where(lsta_da['lat'].values == plat)
        ypos = int(ypos[0])
        #
        #ipdb.set_trace()

        try:
            kernel2, kernel3, cnt, msg_kernel, mcnt = cut_kernel(xpos, ypos, lsta_da.values, dist)
        except TypeError:
            print('Type error for kernel, continue')
            continue

        # if np.nansum(msg_kernel[:,dist::])>=2:   # filter out cases with MCSs at 12
        #     print('Meteosat MCS continue')
        #     continue

        # kernel2_list.append(kernel2)
        # kernel3_list.append(kernel3)
        # cnt_list.append(cnt)

        kernel2_list += kernel2
        kernel3_list += kernel3
        cnt_list += cnt

    if np.nansum(kernel2_list) == 0:
        print('Kernel list is empty, return')
        #ipdb.set_trace()
        return

    # if kernel2_list == []:
    #     return None
    #
    # if len(kernel2_list) == 1:
    #   return None
    # else:

        # kernel2_sum = np.nansum(np.stack(kernel2_list, axis=0), axis=0)
        # kernel3_sum = np.nansum(np.stack(kernel3_list, axis=0), axis=0)
        # cnt_sum = np.nansum(np.stack(cnt_list, axis=0), axis=0)

    print('Returning with kernel')

    del lsta
    del fi
    #ipdb.set_trace()
    return (kernel2_list, kernel3_list, cnt_list) #(kernel2_sum, kernel3_sum, cnt_sum)



def plot(h):
    hour=h

    dic = pkl.load(open(
            cnst.network_data + 'figs/NFLICS/LSTA_stats_study/' + "/composite_new_SC110km2_"+str(2004)+'_h'+str(h).zfill(2)+".p", "rb"))

    def coll(dic, h, year):
        print(h)
        core = pkl.load(open(
            cnst.network_data + 'figs/NFLICS/LSTA_stats_study/' + "/composite_new_SC110km2_"+str(year)+'_h'+str(h).zfill(2)+".p", "rb"))
        for id, k in enumerate(core.keys()):
            try:
                dic[k] = dic[k] + core[k]
            except KeyError:
                dic[k] = core[k]

    for y in range(2005, 2016):
        coll(dic, h, y)

    extent = dic['ano'].shape[1]/2-1

    f = plt.figure(figsize=(14, 7))
    ax = f.add_subplot(231)

    plt.contourf(dic['regional']/100 / dic['cnt'], cmap='RdBu_r',  levels=[-0.5,-0.4,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5], extend='both')
    plt.plot(extent, extent, 'bo')

    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 5) - extent) * 3, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')
    plt.title('Regional anomaly, Nb cores: ' + str(np.max(dic['cnt'])) + '| ' + str(hour).zfill(2) + '00UTC, Jun-Sep',
              fontsize=10)

    ax = f.add_subplot(232)

    plt.contourf((dic['regional']/100 / dic['cnt'])  , cmap='RdBu_r',  levels=[-0.5,-0.4,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5], extend='both')
    plt.plot(extent, extent, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 5) - extent) * 3, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')
    plt.title('Regional anomaly', fontsize=10)

    ax = f.add_subplot(233)

    plt.contourf((dic['ano']/100 / dic['cnt']), cmap='RdBu_r',  levels=[-0.5,-0.4,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5], extend='both') #-(rkernel2_sum / rcnt_sum)
    plt.plot(extent, extent, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 5) - extent) * 3, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')
    plt.title('Seasonal anomaly',
              fontsize=10)

    ax = f.add_subplot(234)

    plt.contourf((dic['ano'] / 100 / dic['cnt']) , cmap='RdBu_r',  levels=[-0.5,-0.4,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5], extend='both')
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
    # plt.savefig(cnst.network_data + "/figs/LSTA-bullshit/AGU/" + +str(hour).zfill(2)+'_allplots.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png')
    # plt.close()

def plot_gewex(h):
    hour=h
    path = cnst.network_data + "figs/LSTA/corrected_LSTA/new/LSTA_only_plots/LSTA_old_vs_new_MCSfilter"
    dic = pkl.load(open(path + "/composite_old_"+str(hour).zfill(2)+".p", "rb"))

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

    plt.savefig(path +'/lsta_old_' + str(hour).zfill(2)+'_single.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png')
    plt.close()


def plot_all():

    hours = [14,16,18,20,22,0,2,4,6]#[15,16,17,18,19,20,21,22,23,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14]

    for h in hours:
        plot_gewex(h)
