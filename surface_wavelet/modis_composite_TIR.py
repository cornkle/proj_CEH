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
from utils import u_met, constants, u_parallelise, u_gis, u_arrays

import scipy.ndimage.interpolation as inter
import pickle as pkl

matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)


def loop():

    for l in np.arange(0,24):
        print('Doing '+str(l))
        composite(l)

def composite(h):


    file = constants.MCS_POINTS_DOM

    hour = h

    msg = xr.open_dataarray(file)
    msg = msg[(msg['time.hour'] == hour) & (msg['time.minute'] == 0) & (
        msg['time.year'] >= 2006) & (msg['time.year'] <= 2010) & (msg['time.month'] >= 6) ]

    msg = msg.sel(lat=slice(10,20), lon=slice(-10,10))


    dic = u_parallelise.run_arrays(8,file_loop,msg,['ano', 'regional', 'cnt'])

    for k in dic.keys():
       dic[k] = np.nansum(dic[k], axis=0)


    pkl.dump(dic, open("/users/global/cornkle/figs/LSTA-bullshit/scales/new/composite_topo"+str(hour).zfill(2)+".p", "wb"))
    extent = dic['ano'].shape[1]/2

    f = plt.figure(figsize=(14, 4))
    ax = f.add_subplot(131)

    plt.contourf(dic['regional'] / dic['cnt'], cmap='RdBu_r',  vmin=-5, vmax=5)
    plt.plot(extent, extent, 'bo')

    # ax.set_xticklabels(np.array((np.linspace(0, extent*2, 5) - extent) * 3, dtype=int))
    # ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')
    plt.title('Regional anomaly, Nb cores: ' + str(np.max(dic['cnt'])) + '| ' + str(hour).zfill(2) + '00UTC, Jun-Sep',
              fontsize=10)

    ax = f.add_subplot(132)

    plt.contourf(dic['cnt'] , cmap='viridis')
    plt.plot(extent, extent, 'bo')
    # ax.set_xticklabels(np.array((np.linspace(0, extent*2, 5) - extent) * 3, dtype=int))
    # ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')
    plt.title('Regional anomaly', fontsize=10)

    ax = f.add_subplot(133)

    plt.contourf((dic['ano'] / dic['cnt']), cmap='viridis_r',  vmin=-75, vmax=-60) #-(rkernel2_sum / rcnt_sum)
    plt.plot(extent, extent, 'bo')
    # ax.set_xticklabels(np.array((np.linspace(0, extent*2, 5) - extent) * 3, dtype=int))
    # ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')
    plt.title('Seasonal anomaly',
              fontsize=10)


    plt.tight_layout()

    plt.savefig('/users/global/cornkle/figs/LSTA-bullshit/scales/new/composites_lsta/'+'TIRtopo_'+str(hour).zfill(2)+'00UTC')
    #plt.close()



def cut_kernel(xpos, ypos, arr):

    dist = 60

    kernel = u_arrays.cut_kernel(arr,xpos, ypos,dist)

    kernel_reg = kernel - np.nanmean(kernel)

    cnt = np.zeros_like(kernel)
    cnt[np.isfinite(kernel)] = 1

    if kernel.shape != (121, 121):
        pdb.set_trace()

    return kernel, kernel_reg, cnt



def file_loop(fi):


    file2 = '/users/global/cornkle/MCSfiles/blob_map_MCSs_-50_JJAS_points_dominant_gt15k.nc'
    mcs = xr.open_dataarray(file2)
    mcs = mcs.sel(lat=slice(10,20), lon=slice(-10,10))
    try:
        mcs = mcs[mcs['time']==fi['time']].squeeze()
    except ValueError:
        return
    mcs.values[mcs.values>-50]=np.nan
    topo = xr.open_dataset(constants.MSG5KM_TOPO)
    topo = topo.sel(lat=slice(10, 20), lon=slice(-10, 10))
    ttopo = topo['h']
    grad = np.gradient(ttopo.values)
    gradsum = abs(grad[0]) + abs(grad[1])

    mcs.values[ttopo.values>=450] = np.nan
    mcs.values[gradsum>30] = np.nan
    print('Doing day: ', fi['time'])
    pos = np.where( (fi.values >= 5) & (fi.values < 65))#(fi.values >= 1) & (fi.values <= 20)) #<-50)#

    if np.sum(pos) == 0:
        print('No blobs found')
        return

    kernel_list = []
    kernelreg_list = []
    cnt_list = []

    for y, x in zip(pos[0], pos[1]):

        try:
            kernel, kernel_reg, cnt = cut_kernel(x, y, mcs)
        except TypeError:
            continue

        kernel_list.append(kernel)
        kernelreg_list.append(kernel_reg)
        cnt_list.append(cnt)

    if kernel_list == []:
        return None
    print(len(kernel_list))
    if len(kernel_list) == 1:
      return None
    else:
        kernel_sum = np.nansum(np.stack(kernel_list, axis=0), axis=0)
        kernel2_sum = np.nansum(np.stack(kernelreg_list, axis=0), axis=0)
        cnt_sum = np.nansum(np.stack(cnt_list, axis=0), axis=0)

    print('Returning')

    return (kernel_sum, kernel2_sum, cnt_sum)

if __name__ == "__main__":
    composite()
