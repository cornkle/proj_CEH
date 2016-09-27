# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 10:15:40 2016

@author: cornkle
"""

import numpy as np
import xarray as xr
from wavelet import util
from utils import u_arrays as ua
import pickle as pkl
from scipy import ndimage
import matplotlib.pyplot as plt
from eod import tm_utils
from collections import defaultdict
import multiprocessing

pool = multiprocessing.Pool(processes=8)


def composite(scale_ind):
    files = ua.locate(".nc", '/users/global/cornkle/MCSfiles/')

    dic = xr.open_dataset(files[0])
    dummy = util.waveletTP(dic['t_lag0'].values, dic['p'].values, 5)

    arr = np.array(dummy['scales'], dtype=str)   # the scales are period/2

    p = dic['p'].values.copy()

    comp_collect = defaultdict(list)
    cnt = 0
    for f in files:
        print('Doing file: '+f)
        dic = xr.open_dataset(f)

        outt = np.array(dic['t_lag0'].values.copy())
        p = np.array(dic['p'].values.copy())

        outt[np.isnan(outt)] = 150
        outt[outt > -40] = 150
        grad = np.gradient(outt)
        outt[outt > -40] = -40
        o2 = outt.copy()
        nok = np.where(abs(grad[0]) > 80)
        d = 2
        i = nok[0]
        j = nok[1]

        for ii, jj in zip(i, j):
            kernel = o2[ii - d:ii + d + 1, jj - d:jj + d + 1]
            #  if not kernel.any():
            #   continue
            #   else:
            o2[ii - d:ii + d + 1, jj - d:jj + d + 1] = ndimage.gaussian_filter(kernel, 3, mode='nearest')

        wav = util.waveletTP(o2, outt, 5)
        wl100 = wav['t'][24, :, :]

        for nb in scale_ind:

            scale = float(arr[nb])

            wl = wav['t'][nb,:,:]
            maxoutt = (wl == ndimage.maximum_filter(wl, scale/5, mode='constant', cval=np.amin(wl) + 1)) # scale/5
            maxoutt = maxoutt.astype(int)
            if scale < 60:
                yy, xx = np.where((maxoutt == 1) & (outt < -55) & (wl100 > 1) & (wl > np.percentile(wl,95)))
            else:
                yy, xx = np.where((maxoutt == 1) & (outt < -55) &  (wl > np.percentile(wl, 95)))

            # f = plt.figure()
            # siz = 5
            # ax = f.add_subplot(1,2,1)
            # plt.imshow(wav['t'][nb,:,:])
            # plt.plot(xx, yy, 'yo', markersize=siz)
            # plt.plot(xx, yy, 'ro', markersize=siz)
            # ax.set_title(str(scale), fontsize=12)
            #
            #
            # ax = f.add_subplot(1,2,2)
            # plt.imshow(o2)
            # plt.plot(xx, yy, 'yo', markersize=siz)
            # plt.plot(xx, yy, 'ro', markersize=siz)
            # ax.set_title(str(scale), fontsize=12)
            # cnt += 1
            #
            # continue

            for y, x in zip(yy,xx):
               # r = np.round(scale/5/2, 0)
                r = 20
               # print(r)
                kernel = tm_utils.cut_kernel(p, x, y, r)

                if kernel.shape != (r*2+1,r*2+1):
                    continue
                #print(scale, kernel.shape)

                comp_collect[scale].append(kernel)


    print('Saved ', cnt)

    keys = comp_collect.keys()
    f = plt.figure()
    siz = 3
    for ind, k in enumerate(keys):
        num = len(keys)
        arr = np.array(comp_collect[k])
        bla = np.nanmean(arr, 0)
        isbigger = arr > 30
        blab = np.nanmean(isbigger, 0)*100

        ax = f.add_subplot(3, num, 1+ind)
        plt.imshow(bla, vmin=0, vmax=5, cmap='viridis')
        plt.title(str(k)+' km')
        plt.plot(20, 20, 'ro', markersize=siz)
        cbar = plt.colorbar()
        cbar.set_label('mm h-1')

        ax = f.add_subplot(3, num, 6 + ind)
        plt.imshow(bla, cmap='viridis')
        plt.title(str(k) + ' km')
        plt.plot(20, 20, 'ro', markersize=siz)
        cbar = plt.colorbar()
        cbar.set_label('mm h-1')

        ax = f.add_subplot(3, num, 11 + ind)
        plt.imshow(blab, vmin=0, vmax=2, cmap='viridis')
        plt.title(str(k) + ' km')
        plt.plot(20, 20, 'ro', markersize=siz)
        cbar = plt.colorbar()
        cbar.set_label('P(>30mm) %')




def file_loop(f):





if __name__ == "__main__":
    #composite()
    llists = composite([0, 12, 24, -1, -13])






