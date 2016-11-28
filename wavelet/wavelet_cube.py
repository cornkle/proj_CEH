# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 10:15:40 2016

@author: cornkle
"""

import numpy as np
import xarray as xr
from wavelet import util
from utils import u_arrays as ua
from scipy import ndimage
import matplotlib.pyplot as plt
from eod import tm_utils
from collections import defaultdict
import multiprocessing
import ipdb


def composite():
    pool = multiprocessing.Pool(processes=6)
    files = ua.locate(".nc", '/users/global/cornkle/MCSfiles/WA15_big_-40_15W-20E/')   # /WA30/
    files = files[0:1000]
    print('Nb files', len(files))
    tt = 'WA15'

    res = pool.map(file_loop, files)
    pool.close()


def file_loop(f):

    ret = []

    print('Doing file: ' + f)
    dic = xr.open_dataset(f)

    outt = dic['tc_lag0'].values
    outp = dic['p'].values
    area = np.sum(outt <= -40)

    tt = np.min(outt[(np.isfinite(outp)) & (np.isfinite(outt))])
    pp = np.max(outp[(np.isfinite(outp)) & (np.isfinite(outt))])

    if (area*25 < 15000) or (pp<0.1) or (pp>200):
        return

    outt[np.isnan(outt)] = 150
    outt[outt >= -40] = 150
    grad = np.gradient(outt)
    o2 = outt.copy()
    nok = np.where(abs(grad[0]) > 80)

    outt[outt >= -40] = -50

    d = 2
    i = nok[0]
    j = nok[1]

    for ii, jj in zip(i, j):
        kern = outt[ii - d:ii + d + 1, jj - d:jj + d + 1]
        out[ii - d:ii + d + 1, jj - d:jj + d + 1] = ndimage.gaussian_filter(kern, 3, mode='nearest')

    wav = util.waveletT(o2, 5)
    wl100 = wav['t'][24, :, :]  # 24 is 60km

    arr = np.array(wav['scales'], dtype=str)

    scale_ind = range(arr.size)#[0, 12, 21, -1, -13]

    for nb in scale_ind:
        ar = []
        orig = float(arr[nb])
        scale = int(np.round(orig))

        wl = wav['t'][nb, :, :]

        maxoutt = (wl == ndimage.maximum_filter(wl,np.round(orig / 5), mode='constant', cval=np.amax(wl) + 1))  # scale/5

        yp, xp = np.where(outp > 20)

        # if scale < 60:
        #     yy, xx = np.where(
        #         (maxoutt == 1) & (outt < -60) & (wl > np.percentile(wl,90)))# & (wl100 > 10))  # & np.isfinite(p) ) # & (wl > np.percentile(wl,98))
        # else:
        yy, xx = np.where((maxoutt == 1) & (outt < -60) & (wl > np.percentile(wl,90)))# & (p>=0))  # & np.isfinite(p)) & (wl > np.percentile(wl,90))

        # if (scale==101) or (scale==60) or (scale==202) or (scale==19):
        #     f = plt.figure()
        #     siz = 5
        #     ax = f.add_subplot(1, 3, 1)
        #     plt.imshow(wav['t'][nb, :, :])
        #     plt.plot(xp, yp, 'yo', markersize=siz)
        #     plt.plot(xx, yy, 'ro', markersize=siz)
        #     ax.set_title(str(scale), fontsize=12)
        #
        #     ax = f.add_subplot(1, 3, 2)
        #     plt.imshow(o2)
        #     plt.plot(xp, yp, 'yo', markersize=siz)
        #     plt.plot(xx, yy, 'ro', markersize=siz)
        #     ax.set_title(str(scale), fontsize=12)
        #
        #     ax = f.add_subplot(1, 3, 3)
        #     plt.imshow(p)
        #     #plt.plot(xp, yp, 'yo', markersize=siz)
        #     plt.plot(xx, yy, 'ro', markersize=siz)
        #     ax.set_title(str(scale), fontsize=12)

        yyy=[]
        xxx=[]
        for y, x in zip(yy, xx):

            r = 20
            kernel = tm_utils.cut_kernel(outp, x, y, r)

            if kernel.shape != (r * 2 + 1, r * 2 + 1):
                continue

            if np.nansum(kernel)==0:
                continue
            xxx.append(x)
            yyy.append(y)
            ar.append(kernel)
        if ar==[]:
            continue

        # if (scale==57) or (scale==60) or (scale==18) or (scale==19):
        #     f = plt.figure()
        #     siz = 5
        #     ax = f.add_subplot(1, 3, 1)
        #     plt.imshow(wav['t'][nb, :, :])
        #     plt.plot(xp, yp, 'yo', markersize=siz)
        #     plt.plot(xxx, yyy, 'ro', markersize=siz)
        #     ax.set_title(str(scale), fontsize=12)
        #
        #     ax = f.add_subplot(1, 3, 2)
        #     plt.imshow(o2)
        #     plt.plot(xp, yp, 'yo', markersize=siz)
        #     plt.plot(xxx, yyy, 'ro', markersize=siz)
        #     ax.set_title(str(scale), fontsize=12)
        #
        #     ax = f.add_subplot(1, 3, 3)
        #     plt.imshow(p)
        #     #plt.plot(xp, yp, 'yo', markersize=siz)
        #     plt.plot(xxx, yyy, 'ro', markersize=siz)
        #     ax.set_title(str(scale), fontsize=12)

        ar = np.asarray(ar)

        isbigger = ar >= 30
        isfinite = np.isfinite(ar)
        isbig = np.nansum(isbigger, 0)
        isfin = np.nansum(isfinite, 0)
        mean = np.nansum(ar, 0)

        #print(scale, ar.shape)

        ret.append((scale, isbig, isfin, mean, ar.shape[0]))

        dic.close()

    return ret


def file_loop_p(f):

    scale_ind = [0,12,21,-1,-13]

    ret = []

    print('Doing file: ' + f)
    dic = xr.open_dataset(f)

    outt = dic['tc_lag0'].values
    p = dic['p'].values
    # p2 = p.copy()
    # p2[~np.isfinite(p2)]=-0.001

    outt[np.isnan(outt)] = 150
    outt[outt >= -40] = 150

    outt[outt >= -40] = -40
    grad = np.gradient(outt)
    o2 = outt.copy()
    nok = np.where(abs(grad[0]) > 80)
    d = 2
    i = nok[0]
    j = nok[1]

    for ii, jj in zip(i, j):
        kern = o2[ii - d:ii + d + 1, jj - d:jj + d + 1]
        o2[ii - d:ii + d + 1, jj - d:jj + d + 1] = ndimage.gaussian_filter(kern, 3, mode='nearest')

    wav = util.waveletT(o2, 5)
    wl100 = wav['t'][21, :, :]  # 24 is 60km

    arr = np.array(wav['scales'], dtype=str)

    for nb in scale_ind:
        orig = float(arr[nb])
        scale = np.round(orig)

        wl = wav['t'][nb, :, :]
        # pl = wav['p'][nb, :, :]
        maxoutt = (wl == ndimage.maximum_filter(wl,np.round(orig / 5), mode='constant', cval=np.amax(wl) + 1))  # scale/5

        #filter = ndimage.minimum_filter(outt, np.round(orig / 5), mode='constant', cval=np.min(outt) - 1)
        #filter[filter == np.max(outt)] = -999
        # filter = ndimage.maximum_filter(wl, np.round(orig / 5), mode='constant', cval=np.max(wl) + 1)
        # filter[filter == np.min(wl)] = -999
        # maxoutt = (wl == filter)

        # yy, xx = np.where((maxoutt == 1) & (outt <-55))

        yp, xp = np.where(p > 30)

        if scale < 31:
            yy, xx = np.where(
                (maxoutt == 1) & (outt < -60) & (wl100 > 20))  # & np.isfinite(p) ) # & (wl > np.percentile(wl,98))
        else:
            yy, xx = np.where(
                (maxoutt == 1) & (outt < -60)  )  # & np.isfinite(p)) & (wl > np.percentile(wl,90))

        # f = plt.figure()
        # siz = 5
        # ax = f.add_subplot(1, 3, 1)
        # plt.imshow(wav['t'][nb, :, :])
        # plt.plot(xp, yp, 'yo', markersize=siz)
        # plt.plot(xx, yy, 'ro', markersize=siz)
        # ax.set_title(str(scale), fontsize=12)
        #
        # ax = f.add_subplot(1, 3, 2)
        # plt.imshow(o2)
        # plt.plot(xp, yp, 'yo', markersize=siz)
        # plt.plot(xx, yy, 'ro', markersize=siz)
        # ax.set_title(str(scale), fontsize=12)
        #
        # ax = f.add_subplot(1, 3, 3)
        # plt.imshow(p)
        # #plt.plot(xp, yp, 'yo', markersize=siz)
        # plt.plot(xx, yy, 'ro', markersize=siz)
        # ax.set_title(str(scale), fontsize=12)

        for y, x in zip(yy, xx):
            # r = np.round(scale/5/2, 0)
            r = 20
            kernel = tm_utils.cut_kernel(p, x, y, r)

            if kernel.shape != (r * 2 + 1, r * 2 + 1):
                continue



            ret.append((scale,kernel))

    return ret




if __name__ == "__main__":
    composite()







