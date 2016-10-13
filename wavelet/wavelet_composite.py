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
    files = ua.locate(".nc", '/users/global/cornkle/MCSfiles/WA30/')   # /WA30/
    print('Nb files', len(files))
    tt = 'WA30'

    comp_collect = {}

    xmatrix = np.array([(np.arange(0, 41) - 20) * 5, ] * 41)
    ymatrix = np.array([(np.arange(0, 41) - 20) * 5, ] * 41).transpose()

    dist = ua.distance(xmatrix, ymatrix,0,0)
    dist = np.round(dist)
    dist = dist.astype(int)

    res = pool.map(file_loop, files)
    pool.close()

    res = [item for sublist in res for item in sublist] # flatten list of lists

    #return res

    for v in res:
        comp_collect[v[0]]={'big': [], 'fin' : [], 'mean' : [], 'shape' : []}  #  big , fin, mean ,shape

    keys = comp_collect.keys()

    for v in res:
        comp_collect[v[0]]['big'].append(v[1])
        comp_collect[v[0]]['fin'].append(v[2])
        comp_collect[v[0]]['mean'].append(v[3])
        comp_collect[v[0]]['shape'].append(v[4])

    for k in keys:
        a = np.asarray(comp_collect[k]['big'])
        comp_collect[k]['big'] = a
        b = np.asarray(comp_collect[k]['fin'])
        comp_collect[k]['fin'] = b
        c = np.asarray(comp_collect[k]['mean'])
        comp_collect[k]['mean'] = c
        d = np.asarray(comp_collect[k]['shape'])
        comp_collect[k]['shape'] = d


    xp = (np.arange(0, 41)-20)*5
    f = plt.figure()
    siz = 3

    keys = comp_collect.keys()
    print(keys)

    #return comp_collect

    ######### 2d plots
    ll = [15,30, 60, 90, 202]#keys #[18,30,50,101,202]
    for ind, k in enumerate(ll):
        num = len(ll)
        arr = comp_collect[k]['mean']
        fin = comp_collect[k]['fin']
        big = comp_collect[k]['big']

        bla = np.nansum(arr, 0) / np.nansum(fin, 0)
        blab = np.nansum(big, 0) / np.nansum(fin, 0) *100

        ax = f.add_subplot(3, num, 1 + ind)
        plt.imshow(bla, vmin=0, vmax=3, cmap='viridis')
        plt.title(str(k) + ' km', fontsize=9)
        plt.plot(20, 20, 'ro', markersize=siz)
        cbar = plt.colorbar()
        cbar.set_label('mm h-1')

        ax = f.add_subplot(3, num, 6 + ind)
        plt.imshow(bla, cmap='viridis')
        plt.title(str(k) + ' km', fontsize=9)
        plt.plot(20, 20, 'ro', markersize=siz)
        cbar = plt.colorbar()
        cbar.set_label('mm h-1')

        ax = f.add_subplot(3, num, 11 + ind)
        plt.imshow(blab, vmin=0, vmax=3, cmap='viridis')
        plt.title(str(k) + ' km', fontsize=9)
        plt.plot(20, 20, 'ro', markersize=siz)
        cbar = plt.colorbar()
        cbar.set_label('P(>30mm) %')


    col = ['r', 'b', 'g', 'y', 'black']

    pplot=[]
    ll = keys
    for ind, k in enumerate(ll):
        num = len(ll)
        arr = comp_collect[k]['mean']
        fin = comp_collect[k]['fin']
        big = comp_collect[k]['big']
        shape = comp_collect[k]['shape']
        bbig = np.nansum(big, 0)
        ffin = np.nansum(fin, 0)
        sshape = np.nansum(shape)

        siz = []
        sum = []
        abssum = []
        for d in range(dist.min(), dist.max()+1):

            pos = np.where(dist == d)
            if not pos[0].any():
                continue

            ssum = np.sum(bbig[pos])
            fsum = np.sum(ffin[pos])
            siz.append(d)
            sum.append(ssum/fsum*100) # probability
            abssum.append(ssum)

        dmax = siz[np.argmax(sum)]  # Distance to point where we find max probability
        ddmax = abssum[np.argmax(sum)] # Number of pcp>30 where we find max probability
        pplot.append((k, dmax, dmax-np.max(sum), sshape, np.max(sum), ddmax, np.max(abssum)))

  #  ipdb.set_trace()

    f = plt.figure()
    for p in pplot:
        plt.scatter(p[0], p[1])
        plt.xlabel('T scales')
        plt.ylabel('Distance Pmax(Pcp>30mm)')
        plt.title(tt)

    f = plt.figure()
    for p in pplot:
        plt.scatter(p[0], p[4])
        plt.xlabel('T scales')
        plt.ylabel('Pmax(Pcp>30mm)')
        plt.title(tt)

    f = plt.figure()
    for p in pplot:
        plt.scatter(p[0], p[5])
        plt.xlabel('T scales')
        plt.ylabel('nb Pcp>30mm(Pmax)')
        plt.title(tt)

    f = plt.figure()
    for p in pplot:
        plt.scatter(p[0], p[3])
        plt.xlabel('T scales')
        plt.ylabel('nb in composite')
        plt.title(tt)



def file_loop(f):

    ret = []


    print('Doing file: ' + f)
    dic = xr.open_dataset(f)

    outt = dic['t_lag0'].values
    p = dic['p'].values

    outt[np.isnan(outt)] = 150
    outt[outt >= -40] = 150

    outt[outt >= -40] = -50
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
    wl100 = wav['t'][24, :, :]  # 24 is 60km

    arr = np.array(wav['scales'], dtype=str)

    scale_ind = range(arr.size)#[0, 12, 21, -1, -13]

    for nb in scale_ind:
        ar = []
        orig = float(arr[nb])
        scale = int(np.round(orig))

        wl = wav['t'][nb, :, :]

        maxoutt = (wl == ndimage.maximum_filter(wl,np.round(orig / 5), mode='constant', cval=np.amax(wl) + 1))  # scale/5

        yp, xp = np.where(p > 20)

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
            kernel = tm_utils.cut_kernel(p, x, y, r)

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







