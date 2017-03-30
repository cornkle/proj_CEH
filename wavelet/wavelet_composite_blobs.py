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
import scipy.ndimage.filters as filters
import matplotlib.pyplot as plt
import matplotlib
from eod import tm_utils
import multiprocessing
import pdb
from collections import OrderedDict
from scipy.ndimage.measurements import label
import pandas as pd
import pickle as pkl
matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)


def composite():
    pool = multiprocessing.Pool(processes=7)
    files = ua.locate(".nc", '/users/global/cornkle/MCSfiles/WA15_big_-40_15W-20E_size/')   # /WA30/
    out = '/users/global/cornkle/C_paper/wavelet/saves/pandas/'
    #files = files[0:100]
    print('Nb files', len(files))
    tt = 'WA15'

    comp_collect = {}
    precip = {}

    res = pool.map(file_loop, files)
    pool.close()
    res = [x for x in res if x is not None]

    nb_sys = len(res)

    print('Number systems: ', nb_sys)

    res = [item for sublist in res for item in sublist]  # flatten list of lists

    for v in res:
        comp_collect[v[2]] = {'p': [], 't': [], 'scale': [], 'hour': [], 'id': []}
        precip[v[2]] = []

        # ret.append((kernel, kernelt, sc, id, dic['time.hour'].values.tolist(),
        #             clat, clon, lat_min, lat_max, lon_min, lon_max, area,
        #             bulk_pmax, bulk_pmean, bulk_tmean, bulk_tmean_p, bulk_tmin_p, bulk_g30,
        #             circle_Tcenter, circle_p, circle_t, circle_valid, circle_sum,
        #             circle_nz, circle_g30, circle_max, circle_p99, circle_p95, circle_p90))

    dic = OrderedDict([('scale', []), ('id', []), ('hour', []),
                       ('clat', []), ('clon', []), ('lat_min', []), ('lat_max', []), ('lon_min', []), ('lon_max', []),
                       ('area', []),
                       ('bulk_pmax', []), ('bulk_pmean', []), ('bulk_tmean', []), ('bulk_tmean_p', []),
                       ('bulk_tmin_p', []), ('bulk_g30', []),
                       ('circle_Tcentre', []), ('circle_p', []), ('circle_t', []),
                       ('circle_val', []), ('circle_sum', []),
                       ('circle_nz', []), ('circle_g30', []), ('circle_max', []), ('circle_p99', []),
                       ('circle_p95', []), ('circle_p90', [])])

    keys = comp_collect.keys()
    print(keys)

    for v in res:

        print(v[2])

        comp_collect[v[2]]['p'].append(v[0])
        comp_collect[v[2]]['t'].append(v[1])
        comp_collect[v[2]]['hour'].append(v[4])
        comp_collect[v[2]]['id'].append(v[3])

        for cnt, kk in enumerate(dic.keys()):
            dic[kk].append(v[cnt + 2])  # omit kernel and kernelt

        precip[v[2]].extend(v[19])

    pkl.dump(dic, open(out + '3dmax_gt15000_blobs.p', 'wb'))

    pkl.dump(precip, open(out + 'precip_3dmax_gt15000_blobs.p', 'wb'))

    pkl.dump(comp_collect, open(out + 'comp_collect_composite_blobs.p', 'wb'))

    # df = pkl.load(open('/users/global/cornkle/C_paper/wavelet/saves/pandas/3dmax_gt15000.p', 'rb'))
    #
    # print(df['scale'])


def file_loop(fi):
    ret = []

    print('Doing file: ' + fi)

    dic = xr.open_dataset(fi)

    id = fi.split('/')[-1]

    outt = dic['tc_lag0'].values
    outp = dic['p'].values

    outt[np.isnan(outt)] = 150
    outt[outt >= -40] = 150
    grad = np.gradient(outt)
    outt[outt == 150] = np.nan
    outp[np.isnan(outt)] = np.nan


    area = np.sum(outt <= -40)
    try:
        bulk_pmax = np.max(outp[(np.isfinite(outp)) & (np.isfinite(outt))])
    except ValueError:
        return ret

    # 2004-07-17_04:09:00_80.nc (127) is a huge circular storm. good example for organisation
    #area gt 3000km2 cause that's 30km radius if circular
    if (area * 25 < 15000) or (area * 25 > 800000) or (bulk_pmax < 1) or (bulk_pmax > 200): #or (np.sum(np.isfinite(outp)) < (np.sum(np.isfinite(outt))*0.1)):
        print(area*25)
        print('throw out')
        return

    print('Area', area*25)

    perc = np.percentile(outt[np.isfinite(outt)], 60)  # 60

    clat = np.min(dic.lat.values) + ((np.max(dic.lat.values) - np.min(dic.lat.values)) * 0.5)
    clon = np.min(dic.lon.values) + ((np.max(dic.lon.values) - np.min(dic.lon.values)) * 0.5)

    lat_min = np.min(dic.lat.values)
    lat_max = np.max(dic.lat.values)
    lon_min = np.min(dic.lon.values)
    lon_max = np.max(dic.lon.values)

    bulk_tmean = np.nanmean(outt)
    lat = dic['lat'].values
    bulk_tmin_p = np.min(outt[(np.isfinite(outp)) & (np.isfinite(outt))])
    bulk_tmean_p = np.mean(outt[(np.isfinite(outp)) & (np.isfinite(outt))])
    bulk_pmean = np.max(outp[(np.isfinite(outp)) & (np.isfinite(outt))])
    bulk_g30 = np.sum(outp[(np.isfinite(outp)) & (np.isfinite(outt))] >= 30)


    o2 = outt.copy()
    o2[np.isnan(o2)] = perc
    nok = np.where(abs(grad[0]) > 80)
    d = 2
    i = nok[0]
    j = nok[1]

    for ii, jj in zip(i, j):
        kern = o2[ii - d:ii + d + 1, jj - d:jj + d + 1]
        o2[ii - d:ii + d + 1, jj - d:jj + d + 1] = ndimage.gaussian_filter(kern, 3, mode='nearest')

    wav = util.waveletT(o2, 5)
    wl = wav['t']  # [nb, :, :]

    # maxout = (
    #     wl == ndimage.maximum_filter(wl, (6, 6, 5), mode='constant',
    #                                  cval=np.amax(wl) + 1))  # (np.round(orig / 5))

    #(3, 3 ,len(wav['scales'])-1)
    arr = np.array(wav['scales'], dtype=str)

    scale_ind = range(arr.size)

    wlperc = wav['t'].copy()
    figure = np.zeros_like(outt)
    #
    # for s in range(wlperc.shape[0]):
    #     wlperc[s,:,:][wlperc[s,:,:] < np.percentile(wlperc[s,:,:][wlperc[s,:,:]>=0.05], 70)] = 0
    #
    # labels, numL = label(wlperc)
    # for s in [0,1]:
    #
    #      f = plt.figure()
    #      plt.imshow(labels[s,:,:])
    #
    # return

    yp, xp = np.where(outp > 30)

    size = [15,20,30,40,50,60,70,80,90,100,120,140,160,180,200,220]

    for sc in wav['scales'][::-1]:
        print(sc)

        mimin = np.argmin(np.abs(size-sc))


        idd = wav['scales'].tolist().index(sc)
        wll = wlperc[idd,:,:].copy()
      #  wll[wll < np.percentile(wll[wll >= 0.01], 90)] = 0
      #  wll[wll < np.percentile(np.arange(wll.min(), wll.max()+0.1,0.1), 5)] = 0
        wll[wll < np.percentile(wll[wll >= 0.05], 75)] = 0
      #  wll[(wll < sc**.5)] = 0


        labels, numL = label(wll)

        figure[np.nonzero(labels)] = size[mimin]
    figure[np.isnan(outt)]=0

    # pdb.set_trace()
    # figure[figure==0]=np.nan
    # f = plt.figure()
    # f.add_subplot(131)
    # plt.imshow(outt)
    # plt.imshow(figure, cmap='viridis')
    # f.add_subplot(132)
    # plt.imshow(figure, cmap='viridis')
    # plt.plot(xp, yp, 'yo', markersize=3)
    # f.add_subplot(133)
    # plt.imshow(outt)
    # plt.show()


    for sc in np.unique(figure):

        if sc == 0:
            continue

        perblob = figure.copy()
        perblob[perblob!=sc] = 0
        labels_blob, numL = label(perblob)

        for blob in np.unique(labels_blob):

            if blob == 0:
                continue

            pos = np.where(labels_blob == blob) #(figure == sc)

            if len(pos[0]) < 5:
                continue

            circle_p = outp[pos]
            circle_t = outt[pos]
            circle_valid = np.sum(np.isfinite(circle_p))

            # if no valid temperature pixels (noise from wavelet outside of storm) or very little valid TRMM pixels
            # (very small overlap of TRMM and MSG), continue!
            if (np.sum(np.isfinite(circle_t)) <= 0 ):
                continue

            if  ((circle_valid) < 3 ):
                continue

            circle_sum = np.nansum(circle_p)
            ## some rain at least
            if circle_sum < 0.1:
                continue

            circle_Tcenter = np.nanmin(outt[pos])

            y, x = np.where((outt == circle_Tcenter) & (labels_blob == blob) & (figure == sc))

            # figure[figure==0]=np.nan
            # f = plt.figure()
            # f.add_subplot(131)
            # plt.imshow(outt)
            # plt.imshow(figure, cmap='viridis')
            # f.add_subplot(132)
            # plt.imshow(figure, cmap='viridis')
            # plt.plot(xp, yp, 'yo', markersize=3)
            # plt.plot(pos[1], pos[0], 'ro', markersize=3)
            # f.add_subplot(133)
            # plt.imshow(outt)
            # plt.show()

            r = 20
            kernel = tm_utils.cut_kernel(outp, x[0], y[0], r)
            kernelt = tm_utils.cut_kernel(outt, x[0], y[0], r)

            if kernel.shape != (r * 2 + 1, r * 2 + 1):
                kernel = np.zeros((41, 41)) + np.nan

            if np.nansum(kernel) < 1:
                continue

            circle_nz = np.nansum(circle_p>0.1)
            circle_g30 = np.nansum(circle_p >= 30)

            try:
                circle_max = np.nanmax(circle_p)
            except ValueError:
                circle_max = np.nan
            try:
                circle_p99 = np.percentile(circle_p[circle_p>=0.1], 99)
            except IndexError:
                circle_p99 = np.nan
            try:
                circle_p95 = np.percentile(circle_p[circle_p>=0.1], 95)
            except IndexError:
                circle_p95 = np.nan
            try:
                circle_p90 = np.percentile(circle_p[circle_p>=0.1], 90)
            except IndexError:
                circle_p90 = np.nan


            #### HOW TO GIVE BACK THE MAX SCALE PER SYSTEM??

            ret.append((kernel, kernelt, sc, id, dic['time.hour'].values.tolist(),
                        clat, clon, lat_min, lat_max, lon_min, lon_max, area,
                        bulk_pmax, bulk_pmean, bulk_tmean, bulk_tmean_p, bulk_tmin_p, bulk_g30,
                        circle_Tcenter, circle_p, circle_t, circle_valid, circle_sum,
                        circle_nz, circle_g30, circle_max, circle_p99, circle_p95, circle_p90))


    figure[figure==0]=np.nan
    f = plt.figure()
    f.add_subplot(133)
    plt.imshow(outt, cmap='inferno')
    plt.imshow(figure, cmap='viridis')
    f.add_subplot(132)
    plt.imshow(figure, cmap='viridis')
    plt.colorbar()
    plt.plot(xp, yp, 'yo', markersize=3)
    f.add_subplot(131)
    plt.imshow(outt, cmap='inferno')
    plt.plot(xp, yp, 'yo', markersize=3)
    plt.show()
    # #
    #
    # f = plt.figure()
    # fcnt = 0
    # vv = 4
    # for s in scale_ind:
    #
    #     pos = np.where((maxout[s, :, :] == 1) & (outt <= -40))
    #
    #     if len(pos[0]) == 0:
    #         continue
    #     fcnt+=1
    #     ax = f.add_subplot(vv,vv,fcnt)
    #     ax.imshow(wl[s,:,:])
    #
    #    # plt.plot(xp, yp, 'yo', markersize=3)
    #
    #     ax.set_title(str(wav['scales'][s]))
    #
    # ax = f.add_subplot(vv, vv, fcnt+1)
    # plt.imshow(figure, cmap='viridis')
    # plt.plot(xp, yp, 'yo', markersize=3)
    # ax = f.add_subplot(vv, vv, fcnt + 2)
    # plt.imshow(outt)
    # ax = f.add_subplot(vv, vv, fcnt + 3)
    # plt.imshow(outp)
    # plt.plot(xp, yp, 'yo', markersize=3)
    # ax = f.add_subplot(vv, vv, fcnt + 4)
    # plt.imshow(figure, cmap='viridis')
    #
    #
    # plt.show()
    #
    # dic.close()

    return ret



if __name__ == "__main__":
    composite()







