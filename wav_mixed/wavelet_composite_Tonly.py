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
from CLOVER import era_geop_t3d as era_geop
from utils import u_gis

def composite():
    pool = multiprocessing.Pool(processes=7)
    files = ua.locate(".nc", '/users/global/cornkle/MCSfiles/WA15_big_-40_15W-20E_size_zR/')   # /WA30/
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

        # ret.append((kernel, kernelt, sc - 1, id, dic['time.hour'].values.tolist(),
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
                       ('circle_p95', []), ('circle_p90', []), ('circle_pc', [])])

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

    pkl.dump(dic, open(out + '3dmax_gt15000_TR_lax.p', 'wb'))

    pkl.dump(precip, open(out + 'precip_3dmax_gt15000_TR_lax.p', 'wb'))

    #pkl.dump(comp_collect, open(out + 'comp_collect_composite_T.p', 'wb'))

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
    outpc = dic['pconv'].values

    clat = np.min(dic.lat.values) + ((np.max(dic.lat.values) - np.min(dic.lat.values)) * 0.5)
    clon = np.min(dic.lon.values) + ((np.max(dic.lon.values) - np.min(dic.lon.values)) * 0.5)

    lat_min = np.min(dic.lat.values)
    lat_max = np.max(dic.lat.values)
    lon_min = np.min(dic.lon.values)
    lon_max = np.max(dic.lon.values)

    lon = dic['lon'].values
    lat = dic['lat'].values

    outt[np.isnan(outt)] = 150
    outt[outt >= -40] = 150

    outt[outt==150] = np.nan
    outp[np.isnan(outt)] = np.nan
    outpc[np.isnan(outt)] = np.nan

    area = np.sum(outt <= -40)
    try:
        bulk_pmax = np.max(outp[(np.isfinite(outp)) & (np.isfinite(outt))])
    except ValueError:
        return ret
    try:
        bulk_pmin = np.min(outp[(np.isfinite(outp)) & (np.isfinite(outt))])
    except ValueError:
        return ret

    # area gt 3000km2 cause that's 30km radius if circular
    if (area * 25 < 15000) or (area * 25 > 800000) or (
        bulk_pmax > 200) or (
        bulk_pmin < 0):  # or (np.sum(np.isfinite(outp)) < (np.sum(np.isfinite(outt))*0.1)): or (bulk_pmax < 1)
        print(area * 25)
        print('throw out')
        return ret

    bulk_pmean = np.max(outp[(np.isfinite(outp)) & (np.isfinite(outt))])
    bulk_g30 = np.sum(outp[(np.isfinite(outp)) & (np.isfinite(outt))] >= 30)
    lat = dic['lat'].values
    bulk_tmin_p = np.min(outt[(np.isfinite(outp)) & (np.isfinite(outt))])
    bulk_tmean_p = np.mean(outt[(np.isfinite(outp)) & (np.isfinite(outt))])
    bulk_tmean = np.nanmean(outt)

    outint = np.round(outt).astype(int)

    bins = np.arange(-95,-38,4)
    cent = (np.arange(-95,-38,4)-2)[1::]

    yp, xp = np.where(outp > 30)

    figure = np.zeros_like(outt)
    for idd, sc in enumerate(bins):

        if idd == 0:
            continue

        filter = np.where((outt>=bins[idd-1]) & (outt<sc))

        figure[filter] = cent[idd-1]
    #
    # f = plt.figure()
    # plt.imshow(figure)
    # plt.plot(xp, yp, 'yo', markersize=3)

    for sc in np.unique(figure):

        if sc == 0:
            continue

        perblob = figure.copy()
        perblob[perblob!=sc] = 0
        labels_blob, numL = label(perblob)

        for blob in np.unique(labels_blob):

            if blob == 0:
                continue

            pos = np.where((labels_blob == blob)) #(figure == sc)
            y = pos[0]
            x = pos[1]

            height = era_geop.era_Tlapse(int(dic['time.month']), sc, lon[y[0],x[0]], lat[y[0],x[0]])  # height in meters
            lx, ly = u_gis.parallax_corr_msg(0, 0,lon[y[0], x[0]], lat[y[0]],x[0]], height / 1000)

            lx = int(np.round(lx / 5.))
            ly = int(np.round(ly / 5.))  # km into pixels

            circle_p = outp[y-ly, x-lx]
            circle_pc = outpc[y-ly, x-lx]
            circle_t = outt[pos]
            circle_valid = np.sum(np.isfinite(circle_p))

            # if no valid temperature pixels (noise from wavelet outside of storm) or very little valid TRMM pixels
            # (very small overlap of TRMM and MSG), continue!
            if (np.sum(np.isfinite(circle_t)) <= 3 ):
                continue

            if  ((circle_valid) <= 3 ):
                continue

            circle_sum = np.nansum(circle_p)
            # # ## some rain at least
            # if circle_sum < 0.1:
            #     continue

            circle_Tcenter = np.nanmin(outt[pos])

            # y, x = np.where((outt == circle_Tcenter) & (labels_blob == blob))
            #
            # # figure[figure==0]=np.nan
            # # f = plt.figure()
            # # f.add_subplot(131)
            # # plt.imshow(outt)
            # # plt.imshow(figure, cmap='viridis')
            # # f.add_subplot(132)
            # # plt.imshow(figure, cmap='viridis')
            # # plt.plot(xp, yp, 'yo', markersize=3)
            # # plt.plot(pos[1], pos[0], 'ro', markersize=3)
            # # f.add_subplot(133)
            # # plt.imshow(outt)
            # # plt.show()
            #
            # r = 20
            # kernel = tm_utils.cut_kernel(outp, x[0], y[0], r)
            # kernelt = tm_utils.cut_kernel(outt, x[0], y[0], r)
            #
            # if kernel.shape != (r * 2 + 1, r * 2 + 1):
            #     kernel = np.zeros((41, 41)) + np.nan
            #
            # if np.nansum(kernel) < 1:
            #     continue

            kernel = 0
            kernelt = 0

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

            ret.append((kernel, kernelt, sc-1, id, dic['time.hour'].values.tolist(),
                        clat, clon, lat_min, lat_max, lon_min, lon_max, area,
                        bulk_pmax, bulk_pmean, bulk_tmean, bulk_tmean_p, bulk_tmin_p, bulk_g30,
                        circle_Tcenter, circle_p, circle_t, circle_valid, circle_sum,
                        circle_nz, circle_g30, circle_max, circle_p99, circle_p95, circle_p90, circle_pc))


    # figure[figure==0]=np.nan
    # f = plt.figure()
    # f.add_subplot(133)
    # plt.imshow(outt, cmap='inferno')
    # plt.imshow(figure, cmap='viridis')
    # f.add_subplot(132)
    # plt.imshow(figure, cmap='viridis')
    # plt.colorbar()
    # plt.plot(xp, yp, 'yo', markersize=3)
    # f.add_subplot(131)
    # plt.imshow(outt, cmap='inferno')
    # plt.plot(xp, yp, 'yo', markersize=3)
    # plt.show()
    #
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







