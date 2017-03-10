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
import ipdb
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
    files = files[0:100]
    print('Nb files', len(files))
    tt = 'WA15'

    comp_collect = {}
    precip = {}

    res = pool.map(file_loop, files)
    pool.close()
    res = [x for x in res if x is not None]

    nb_sys = len(res)

    print('Number systems: ', nb_sys)

    res = [item for sublist in res for item in sublist] # flatten list of lists

    for v in res:

        comp_collect[v[2]]={'p': [], 't' : [], 'scale':[], 'hour':[], 'id' : []}
        precip[v[2]]=[]

     # ret.append((kernel, kernelt, scale, id, dic['time.hour'].values.tolist(),
     #                    clat, clon, lat_min, lat_max, lon_min, lon_max, area,
     #                    bulk_pmax, bulk_pmean, bulk_tmean, bulk_tmean_p, bulk_tmin_p, bulk_g30,
     #                    len(ycirc), circle_Tcenter, circle_p, circle_valid, circle_sum,
     #                    circle_nz, circle_g30, circle_max, circle_p99, circle_p95, circle_p90))

    dic = OrderedDict([('scale', []), ('id' , []), ('hour' , []),
           ('clat',[]), ('clon',[]),('lat_min',[]), ('lat_max' , []), ('lon_min' , []), ('lon_max' , []), ('area' , []),
           ('bulk_pmax' , []), ('bulk_pmean' ,[]), ('bulk_tmean',[]), ('bulk_tmean_p',[]), ('bulk_tmin_p',[]), ('bulk_g30',[]),
           ('circle_pix' , []), ('circle_Tcentre', []), ('circle_p' , []), ('circle_val' , []), ('circle_sum' , []),
           ('circle_nz' , []), ('circle_g30' , []), ('circle_max' , []), ('circle_p99' , []), ('circle_p95' , []), ('circle_p90' , [])])

    keys = comp_collect.keys()
    print(keys)

    for v in res:

        print(v[2])

        comp_collect[v[2]]['p'].append(v[0])
        comp_collect[v[2]]['t'].append(v[1])
        comp_collect[v[2]]['hour'].append(v[4])
        comp_collect[v[2]]['id'].append(v[3])

        for cnt, kk in enumerate(dic.keys()):

            dic[kk].append(v[cnt+2])  # omit kernel and kernelt

        precip[v[2]].extend(v[20])


    pkl.dump(dic, open(out+'3dmax_gt15000_blobs.p','wb'))

    pkl.dump(precip, open(out+'precip_3dmax_gt15000_blobs.p','wb'))

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
    #*0+1   ## ATTENTION CHANGED RAINFALL

    bulk_tmean = np.nanmean(outt)

    outp[np.isnan(outt)]=np.nan

    area = np.sum(outt <= -40)
    lat = dic['lat'].values
    bulk_tmin_p = np.min(outt[(np.isfinite(outp)) & (np.isfinite(outt))])
    bulk_tmean_p = np.mean(outt[(np.isfinite(outp)) & (np.isfinite(outt))])
    bulk_pmax = np.max(outp[(np.isfinite(outp)) & (np.isfinite(outt))])
    bulk_pmean = np.max(outp[(np.isfinite(outp)) & (np.isfinite(outt))])
    bulk_g30 = np.sum(outp[(np.isfinite(outp)) & (np.isfinite(outt))]>=30)

    clat = np.min(dic.lat.values) + ((np.max(dic.lat.values) - np.min(dic.lat.values)) * 0.5)
    clon = np.min(dic.lon.values) + ((np.max(dic.lon.values) - np.min(dic.lon.values)) * 0.5)

    lat_min = np.min(dic.lat.values)
    lat_max = np.max(dic.lat.values)
    lon_min = np.min(dic.lon.values)
    lon_max = np.max(dic.lon.values)


    #area gt 3000km2 cause that's 30km radius if circular
    if (area * 25 < 15000) or (area * 25 > 800000) or (bulk_pmax < 1) or (bulk_pmax > 200): #or (np.sum(np.isfinite(outp)) < (np.sum(np.isfinite(outt))*0.1)):
        print(area*25)
        print('throw out')
        return

    perc = np.percentile(outt[np.isfinite(outt)], 70)  # 60

    outt[np.isnan(outt)] = 150
    outt[outt >= -40] = 150
    grad = np.gradient(outt)
    outt[outt==150] = np.nan



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

    maxout = (
        wl == ndimage.maximum_filter(wl, (6, 6, 5), mode='constant',
                                     cval=np.amax(wl) + 1))  # (np.round(orig / 5))

    #(3, 3 ,len(wav['scales'])-1)
    arr = np.array(wav['scales'], dtype=str)

    scale_ind = range(arr.size)

    wlperc = wav['t'].copy()
    figure = np.zeros_like(outt)

    for s in range(wlperc.shape[0]):
        wlperc[s,:,:][wlperc[s,:,:] < np.percentile(wlperc[s,:,:][wlperc[s,:,:]>=0.05], 70)] = 0

    labels, numL = label(wlperc)

    yp, xp = np.where(outp > 30)

    size = np.array(np.arange(10,200, 10))

    for sc in wav['scales'][::-1]:
        print(sc)

        mimin = np.argmin(np.abs(size-sc))

        id = wav['scales'].tolist().index(sc)

        figure[np.nonzero(labels[id,:,:])] = size[mimin]
    figure[figure==0]=np.nan
    f = plt.figure()
    f.add_subplot(131)
    plt.imshow(outt)
    plt.imshow(figure, cmap='viridis')
    f.add_subplot(132)
    plt.imshow(figure, cmap='viridis')
    plt.plot(xp, yp, 'yo', markersize=3)
    f.add_subplot(133)
    plt.imshow(outt)
    plt.show()

    return

    for nb in scale_ind:


        orig = float(arr[nb])
        scale = int(np.round(orig))

        maxoutt=maxout[nb, :, :]

        yp, xp = np.where(outp > 30)

        try:
            yy, xx = np.where((maxoutt == 1) & (outt <= -40))#  & (wlperc > orig**.5))# & (wlperc > np.percentile(wlperc[wlperc>=0.1], 80)))# & (wlperc > np.percentile(wlperc[wlperc>=0.1], 80) ))  # & (wl100 > 5)
        except IndexError:
            continue

        #print(scale, yy,xx)
        # #Example why we miss loads of extreme rain:
        # # Doing file: /users/global/cornkle/MCSfiles/WA15_big_-40_15W-20E_size/2004-06-07_10:15:00_11.nc
        if (len(yy) != 0): # & (scale == 15): # (len(yy) != 0)
            inds = np.ceil(orig /5/2)
            f = plt.figure(figsize=(10,3))
            siz = 3
            ax = f.add_subplot(1, 3, 1)
            plt.imshow(wav['t'][nb, :, :], cmap='jet')
            plt.plot(xp, yp, 'yo', markersize=siz, label='Intense rain')
            plt.plot(xx, yy, 'ro', markersize=siz, label='Wavelet max')
            ax.set_title(str(scale)+'km: wavelet power', fontsize=10)
            plt.legend(fontsize=10)

            ax = f.add_subplot(1, 3, 2)
            plt.imshow(outt, cmap='jet')
            #plt.plot(xp, yp, 'yo', markersize=siz)
            plt.plot(xx, yy, 'ro', markersize=siz+1)
            plt.colorbar()
         #   plt.plot(txx, tyy, 'go', markersize=siz)
            ax.set_title( 'MSG cloud top temperature', fontsize=10)#

            ax = f.add_subplot(1, 3, 3)
            plt.imshow(outp, cmap='jet')
            #plt.plot(xp, yp, 'yo', markersize=siz)
            plt.plot(xx, yy, 'ro', markersize=siz)
            plt.plot(xx+inds, yy, 'go', markersize=siz)
            plt.plot(xx - inds, yy, 'go', markersize=siz)
            plt.plot(xx, yy+inds, 'go', markersize=siz)
            plt.plot(xx, yy - inds, 'go', markersize=siz)
            ax.set_title('TRMM rainfall (cropped)', fontsize=10)
            plt.colorbar()
            plt.tight_layout()
           # plt.savefig('/users/global/cornkle/C_paper/wavelet/figs/chris_presi/wav_example'+str(scale)+'_'+str(area)+'.png')
          #  plt.close('all')
        #     #plt.show()

        for y, x in zip(yy, xx):

            if figure[nb, y,x] == 0:
                continue

            r = 20
            kernel = tm_utils.cut_kernel(outp, x, y, r)
            kernelt = tm_utils.cut_kernel(outt, x, y, r)

            if kernel.shape != (r * 2 + 1, r * 2 + 1):
                kernel = np.zeros((41,41))+ np.nan

            if np.nansum(kernel) < 1:
                continue

            circle_Tcenter = outt[y, x]


            ycirc, xcirc = np.where(figure[nb, :, :] == figure[nb, y, x])

            circle_p = outp[ycirc, xcirc]

            circle_valid = np.sum(np.isfinite(circle_p))

            if ((circle_valid) < 3 ):   # or (tmin > -70):
                continue

            circle_sum = np.nansum(circle_p)
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

            ret.append((kernel, kernelt, scale, id, dic['time.hour'].values.tolist(),
                        clat, clon, lat_min, lat_max, lon_min, lon_max, area,
                        bulk_pmax, bulk_pmean, bulk_tmean, bulk_tmean_p, bulk_tmin_p, bulk_g30,
                        len(ycirc), circle_Tcenter, circle_p, circle_valid, circle_sum,
                        circle_nz, circle_g30, circle_max, circle_p99, circle_p95, circle_p90))


    f = plt.figure()
    fcnt = 0
    vv = 7
    for s in scale_ind:

        pos = np.where((maxout[s, :, :] == 1) & (outt <= -40))

        if len(pos[0]) == 0:
            continue
        fcnt+=1
        ax = f.add_subplot(vv,vv,fcnt)
        ax.imshow(wl[s,:,:])

       # plt.plot(xp, yp, 'yo', markersize=3)

        ax.set_title(str(wav['scales'][s]))

    ax = f.add_subplot(vv, vv, fcnt+1)
    plt.imshow(figure)
    plt.plot(xp, yp, 'yo', markersize=3)
    ax = f.add_subplot(vv, vv, fcnt + 2)
    plt.imshow(outt)
    ax = f.add_subplot(vv, vv, fcnt + 3)
    plt.imshow(outp)
    plt.plot(xp, yp, 'yo', markersize=3)
    ax = f.add_subplot(vv, vv, fcnt + 4)
    plt.imshow(figure)


    plt.show()

    dic.close()

    return ret



if __name__ == "__main__":
    composite()







