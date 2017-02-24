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
import pandas as pd
import pickle as pkl
matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)


def composite():
    pool = multiprocessing.Pool(processes=7)
    files = ua.locate(".nc", '/users/global/cornkle/MCSfiles/WA15_big_-40_15W-20E_size/')   # /WA30/
    out = '/users/global/cornkle/C_paper/wavelet/saves/pandas/'
    #files = files[0:800]
    print('Nb files', len(files))
    tt = 'WA15'

    comp_collect = {}
    precip = {}

    res = pool.map(file_loop, files)
    pool.close()

    ipdb.set_trace()

    res = [x for x in res if x is not None]

    nb_sys = len(res)


    print('Number systems: ', nb_sys)

    res = [item for sublist in res for item in sublist] # flatten list of lists

    # for iis in range(1, iscale):
    #     yrcirc, xrcirc = ua.draw_ring(x, y, iis, iis + 1, outp)
    #     psum_r = np.nansum(outp[yrcirc, xrcirc])
    #     psum_valid_r = np.sum(np.isfinite(outp[yrcirc, xrcirc]))
    #     psum_valid_nozero_r = np.sum(outp[yrcirc, xrcirc] > 0.1)
    #     pmax_circ.append(psum_r / psum_valid_r)
    #     pmaxnz_circ.append(psum_r / psum_valid_nozero_r)
    #     maxscale.append(iis + 1)
    #
    #     # tarr = np.zeros_like(outp)
    #     # tarr[yrcirc,xrcirc]=1
    #     # if (scale == 15) or (scale == 64) or (scale == 202) or (scale == 101):
    #     #     plt.figure()
    #     #     plt.imshow(tarr)
    #     #     plt.title(str(iscale)+' '+str(orig))
    #     #     plt.plot(x, y, 'bo', markersize=1)
    #     # ring
    #     try:
    #         pmax_ind = np.nanargmax(pmax_circ)
    #         pmax.append(pmax_circ[pmax_ind])
    #         dist.append(maxscale[pmax_ind])
    #     except ValueError:
    #         pmax.append(np.nan)
    #         dist.append(np.nan)
    #     try:
    #         pmaxnz_ind = np.nanargmax(pmaxnz_circ)
    #         pmax_nz.append(pmaxnz_circ[pmaxnz_ind])
    #         distnz.append(maxscale[pmaxnz_ind])
    #     except ValueError:
    #         pmax_nz.append(np.nan)
    #         distnz.append(np.nan)

    #return res

    for v in res:
        comp_collect[v[0]]={'p': [], 't' : [], 'id' : [], 'hour': []}
        precip[v[0]]=[]

    #ret.append((scale, id, kernel, kernelt, tmin, p, area, b_tmin, b_pmax, b_tmean, clat, clon, lat_min, lat_max, lon_min, lon_max, dic['time.hour']))
    dic = {'scale' : [], 'id' : [], 'tmin' : [], 'p' : [], 'b_tmin':[], 'b_pmax':[], 'b_tmean':[], 'hour':[],
           'clat':[], 'clon' : [], 'lat_min' : [], 'lat_max':[], 'lon_min':[], 'lon_max':[], 'area' : [], 'pix_per_c' : []}


    keys = comp_collect.keys()
    print(keys)

    for v in res:

        print(v[0])
     #   ipdb.set_trace()
        comp_collect[v[0]]['p'].append(v[2])
        comp_collect[v[0]]['t'].append(v[3])
        comp_collect[v[0]]['id'].append(v[1])
        comp_collect[v[0]]['hour'].append(v[16])

        dic['scale'].append(v[0])
        dic['id'].append(v[1])
        dic['tmin'].append(v[4])
        dic['p'].append(v[5])
        dic['area'].append(v[6])
        dic['b_tmin'].append(v[7])
        dic['b_pmax'].append(v[8])
        dic['b_tmean'].append(v[9])
        dic['clat'].append(v[10])
        dic['clon'].append(v[11])
        dic['lat_min'].append(v[12])
        dic['lat_max'].append(v[13])
        dic['lon_min'].append(v[14])
        dic['lon_max'].append(v[15])
        dic['hour'].append(v[16])

        precip[v[0]].extend(v[5])

    pkl.dump(dic, open(out + '3dmax_gt15000_percircle.p', 'wb'))

    pkl.dump(precip, open(out+'precip_3dmax_gt15000_percircle_fake.p','wb'))

    pkl.dump(comp_collect, open(out + 'comp_collect_composite_percircle_fake.p', 'wb'))


def file_loop(fi):
    ret = []

    print('Doing file: ' + fi)

    dic = xr.open_dataset(fi)

    id = fi.split('/')[-1]

    outt = dic['tc_lag0'].values
    outp = dic['p'].values

    b_tmean = np.nanmean(outt)

    outp[np.isnan(outt)]=np.nan

    area = np.sum(outt < -40)
    lat_min = np.min(dic['lat'].values)
    lat_max = np.max(dic['lat'].values)
    lon_min = np.min(dic['lon'].values)
    lon_max = np.max(dic['lon'].values)
    b_tmin = np.min(outt[(np.isfinite(outp)) & (np.isfinite(outt))])
    b_pmax = np.max(outp[(np.isfinite(outp)) & (np.isfinite(outt))])
    clat = np.min(dic.lat.values) + ((np.max(dic.lat.values) - np.min(dic.lat.values)) * 0.5)
    clon = np.min(dic.lon.values) + ((np.max(dic.lon.values) - np.min(dic.lon.values)) * 0.5)

    #area gt 3000km2 cause that's 30km radius if circular
    if (area * 25 < 15000) or (area * 25 > 800000) or (b_pmax < 1) or (b_pmax > 200): #or (np.sum(np.isfinite(outp)) < (np.sum(np.isfinite(outt))*0.1)):
        print(area*25)
        print('throw out')
        return

    perc = np.percentile(outt[np.isfinite(outt)], 60)  # 60

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
    wl = wav['t']

    maxout = (
        wl == ndimage.maximum_filter(wl, (5, 5, 5), mode='constant', cval=np.amax(wl) + 1))  # (np.round(orig / 5))

    arr = np.array(wav['scales'], dtype=str)

    scale_ind = range(arr.size)

    # ######shuffled rainfall for test purposes!!
    #
    # pos = np.where(np.isfinite(outp))
    # shuffle = outp[pos]
    # np.random.shuffle(shuffle)
    # outp[pos] = shuffle
    #
    # ######

    for nb in scale_ind:

        orig = float(arr[nb])
        scale = int(np.round(orig))

        maxoutt=maxout[nb, :, :]

        yp, xp = np.where(outp > 30)

        try:
            yy, xx = np.where((maxoutt == 1) & (outt <= -40))#  & (wlperc > orig**.5))# & (wlperc > np.percentile(wlperc[wlperc>=0.1], 80)))# & (wlperc > np.percentile(wlperc[wlperc>=0.1], 80) ))  # & (wl100 > 5)
        except IndexError:
            continue
        print(scale, yy,xx)
        # if len(yy)==1:
        #     max_sc = orig

        # if (len(yy) != 0) & (area == 3177): # & (scale == 15): # (len(yy) != 0)
        #     inds = np.ceil(orig /5/2)
        #     f = plt.figure(figsize=(10,3))
        #     siz = 3
        #     ax = f.add_subplot(1, 3, 1)
        #     plt.imshow(wav['t'][nb, :, :], cmap='jet')
        #     plt.plot(xp, yp, 'yo', markersize=siz, label='Intense rain')
        #     plt.plot(xx, yy, 'ro', markersize=siz, label='Wavelet max')
        #     ax.set_title(str(scale)+'km: wavelet power', fontsize=10)
        #     plt.legend(fontsize=10)
        #
        #     ax = f.add_subplot(1, 3, 2)
        #     plt.imshow(outt, cmap='jet')
        #     plt.plot(xp, yp, 'yo', markersize=siz)
        #     plt.plot(xx, yy, 'ro', markersize=siz+1)
        #     plt.colorbar()
        #  #   plt.plot(txx, tyy, 'go', markersize=siz)
        #     ax.set_title( 'MSG cloud top temperature', fontsize=10)#
        #
        #     ax = f.add_subplot(1, 3, 3)
        #     plt.imshow(outp, cmap='jet')
        #     #plt.plot(xp, yp, 'yo', markersize=siz)
        #     plt.plot(xx, yy, 'ro', markersize=siz)
        #     plt.plot(xx+inds, yy, 'go', markersize=siz)
        #     plt.plot(xx - inds, yy, 'go', markersize=siz)
        #     plt.plot(xx, yy+inds, 'go', markersize=siz)
        #     plt.plot(xx, yy - inds, 'go', markersize=siz)
        #     ax.set_title('TRMM rainfall (cropped)', fontsize=10)
        #     plt.colorbar()
        #     plt.tight_layout()
        #     plt.savefig('/users/global/cornkle/C_paper/wavelet/figs/chris_presi/wav_example'+str(scale)+'_'+str(area)+'.pdf')
        #     plt.close('all')
        # #     #plt.show()



        for y, x in zip(yy, xx):

            #### kernel

            r = 20
            kernel = tm_utils.cut_kernel(outp, x, y, r)
            kernelt = tm_utils.cut_kernel(outt, x, y, r)

            if kernel.shape != (r * 2 + 1, r * 2 + 1):
                kernel = np.zeros((41,41))+ np.nan

            if np.nansum(kernel) < 1:
                continue

            #### circle

            ss = orig
            #ss = 30   # just looking at a fixed radius surrounding points defined by wavelet
            iscale = (np.ceil(ss / 2. / 5.)).astype(int)
            tmin = outt[y, x]

            ycirc, xcirc = ua.draw_cut_circle(x, y, iscale, outp)

            p = outp[ycirc, xcirc]
            psum_valid = np.sum(np.isfinite(p))

            if (psum_valid < 3):   # or (tmin > -70):
                continue

            # if (scale == 57) or (scale == 80) or (scale == 107) or (scale == 21):
            #     f = plt.figure()
            #     siz = 5
            #     ax = f.add_subplot(1, 3, 1)
            #     plt.imshow(wav['t'][nb, :, :], cmap='jet')
            #     plt.plot(xp, yp, 'yo', markersize=siz)
            #     plt.plot(xxx, yyy, 'ro', markersize=siz)
            #     ax.set_title('k'+str(scale), fontsize=12)
            #
            #     ax = f.add_subplot(1, 3, 2)
            #     plt.imshow(outt, cmap='jet')
            #     plt.plot(xp, yp, 'yo', markersize=siz)
            #     plt.plot(xxx, yyy, 'ro', markersize=siz)
            #     ax.set_title('k'+str(scale), fontsize=12)
            #
            #     ax = f.add_subplot(1, 3, 3)
            #     plt.imshow(outp, cmap='jet')
            #     #plt.plot(xp, yp, 'yo', markersize=siz)
            #     plt.plot(xxx, yyy, 'ro', markersize=siz)
            #     ax.set_title('k'+str(scale), fontsize=12)
            #     plt.show()

            kernel = np.array(kernel)
            kernelt = np.array(kernelt)

            ar = np.array(p)
            ar_t = np.array(ar_t)
            pvalues = np.array(pval)
            msum = np.nansum(msum)
            msum_valid = np.nansum(np.array(msum_valid))
            msum_nz = np.nansum(np.array(msum_nz))
            msum_nf = np.nansum(np.array(msum_nf))
            msum_30 = np.nansum(np.array(msum_30))
            try:
                pmax = np.nanmax(pvalues)
            except ValueError:
                pmax = np.nan
            try:
                p99 = np.percentile(pvalues, 99)
            except ValueError:
                p99 = np.nan


            pix_per_c = len(ycirc)



            ret.append((scale, id, kernel, kernelt, tmin, p, area, b_tmin, b_pmax, b_tmean, clat, clon, lat_min, lat_max, lon_min, lon_max, dic['time.hour'].values.tolist()), pix_per_c)

    dic.close()

    return ret



if __name__ == "__main__":
    composite()







