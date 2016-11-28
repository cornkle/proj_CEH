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
import multiprocessing
import ipdb
import pandas as pd


def composite():
    pool = multiprocessing.Pool(processes=7)
    files = ua.locate(".nc", '/users/global/cornkle/MCSfiles/WA15_big_-40_15W-20E/')   # /WA30/
    out = '/users/global/cornkle/C_paper/wavelet/saves/pandas/'
    #files = files[0:100]
    print('Nb files', len(files))
    tt = 'WA15'

    comp_collect = {}

    res = pool.map(file_loop, files)
    pool.close()
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
        comp_collect[v[0]]={'big': [], 'fin' : [], 'mean':[], 'isnz':[]}


    dic = {'scale' : [], 'sum' : [], 'sumvalid' : [], 'shape' : [],  'tmin' :[],  'sumnz':[], 'hour':[],
           'clat':[], 'clon' : [], 'area' : [], 'sum30' : [], 'max_sc' : [], 'tmean' : [] }  #  big , fin,shape, sum, sumvalid, tmin

    keys = comp_collect.keys()
    print(keys)

    for v in res:
        comp_collect[v[0]]['big'].append(v[1])
        comp_collect[v[0]]['fin'].append(v[2])
        comp_collect[v[0]]['mean'].append(v[4])
        comp_collect[v[0]]['isnz'].append(v[5])

        dic['scale'].append(v[0])
        dic['shape'].append(v[3])
        dic['sum'].append(v[6])
        dic['sumvalid'].append(v[7])
        dic['tmin'].append(v[8])
        dic['sumnz'].append(v[9])
        dic['hour'].append(v[10])
        dic['clat'].append(v[11])
        dic['clon'].append(v[12])
        dic['area'].append(v[13])
        dic['sum30'].append(v[14])
        dic['max_sc'].append(v[15])
        dic['tmean'].append(v[16])


    for k in keys:
        a = np.asarray(comp_collect[k]['big'])
        comp_collect[k]['big'] = a
        b = np.asarray(comp_collect[k]['fin'])
        comp_collect[k]['fin'] = b
        d = np.asarray(comp_collect[k]['mean'])
        comp_collect[k]['mean'] = d
        e = np.asarray(comp_collect[k]['isnz'])
        comp_collect[k]['isnz'] = e


    f = plt.figure()
    siz = 3

    keys = comp_collect.keys()
    print(keys)

    df = pd.DataFrame(dic)
   # df.to_pickle(out+'3dmax_gt15000.pkl')

    #return comp_collect

    ######### 2d plots
    ll = [18,30,60,90,101]#keys
    for ind, k in enumerate(ll):
        num = len(ll)
        arr = comp_collect[k]['mean']
        fin = comp_collect[k]['fin']
        big = comp_collect[k]['big']
        nz = comp_collect[k]['isnz']

        bla = np.nansum(arr, 0) / np.nansum(nz, 0)
        bla1 = np.nansum(fin, 0)
        blab = (np.nansum(big, 0) / np.nansum(nz, 0)) *100

        ax = f.add_subplot(3, num, 1 + ind)
        plt.imshow(bla, vmin=0, vmax=5, cmap='viridis')
        plt.title(str(k) + ' km', fontsize=9)
        plt.plot(20, 20, 'ro', markersize=siz)
        cbar = plt.colorbar()
        cbar.set_label('mm h-1')

        ax = f.add_subplot(3, num, 6 + ind)
        plt.imshow(bla1, cmap='viridis')
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
    plt.show()

    col = ['r', 'b', 'g', 'y', 'black']

    pplot=[]
    ll = keys
    for k in ll:

        pos = np.where(np.array(dic['scale'])==k)

        csum = np.array(dic['sum'])[pos]
        csumval = np.array(dic['sumvalid'])[pos]
        csumnz = np.array(dic['sumnz'])[pos]
        shape = np.array(dic['shape'])[pos]

        allsum = np.nansum(csum)
        allval = np.nansum(csumval)
        allnz = np.nansum(csumnz)
        sshape = np.nansum(shape)

        circmean = allsum / allval
        circnz = allsum / allnz

        pplot.append((k, sshape,  circmean, allsum, allval, circnz, allnz))

    f = plt.figure()
    for p in pplot:
        plt.scatter(p[0], p[1])
        plt.xlabel('T scales')
        plt.ylabel('nb in composite')
        plt.title(tt)
    plt.show()

    f = plt.figure()
    for p in pplot:
        plt.scatter(p[0], p[2] )
        plt.xlabel('T scales')
        plt.ylabel('Radial mean Pcp')
        plt.title(tt)
    plt.show()

    f = plt.figure()

    for p in pplot:
        plt.scatter(p[0], p[3])
        plt.xlabel('T scales')
        plt.ylabel('Radial Pcp sum')
        plt.title(tt)
        plt.show()

    f = plt.figure()
    for p in pplot:
        plt.scatter(p[0], p[4])
        plt.xlabel('T scales')
        plt.ylabel('Radial valid sum')
        plt.title(tt)
        plt.show()

    f = plt.figure()
    for p in pplot:
        plt.scatter(p[0], p[6])
        plt.xlabel('T scales')
        plt.ylabel('Radial valid nz sum')
        plt.title(tt)
        plt.show()

    f = plt.figure()
    for p in pplot:
        plt.scatter(p[0], p[5])
        plt.xlabel('T scales')
        plt.ylabel('Radial mean P nozero')
        plt.title(tt)
        plt.show()





def file_loop(f):
    ret = []
    max_sc = 202
    print('Doing file: ' + f)

    dic = xr.open_dataset(f)

    outt = dic['tc_lag0'].values
    outp = dic['p'].values

    tmean = np.nanmean(outt)

    outp[np.isnan(outt)]=np.nan

    area = np.sum(outt < -40)
    lat = dic['lat'].values
    tt = np.min(outt[(np.isfinite(outp)) & (np.isfinite(outt))])
    pp = np.max(outp[(np.isfinite(outp)) & (np.isfinite(outt))])
    clat = np.min(dic.lat.values) + ((np.max(dic.lat.values) - np.min(dic.lat.values)) * 0.5)
    clon = np.min(dic.lon.values) + ((np.max(dic.lon.values) - np.min(dic.lon.values)) * 0.5)

    #area gt 3000km2 cause that's 30km radius if circular
    if (area * 25 < 15000) or (area * 25 > 500000) or (pp < 1) or (pp > 200): #or (np.sum(np.isfinite(outp)) < (np.sum(np.isfinite(outt))*0.1)):
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

    arr = np.array(wav['scales'], dtype=str)

    scale_ind = range(arr.size)  # [0, 12, 21, -1, -13]

    for nb in scale_ind:

        ar = []
        msum = []
        msum_valid=[]
        msum_nz = []
        msum_30 = []

        ttmin = []

        orig = float(arr[nb])
        scale = int(np.round(orig))

        wl100 = wav['t'][-13, :, :]   # just using points in the middle

        wl = wav['t'][nb, :, :]

        maxoutt = (
        wl == ndimage.maximum_filter(wl, (5,5) , mode='constant', cval=np.amax(wl) + 1))  #(np.round(orig / 5))

        yp, xp = np.where(outp > 30)

        try:
            yy, xx = np.where((maxoutt == 1) & (outt <= -40) & (wl > np.percentile(wl[wl>=0.1], 80) ) & (wl100 > 5) )# & (wlperc > np.percentile(wlperc[wlperc>=0.1], 80) ))  # & (wl100 > 5)
        except IndexError:
            continue
        #print(scale, yy,xx)
        if len(yy)==1:
            max_sc = orig

        # if (scale == 17) or (scale == 18) or (scale == 25) or (scale == 85):
        #     inds = np.ceil(orig /5/2)
        #     f = plt.figure()
        #     siz = 3
        #     ax = f.add_subplot(1, 3, 1)
        #     plt.imshow(wav['t'][nb, :, :], cmap='jet')
        #     plt.plot(xp, yp, 'yo', markersize=siz)
        #     plt.plot(xx, yy, 'ro', markersize=siz)
        #     ax.set_title(str(scale), fontsize=12)
        #
        #     ax = f.add_subplot(1, 3, 2)
        #     plt.imshow(outt, cmap='jet')
        #     plt.plot(xp, yp, 'yo', markersize=siz)
        #     plt.plot(xx, yy, 'ro', markersize=siz+1)
        #  #   plt.plot(txx, tyy, 'go', markersize=siz)
        #     ax.set_title(str(scale), fontsize=12)
        #
        #     ax = f.add_subplot(1, 3, 3)
        #     plt.imshow(outp, cmap='jet')
        #     #plt.plot(xp, yp, 'yo', markersize=siz)
        #     plt.plot(xx, yy, 'ro', markersize=siz)
        #     plt.plot(xx+inds, yy, 'bo', markersize=siz)
        #     plt.plot(xx - inds, yy, 'bo', markersize=siz)
        #     plt.plot(xx, yy+inds, 'bo', markersize=siz)
        #     plt.plot(xx, yy - inds, 'bo', markersize=siz)
        #     ax.set_title(str(scale), fontsize=12)
        #     plt.show()

        yyy = []
        xxx = []
        for y, x in zip(yy, xx):

            r = 20
            kernel = tm_utils.cut_kernel(outp, x, y, r)

            if kernel.shape != (r * 2 + 1, r * 2 + 1):
                continue

            if np.nansum(kernel) < 1: #== 0:
                continue

            ar.append(kernel)
        if ar == []:
            continue

        for y, x in zip(yy, xx):

            ss = orig
            #ss = 15   # just looking at a fixed radius surrounding points defined by wavelet
            iscale = (np.ceil(ss / 2. / 5.)).astype(int)
            tmin = outt[y, x]

            ycirc, xcirc = ua.draw_cut_circle(x, y, iscale, outp)

            psum = np.nansum(outp[ycirc, xcirc])

            psum_valid = np.sum(np.isfinite(outp[ycirc, xcirc]))

            if psum_valid < 3:
                continue
            #
            # #only take circles that are 80% overlap
            # if psum_valid < len(outp[ycirc, xcirc])*0.8:
            #     continue

            # if (scale == 21) or  (scale == 202):
            #     bla = np.zeros_like(outp)
            #
            #     bla[ycirc, xcirc]=100
            #     bla[np.isnan(outp)]=0
            #     plt.figure()
            #     plt.title(str(scale))
            #     plt.imshow(bla)
            #     plt.show()


            xxx.append(x)
            yyy.append(y)

            psum_valid_nozero = np.sum(outp[ycirc,xcirc]>0)
            psum_30 = np.sum(outp[ycirc, xcirc] >= 30)


            #circle
            msum.append(psum)
            msum_valid.append(psum_valid)
            msum_nz.append(psum_valid_nozero)
            ttmin.append(tmin)
            msum_30.append(psum_30)


        # if (scale == 21) or (scale == 64) or (scale == 202) or (scale == 101):
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

        ar = np.array(ar)
        msum = np.sum(np.array(msum))
        msum_valid = np.sum(np.array(msum_valid))
        msum_nz = np.sum(np.array(msum_nz))
        msum_30 = np.sum(np.array(msum_30))

        isbigger = ar >= 30
        isfinite = np.isfinite(ar)
        isnozero = ar > 0

        isbig = np.nansum(isbigger, 0)
        isfin = np.nansum(isfinite, 0)
        nmean = np.nansum(ar,0)
        isnz = np.nansum(isnozero, 0)



        #### HOW TO GIVE BACK THE MAX SCALE PER SYSTEM??

        ret.append((scale, isbig, isfin, ar.shape[0],nmean, isnz, msum, msum_valid, tmin, msum_nz, dic['time.hour'].values.tolist(), clat, clon, area, msum_30, max_sc, tmean))

        dic.close()

    # print(max_sc)
    # print(area)
    return ret


if __name__ == "__main__":
    composite()

