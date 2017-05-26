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
from utils import u_plot
import matplotlib as mpl
import multiprocessing
import skimage
import ipdb
import pdb
from collections import OrderedDict
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import pandas as pd
import pickle as pkl
matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)


def composite():
    pool = multiprocessing.Pool(processes=7)
    files = ua.locate(".nc", '/users/global/cornkle/MCSfiles/WA15_big_-40_15W-20E_size_zR/')   # /WA30/
    out = '/users/global/cornkle/C_paper/wavelet/saves/pandas/'
    #files = files[0:400]
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

    # ret.append((kernel, kernelt, sc, id, dic['time.hour'].values.tolist(),
    #             clat, clon, lat_min, lat_max, lon_min, lon_max, area,
    #             bulk_pmax, bulk_pmean, bulk_tmean, bulk_tmean_p, bulk_tmin_p, bulk_g30,
    #             circle_Tcenter, circle_p, circle_t, circle_valid, circle_sum,
    #             circle_nz, circle_g30, circle_max, circle_p99, circle_p95, circle_p90))

    dic = OrderedDict([('scale', []), ('id' , []), ('hour' , []),
           ('clat',[]), ('clon',[]),('lat_min',[]), ('lat_max' , []), ('lon_min' , []), ('lon_max' , []), ('area' , []),
           ('bulk_pmax' , []), ('bulk_pmean' ,[]), ('bulk_tmean',[]), ('bulk_tmean_p',[]), ('bulk_tmin_p',[]), ('bulk_g30',[]),
           ('circle_pix' , []), ('circle_Tcentre', []), ('circle_p' , []), ('circle_t' , []), ('circle_val' , []), ('circle_sum' , []),
           ('circle_nz' , []), ('circle_g30' , []), ('circle_max' , []), ('circle_p99' , []), ('circle_p95' , []), ('circle_p90' , []), ('circle_val_all', []), ('circle_pc', [])])

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


    pkl.dump(dic, open(out+'3dmax_gt15000_noR.p','wb'))

    pkl.dump(precip, open(out+'precip_3dmax_gt15000_noR.p','wb'))

    pkl.dump(comp_collect, open(out + 'comp_collect_composite_noR.p', 'wb'))

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

    outplot = outp.copy()
    #*0+1   ## ATTENTION CHANGED RAINFALL

    outt[np.isnan(outt)] = 150
    outt[outt >= -40] = 150
    grad = np.gradient(outt)
    outt[outt == 150] = np.nan
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

    #area gt 3000km2 cause that's 30km radius if circular
    if (area * 25 < 15000) or (area * 25 > 800000)  or (bulk_pmax > 200) or (bulk_pmin < 0): # or (bulk_pmax < 1) #or (np.sum(np.isfinite(outp)) < (np.sum(np.isfinite(outt))*0.1)):
        print(area*25)
        print('throw out')
        return

    perc = np.percentile(outt[np.isfinite(outt)], 60)  # 60

    clat = np.min(dic.lat.values) + ((np.max(dic.lat.values) - np.min(dic.lat.values)) * 0.5)
    clon = np.min(dic.lon.values) + ((np.max(dic.lon.values) - np.min(dic.lon.values)) * 0.5)

    if (clon > 28) or (clon < -17.2) or (clat < 4.1):
        return

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

    figure = np.zeros_like(outt)

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
    arrf = np.array(wav['scales'], dtype=float)

    scale_ind = range(arr.size)

    yp, xp = np.where(outp > 30)

    figure = np.zeros_like(outt)


    wll = wav['t']#[nb, :, :]
    maxs = np.zeros_like(wll)

    # maxoutt = (
    #     wll == ndimage.maximum_filter(wll, (1,5, 5), mode='reflect', cval=np.amax(wll) + 1))  # (np.round(orig / 5))

    yyy=[]
    xxx=[]
    scal=[]
    for nb in scale_ind[::-1]:

        orig = float(arr[nb])
        scale = int(np.round(orig))

        print(np.round(orig))

        wl = wll[nb, :, :]
        #maxout = maxoutt[nb, :, :]

        maxout = (
            wl == ndimage.maximum_filter(wl, (5,5), mode='constant', cval=np.amax(wl) + 1))  # (np.round(orig / 5))

        try:
            yy, xx = np.where((maxout == 1) & (outt <= -40)  & ((wl >= np.percentile(wl[wl >= 0.5], 90)) & (wl > orig**.5) )) #(wl >= np.percentile(wl[wl >= 0.5], 90)))# & (wl > orig**.5))#& (wl >= np.percentile(wl[wl >= 0.5], 90))) #)& (wl > orig**.5) (wl >= np.percentile(wl[wl >= 0.1], 90)) )#(wl > orig**.5))#  & (wlperc > orig**.5))# & (wlperc > np.percentile(wlperc[wlperc>=0.1], 80)))# & (wlperc > np.percentile(wlperc[wlperc>=0.1], 80) ))  # & (wl100 > 5)
        except IndexError:
            continue



        for y, x in zip(yy, xx):

            ss = orig
            iscale = (np.ceil(ss / 2. / 5.)).astype(int)
            if ss <= 20:
                iscale = iscale+1

            ycirc, xcirc = ua.draw_cut_circle(x, y, iscale, outp)

            figure[ycirc, xcirc] = np.round(orig)
            xxx.append(x)
            yyy.append(y)
            scal.append(orig)

    figure[np.isnan(outt)]=0

    circle_val_all = np.sum(figure>0)


    for y, x, sc in zip(yyy, xxx, scal):

        #print(sc, y, x)
        dummy = figure.copy()

        posi =  np.where(arrf == sc)

        int_sc = np.round(sc)
        radius = sc

        r = 20
        kernel = tm_utils.cut_kernel(outp, x, y, r)
        kernelt = tm_utils.cut_kernel(outt, x, y, r)
        kernelf = tm_utils.cut_kernel(figure, x, y, r)

        if kernel.shape != (r * 2 + 1, r * 2 + 1):
            kernel = np.zeros((41, 41)) + np.nan
        if kernelt.shape != (r * 2 + 1, r * 2 + 1):
            kernelt = np.zeros((41, 41)) + np.nan


        #pdb.set_trace()

        # if np.nansum(kernel) < 1:
        #     continue

        # npos = np.where(kernelf != int_sc)
        # kernel[npos] = np.nan
        # kernelt[npos] = np.nan
        #
        # f = plt.figure()
        # plt.imshow(kernelt, cmap='inferno')

        circle_Tcenter = outt[y, x]

        iscale = (np.ceil(radius / 2. / 5.)).astype(int)
        if ss <= 20:
            iscale = iscale + 1
        ycircf, xcircf = ua.draw_cut_circle(x, y, iscale, outp)


        pos = np.where((figure[ycircf, xcircf] == int_sc))

        if len(pos[0]) <= 3:
            continue

        # f = plt.figure()
        # dummy[ycircf[pos], xcircf[pos]] = 1000
        # plt.imshow(dummy)
        # plt.plot(xp, yp, 'ro', markersize=3)
        # plt.title(str(int_sc))

        circle_p = outp[ycircf[pos], xcircf[pos]]
        circle_pc = outpc[ycircf[pos], xcircf[pos]]
        circle_t = outt[ycircf[pos], xcircf[pos]]
        circle_valid = np.sum(np.isfinite(circle_p))


        if ((circle_valid) < 3 ):   # or (tmin > -70):
            continue

        # ## some rain
        # if np.nansum(circle_p) < 0.1:
        #     continue

        circle_sum = np.nansum(circle_p)
        circle_nz = np.nansum(circle_p>0.1)
        circle_g30 = np.nansum(circle_p >= 30)

        #print('circle_g30', circle_g30)

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

        maxs[posi, y, x] = 1

        #### HOW TO GIVE BACK THE MAX SCALE PER SYSTEM??

        #print('GOING BACK')

        ret.append((kernel, kernelt, int_sc, id, dic['time.hour'].values.tolist(),
                    clat, clon, lat_min, lat_max, lon_min, lon_max, area,
                    bulk_pmax, bulk_pmean, bulk_tmean, bulk_tmean_p, bulk_tmin_p, bulk_g30,
                    len(ycircf), circle_Tcenter, circle_p, circle_t, circle_valid, circle_sum,
                    circle_nz, circle_g30, circle_max, circle_p99, circle_p95, circle_p90, circle_val_all, circle_pc))

    # figure[figure == 0] = np.nan
    # f = plt.figure()
    # ax1 = f.add_subplot(131)
    # plt.imshow(outt)
    # plt.imshow(figure, cmap='viridis')
    #
    # ax2 = f.add_subplot(132)
    # plt.imshow(figure, cmap='viridis')
    # plt.plot(xp, yp, 'yo', markersize=3)
    # ax3 = f.add_subplot(133)
    # plt.imshow(outt)
    # ax1.invert_yaxis()
    # ax2.invert_yaxis()
    # ax3.invert_yaxis()
    #
    # plt.show()
    ##file 130!!! nR

    # spos = np.where(np.array(scal, dtype=int) == 15)
    # figure[figure == 0] = np.nan
    # f = plt.figure(figsize = (7,6), dpi=300)
    # ax2 = f.add_subplot(111)
    # #plt.imshow(outt, cmap='inferno')
    # lev = np.arange(-90, -39, 2)
    # ax2.imshow(outt, cmap='Greys', vmax=-40)
    # #plt.imshow(figure, cmap='viridis')
    # lev = arr
    # mt = ax2.imshow(figure, cmap='viridis')
    # plt.plot(np.array(xxx)[spos], np.array(yyy)[spos], 'wo', markersize=3, label='Wavelet power maximum')
    # #plt.plot(xp, yp, 'o', markersize=3)
    # ax2.invert_yaxis()
    # ax2.set_xlim(20, 140)
    # ax2.set_ylim(20, 140)
    # ax2.set_xticklabels(np.arange(100, 800, 100))
    # ax2.set_yticklabels(np.arange(100, 800, 100))
    # ax2.plot(xp , yp , 'o', markersize=3, label='Rain > 30mm h$^{-1}$')
    # ax2.set_xlabel('Spatial extent (km)')
    # ax2.set_ylabel('Spatial extent (km)')
    # # f.add_subplot(132)
    # # plt.imshow(figure, cmap='viridis')
    # # plt.plot(xp, yp, 'ro', markersize=3)
    # plt.colorbar(mt, label = 'Scale (km)')
    # #plt.legend(loc=4)
    # # ax1 = f.add_subplot(121)
    # # plt.imshow(outt, cmap='inferno')
    # # plt.plot(xp, yp, 'o', markersize=3)
    # # plt.plot(xxx, yyy, 'wo', markersize=3)
    # # ax1.invert_yaxis()
    #
    # plt.tight_layout()
    # plt.show()
    # spath = '/users/global/cornkle/C_paper/wavelet/figs/paper/'
    # plt.savefig(spath + '/method_circles2.png', dpi=300)
    #
    #
    #
    #from utils import u_arrays as ua
    #files = ua.locate(".nc", '/users/global/cornkle/MCSfiles/WA15_big_-40_15W-20E_size_zR/')
    #from wavelet import wavelet_composite_no_overlap as wcno
    #bla = wcno.file_loop(files[130])


    # f = plt.figure(figsize = (6.5,11), dpi=300)
    #
    # gridspec.GridSpec(4,1)
    # posi = 71#116 ## 118
    # ax1 = plt.subplot2grid((4,1),(0,0),rowspan=2)
    # ax2 = plt.subplot2grid((4,1),(2,0))
    # ax3 = plt.subplot2grid((4,1),(3,0))
    #
    # lev = np.arange(-90,-39,4)
    # ax1.contourf(np.arange(wll.shape[2]) * 5, np.arange(wll.shape[1]) * 5, outplot, cmap='Blues')
    # mt = ax1.contourf(np.arange(wll.shape[2])*5,np.arange(wll.shape[1])*5,outt, cmap='Greys', vmax=-40, levels = lev)
    # ax1.plot(np.arange(wll.shape[2])*5, [posi*5]*len(np.arange(wll.shape[2])*5), linestyle = '--', linewidth=2, color = 'black')
    # ax1.invert_yaxis()
    # ax1.set_xlim(100,700)
    # ax1.set_ylim(100, 700)
    #
    # ax1.plot(xp*5, yp*5, 'o', markersize=3, label='Rain > 30mm h$^{-1}$')
    # ax1.legend(loc=4)
    # ax1.set_ylabel('Spatial extent (km)')
    # ax1.set_title(str(dic['time.year'].values)+'-'+str(dic['time.month'].values)+'-'+str(dic['time.day'].values)+' '+str(dic['time.hour'].values)+':'+str(dic['time.minute'].values)+'UTC')
    #
    # colors = cm.viridis(np.linspace(0, 1, len([0,1, 2,5,10,20,40,60,80])))
    #
    # ax2.plot(np.arange(wll.shape[2])*5, outt[posi,:], color='r')  #118
    # ax2.set_xlim(100, 700)
    #
    # ax2.set_ylabel('Cloud-top temperature ($^{\circ}$C)')
    # ax22 = ax2.twinx()
    # ax22.set_xlim(100, 700)
    # ax22.plot(np.arange(wll.shape[2])*5, outp[posi,:])
    # ax22.set_ylabel('Rain (mm h$^{-1}$)')
    # print(np.nanmax(wll[:,posi,:]))
    # #arrt = np.array(arr, dtype=int)
    #
    # mp = ax3.contourf(np.arange(wll.shape[2])*5, arr,wll[:,posi,:], levels=[0,1, 2,5,10,20,40,80,100], colors=colors)
    # maxs = np.mean(maxs[:,posi-1:posi+2, :], 1) # -1, +2
    # #ax3.contour(np.arange(wll.shape[2])*5, arr,maxs, cmap='Greys_r')
    #
    # ppos = np.where(maxs)
    #
    # for p1, p2 in zip(ppos[1], ppos[0]):
    #     ax3.errorbar((np.arange(wll.shape[2])*5)[p1], arrf[p2], xerr=arrf[p2]/2, fmt='o', ecolor='white', color='white', capthick=3, ms=3, elinewidth=0.7)
    # ax3.set_xlim(100,700)
    # ax3.set_ylim(15, 180)
    # ax3.set_xlabel('Spatial extent (km)')
    # ax3.set_ylabel('Length scale (km)')
    #
    # plt.tight_layout()
    #
    # f.subplots_adjust(right=0.86)
    #
    # cax = f.add_axes([0.87, 0.545, 0.025, 0.415])
    # cb = plt.colorbar(mt, cax=cax, label='Cloud-top temperature ($^{\circ}$C)')
    # cb.ax.tick_params(labelsize=12)
    #
    # cax = f.add_axes([0.87, 0.065, 0.025, 0.175])
    # cb = plt.colorbar(mp, cax=cax, label='Wavelet power')
    # cb.ax.tick_params(labelsize=12)
    #
    # fsiz = 14
    # x = 0.02
    # plt.annotate('a)', xy=(x, 0.96), xytext=(0, 4), size=fsiz, xycoords=('figure fraction', 'figure fraction'),
    #              textcoords='offset points')
    # plt.annotate('b)', xy=(x, 0.51), xytext=(0, 4), size=fsiz, xycoords=('figure fraction', 'figure fraction'),
    #              textcoords='offset points')
    # plt.annotate('c)', xy=(x, 0.245), xytext=(0, 4), size=fsiz, xycoords=('figure fraction', 'figure fraction'),
    #              textcoords='offset points')
    #
    # plt.show()
    # spath = '/users/global/cornkle/C_paper/wavelet/figs/paper/'
    # plt.savefig(spath+'/method2.png', dpi=300)

    dic.close()

    #plt.close('all')

    return ret



if __name__ == "__main__":
    composite()







