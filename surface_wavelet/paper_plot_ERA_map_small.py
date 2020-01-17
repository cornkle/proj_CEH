# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 10:15:40 2016

@author: cornkle
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib
import ipdb
import pandas as pd
from collections import OrderedDict
import salem
from utils import u_met, u_parallelise, u_arrays as ua, constants as cnst, u_darrays
from scipy.interpolate import griddata
import multiprocessing
import glob

import pickle as pkl


matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)

def plot_eh():
    for eh in range(-38,1,3):
        plot_doug_all(17,eh)


def plot_noon():
    afternoon = list(range(15,24))
    night = list(range(0,8))
    all = afternoon + night

    hlist = []
    for hh in all:
        if hh >= 15:
            hlist.append((hh,12-hh))
        else:
            hlist.append((hh, 12-(hh+24)))

    for h in hlist:
        plot_doug_all(h[0],h[1])



def plot_doug_timeseries():

    x = len(list(range(-38, 4, 3)))
    y = 401

    # outticks = list(range(-30, 1, 5))
    # ranges = np.arange(-30,1,3)
    #
    # outticks = [12,17,22,3,8,13,18]
    #outticks = [1, 6, 11, 16, 21, 2, 7, 12, 17]
    ranges = np.arange(-38, 4, 3)

    h=17

    outdic = {}
    dummyfile = glob.glob(cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/ERA5_composite_cores_*_small_cores.p")
    dummy = pkl.load(open(dummyfile[0], "rb"))

    for k in dummy.keys():
        outdic[k] = np.zeros((y, x))

    for ids, eh in enumerate(range(-38,4,3)):

        dic = {}

        def coll(dic, h, eh, year):
            print(h)
            core = pkl.load(open(
                cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/ERA5_composite_cores_LSTA_500w04_15k_"+str(eh) + "UTCERA"+str(h).zfill(2)+'_'+str(year)+"_small_cores.p", "rb"))
            for id, k in enumerate(core.keys()):
                try:
                    dic[k] = dic[k] + core[k]
                except KeyError:
                    dic[k] = core[k]


        for y in range(2006,2011):
            coll(dic, h, eh, y)

        for k in dic.keys():
            outdic[k][:,ids] = dic[k][:,190:211].mean(axis=1)

    print(outdic.keys())
    f = plt.figure(figsize=(15,10))
    ax = f.add_subplot(321)

    plt.contourf(ranges,np.arange(401), (outdic['q'])*1000/ outdic['cntp'],levels=np.linspace(-0.5,0.5,14), cmap='RdBu', extend='both')
    plt.colorbar(label=r'g kg$^{-1}$')
    contours = plt.contour(ranges, np.arange(401), (outdic['v650']) / outdic['cntp'], levels=np.linspace(-0.8,0.8,11), cmap='RdBu')
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    ax.set_yticklabels(np.array((np.linspace(0, 200 * 2, 9) - 200) * 3, dtype=int))
    #ax.set_xticklabels(outticks)

    plt.title('Shading: q-anomaly, contours: 650hpa v-wind anomaly')
    plt.xlabel('Hour relative to convective core at 1700UTC')
    plt.ylabel('North-South distance from core (km)')

    plt.hlines(200,xmin=ranges[0], xmax=ranges[-1], linestyle='dashed')


    ax = f.add_subplot(324)

    plt.contourf(ranges,np.arange(401), (outdic['t'] / outdic['cntp']), extend='both',levels=[ -0.8,-0.7, -0.6,-0.5,-0.4,-0.2,-0.1, 0.1,0.2,0.4,0.5,0.6, 0.7, 0.8], cmap='RdBu_r')
    plt.colorbar(label=r'K')
    contours = plt.contour(ranges, np.arange(401), (outdic['u650_orig'] / outdic['cntp']), extend='both', cmap='viridis', levels=np.arange(-15,-8,0.5)) #levels=np.arange(-15,-8,0.5)) #

    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    ax.set_yticklabels(np.array((np.linspace(0, 200 * 2, 9) - 200) * 3, dtype=int))
    #ax.set_xticklabels(outticks)

    plt.title('Shading: t-anomaly, contours: 650hPa u-wind')
    plt.hlines(200,xmin=ranges[0], xmax=ranges[-1], linestyle='dashed')
    plt.xlabel('Hour of day')
    plt.ylabel('North-South distance from core (km)')


    ax = f.add_subplot(323)
    plt.contourf(ranges,np.arange(401), (outdic['t'])/ outdic['cntp'], extend='both',  cmap='RdBu_r',levels=[ -0.8,-0.7, -0.6,-0.5,-0.4,-0.2,-0.1, 0.1,0.2,0.4,0.5,0.6, 0.7, 0.8])
    plt.colorbar(label=r'K')
    contours = plt.contour(ranges, np.arange(401), (outdic['u650'] / outdic['cntp']), extend='both', cmap='RdBu', levels=np.linspace(-0.7,0.7,9))
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    ax.set_yticklabels(np.array((np.linspace(0, 200 * 2, 9) - 200) * 3, dtype=int))
    #ax.set_xticklabels(outticks)

    plt.title('Shading:t-anomaly, contours: 650hpa u-wind anomaly')
    plt.hlines(200,xmin=ranges[0], xmax=ranges[-1], linestyle='dashed')
    plt.xlabel('Hour of day')
    plt.ylabel('North-South distance from core (km)')



    ax = f.add_subplot(322)
    plt.contourf(ranges,np.arange(401), (outdic['q'])*1000/ outdic['cntp'],levels=np.linspace(-0.5,0.5,14), cmap='RdBu', extend='both')
    plt.colorbar(label=r'g kg$^{-1}$')
    contours = plt.contour(ranges, np.arange(401), (outdic['u925'] / outdic['cntp']), extend='both', cmap='RdBu', levels=np.linspace(-0.7,0.7,9))
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    ax.set_yticklabels(np.array((np.linspace(0, 200 * 2, 9) - 200) * 3, dtype=int))
    #ax.set_xticklabels(outticks)
    plt.title('Shading: q-anomaly, contours: 925hpa u-wind anomaly')
    plt.hlines(200,xmin=ranges[0], xmax=ranges[-1], linestyle='dashed')
    plt.xlabel('Hour of day')
    plt.ylabel('North-South distance from core (km)')

    ax = f.add_subplot(325)
    plt.contourf(ranges,np.arange(401), (outdic['theta_e'] / outdic['cntp']), extend='both', levels=np.linspace(-2,2,12), cmap='PuOr_r')
    plt.colorbar(label=r'K')
    contours = plt.contour(ranges, np.arange(401), (outdic['v925'] / outdic['cntp']), extend='both', cmap='RdBu', levels=np.linspace(-0.7,0.7,9))
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    ax.set_yticklabels(np.array((np.linspace(0, 200 * 2, 9) - 200) * 3, dtype=int))
    #ax.set_xticklabels(outticks)

    plt.title(r'Shading:$\Delta \theta_{e} anomaly$, contours: 925hpa v-wind anomaly')
    plt.hlines(200,xmin=ranges[0], xmax=ranges[-1], linestyle='dashed')
    plt.xlabel('Hour of day')
    plt.ylabel('North-South distance from core (km)')

    ax = f.add_subplot(326)
    plt.contourf(ranges,np.arange(401), (outdic['div'])/ outdic['cntp']*100, extend='both',  cmap='RdBu', levels=np.linspace(-0.5,0.5,10))
    plt.colorbar(label=r'10$^{-2}$ s$^{-1}$')
    contours = plt.contour(ranges, np.arange(401),(outdic['v925_orig'] / outdic['cntp']), extend='both',levels=np.linspace(-3,3,11), cmap='RdBu')
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    ax.set_yticklabels(np.array((np.linspace(0, 200 * 2, 9) - 200) * 3, dtype=int))
    #ax.set_xticklabels(outticks)
    plt.title(r'Shading:Divergence, contours: 925hpa v-wind')
    plt.hlines(200,xmin=ranges[0], xmax=ranges[-1], linestyle='dashed')
    plt.xlabel('Hour of day')


    plt.tight_layout()
    plt.show()
    #plt.savefig(cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/plots/ERA5_"+str(h).zfill(2)+'_timeseries_short.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png)
    #plt.close()

def plot_timeseries_small():

    x = len(list(range(-38, 4, 3)))
    y = 401

    # outticks = list(range(-30, 1, 5))
    # ranges = np.arange(-30,1,3)
    #
    # outticks = [12,17,22,3,8,13,18]
    #outticks = [6, 11, 16, 21, 2, 7, 12, 17]
    ranges = np.arange(-38, 4, 3)

    h = 17

    outdic = {}
    dummyfile = glob.glob(
        cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/ERA5_composite_cores_*_small_cores.p")
    dummy = pkl.load(open(dummyfile[0], "rb"))

    for k in dummy.keys():
        outdic[k] = np.zeros((y, x))

    for ids, eh in enumerate(range(-38, 4, 3)):

        dic = {}

        def coll(dic, h, eh, year):
            print(h)
            core = pkl.load(open(
                cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/ERA5_composite_cores_LSTA_500w04_15k_" + str(
                    eh) + "UTCERA" + str(h).zfill(2) + '_' + str(year) + "_small_cores.p", "rb"))

            for id, k in enumerate(core.keys()):
                try:
                    dic[k] = dic[k] + core[k]
                except KeyError:
                    dic[k] = core[k]

        for y in range(2006, 2011):
            coll(dic, h, eh, y)

        for k in dic.keys():
            outdic[k][:, ids] = dic[k][:, 190:211].mean(axis=1)
    #
    # pos = np.isclose(outdic['v925'], np.zeros_like(outdic['v925']), atol=0.7)
    # #ipdb.set_trace()
    # ppos = np.where(pos)
    #
    # for pp in ppos:
    #     if pp[0] < 150:
    #         outdic['v925'][pp] = np.nan


    # plt.figure()
    # plt.pcolormesh(outdic['v925'], vmin=0)
    # plt.colorbar()

    outdic['v925'][50:90,5:9] = np.nan


    print(outdic.keys())
    f = plt.figure(figsize=(6, 8), dpi=200)


    ax = f.add_subplot(211)
    plt.contourf(ranges, np.arange(401), (outdic['t']) / outdic['cntp'], extend='both', cmap='RdBu_r',
                levels=[-0.8, -0.7, -0.6, -0.5, -0.4, -0.2, -0.1, 0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8])
    plt.colorbar(label=r'K')
    plt.contour(ranges, np.arange(401), (outdic['t']) / outdic['cntp'], extend='both', colors='k', linewidths=0.1,
                 levels=[-0.8, -0.7, -0.6, -0.5, -0.4, -0.2, -0.1, 0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8])

    contours = plt.contour(ranges, np.arange(401), (outdic['u650'] / outdic['cntp']), extend='both', colors='k',
                           levels=np.arange(-1,0,0.2), linewidths=1)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')


    contours = plt.contour(ranges, np.arange(401),(outdic['v925']/ outdic['cntp']), extend='both',colors='k', linewidths=5, levels=[-50,0,50])  #np.arange(-15,-10,0.5)
    contours = plt.contour(ranges, np.arange(401),(outdic['v925_orig']/ outdic['cntp']), extend='both',colors='k', linewidths=5, levels=[-50,0,50])

    contours = plt.contour(ranges, np.arange(401),(outdic['v925']/ outdic['cntp']), extend='both',colors='k', linewidths=3, levels=[-50,0,50])  #np.arange(-15,-10,0.5)
    contours = plt.contour(ranges, np.arange(401),(outdic['v925_orig']/ outdic['cntp']), extend='both',colors='r', linewidths=3, levels=[-50,0,50])

    contours = plt.contour(ranges, np.arange(401),(outdic['v650_orig']/ outdic['cntp']), extend='both',colors=['r','r','white','b', 'b'], linewidths=0.5, levels=[-2,-1,0,0.25,0.5])  #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')

    plt.axvline(x=-5, color='slategrey')
    plt.axvline(x=-29, color='slategrey')
    plt.axhline(y=200, color='slategrey')
    plt.plot(-5, 200, 'ko')
    plt.plot(0, 200, 'ro')

    plt.text(0.02,0.1, 'ITD 0-line', color='red', fontsize=14, transform=ax.transAxes)
    plt.text(0.02, 0.03, 'v-wind anomaly 0-line', color='k', fontsize=14, transform=ax.transAxes)


    ax.set_yticklabels(np.array((np.linspace(0, 200 * 2, 9) - 200) * 3, dtype=int))
    # ax.set_xticklabels(outticks)

    plt.title('Shading:t-anomaly, contours: 650hpa u-wind anomaly')
    #plt.hlines(200, xmin=ranges[0], xmax=ranges[-1], linewidth=1)
    plt.xlabel('Hour relative to t0 [1700UTC]')
    plt.ylabel('North-South distance from core (km)')

    ax = f.add_subplot(212)


    plt.contourf(ranges, np.arange(401), (outdic['q']) * 1000 / outdic['cntp'], levels=[-0.8, -0.7, -0.6, -0.5, -0.4, -0.2, -0.1, 0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8],
                 cmap='RdBu', extend='both')
    plt.colorbar(label=r'g kg$^{-1}$')

    plt.contour(ranges, np.arange(401), (outdic['q']) * 1000 / outdic['cntp'], levels=[-0.8, -0.7, -0.6, -0.5, -0.4, -0.2, -0.1, 0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8],
     colors = 'k', linewidths = 0.1)

    contours = plt.contour(ranges, np.arange(401), outdic['theta_e'] / outdic['cntp'],colors='white', levels=np.arange(1,3,0.5), linewidths=1.5)

    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')


    contours = plt.contour(ranges, np.arange(401), (outdic['v925'] / outdic['cntp']), extend='both', colors='k',
                           linewidths=5, levels=[-50, 0, 50])  # np.arange(-15,-10,0.5)
    contours = plt.contour(ranges, np.arange(401), (outdic['v925_orig'] / outdic['cntp']), extend='both', colors='k',
                           linewidths=5, levels=[-50, 0, 50])

    contours = plt.contour(ranges, np.arange(401), (outdic['v925'] / outdic['cntp']), extend='both', colors='k',
                           linewidths=3, levels=[-50, 0, 50])  # np.arange(-15,-10,0.5)
    contours = plt.contour(ranges, np.arange(401), (outdic['v925_orig'] / outdic['cntp']), extend='both', colors='r',
                           linewidths=3, levels=[-50, 0, 50])

    plt.axvline(x=-5, color='slategrey')
    plt.axvline(x=-29, color='slategrey')
    plt.axhline(y=200, color='slategrey')
    plt.plot(-5, 200, 'ko')
    plt.plot(0, 200, 'ro')

    ax.set_yticklabels(np.array((np.linspace(0, 200 * 2, 9) - 200) * 3, dtype=int))
    # ax.set_xticklabels(outticks)

    plt.title(r'Shading:q-anomaly, contours: $\Delta \theta_{e}$ anomaly')
    # plt.hlines(200, xmin=ranges[0], xmax=ranges[-1], linewidth=1)
    plt.xlabel('Hour relative to t0 [1700UTC]')
    plt.ylabel('North-South distance from core (km)')

    plt.tight_layout()
    #plt.show()
    plt.savefig(cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/plots/ERA5_"+str(h).zfill(2)+'_timeseries_SMALL_DRY.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png)
    plt.close()



def plot_doug_all_diff(h, eh):

    dic = {}
    dic2 = {}

    name = "ERA5_composite_cores_LSTA_500w04_15k_LSTAcp90"#"ERA5_composite_cores_AMSRE_w1_15k_minusMean"
    name2 = "ERA5_composite_cores_LSTA_500w04_15k_LSTAcp10"

    def coll(dic, h, eh, year,name):
        print(h)
        core = pkl.load(open(
            cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/"+name+str(eh) + "UTCERA"+str(h).zfill(2)+'_'+str(year)+"_small_cores.p", "rb"))
        for id, k in enumerate(core.keys()):
            try:
                dic[k] = dic[k] + core[k]
            except KeyError:
                dic[k] = core[k]


    for y in range(2006,2011):
        coll(dic, h, eh, y,name)

    for y in range(2006,2011):
        coll(dic2, h, eh, y,name2)

    extent = (dic['lsta'].shape[1]-1)/2
    xlen = dic['lsta'].shape[1]
    ylen = dic['lsta'].shape[0]

    xv, yv = np.meshgrid(np.arange(ylen), np.arange(xlen))
    st=30
    xquiv = xv[4::st, 4::st]
    yquiv = yv[4::st, 4::st]

    u = (dic['u925']/ dic['cntp'])[4::st, 4::st]
    v = (dic['v925']/ dic['cntp'])[4::st, 4::st]

    u600 = (dic['u650_orig']/ dic['cntp'])[4::st, 4::st]
    v600 = (dic['v650_orig']/ dic['cntp'])[4::st, 4::st]

    u_orig = (dic['v925_orig']/ dic['cntp'])[4::st, 4::st]
    v_orig = (dic['v925_orig']/ dic['cntp'])[4::st, 4::st]


    f = plt.figure(figsize=(15,7))
    ax = f.add_subplot(231)


    plt.contourf((dic['plsta'] / dic['plcnt'])-((dic2['plsta'] / dic2['plcnt'])), cmap='RdBu_r', extend='both', levels=np.linspace(-5,5,12)) #-(rkernel2_sum / rcnt_sum)
    plt.plot(extent, extent, 'bo')
    plt.colorbar(label='K', format='%1.2f')
    plt.text(0.02,0.08, 'ITD 0-line', color='turquoise', fontsize=12, transform=ax.transAxes)
    plt.text(0.02, 0.03, 'ITD anomaly 0-line', color='k', fontsize=12, transform=ax.transAxes)

    # plt.annotate('ITD 0-line', xy=(0.04, 0.1), xytext=(0, 4), size=15, color='turquoise', xycoords=('figure fraction', 'figure fraction'))
    #              #             textcoords='offset points')   #transform=ax.transAxes
    #pdb.set_trace()

    contours = plt.contour((dic['t'] / dic['cntp'])-(dic2['t'] / dic2['cntp']), extend='both', cmap='PuOr_r', linewidths=2, levels=np.linspace(-1,1,8)) # #, levels=np.arange(1,5, 0.5)
    plt.clabel(contours, inline=True, fontsize=11, fmt='%1.1f')
    #contours2 = plt.contour((dic['v925']) / dic['cntp'], extend='both', cmap='RdBu', levels=np.linspace(-1, 1,9))  # , levels=np.linspace(-1,1,10)#(dic['probc']/ dic['cntc'])*100, extend='both', levels=np.arange(15,70,12), cmap='jet') # #, levels=np.arange(1,5, 0.5)
    #plt.clabel(contours2, inline=True, fontsize=11, fmt='%1.0f')

    contours = plt.contour((dic['v925']/ dic['cntc'])-(dic2['v925']/ dic2['cntc']), extend='both',colors='k', linewidths=4, levels=[-50,0,50])  #np.arange(-15,-10,0.5)
    #plt.clabel(contours, inline=True, fontsize=9, fmt='%1.1f')

    contours = plt.contour((dic['v925_orig']/ dic['cntc'])-(dic2['v925_orig']/ dic2['cntc']), extend='both',colors='turquoise', linewidths=4, levels=[-50,0,50])  #np.arange(-15,-10,0.5)
    #plt.clabel(contours, inline=True, fontsize=9, fmt='%1.1f')

    plt.plot(extent, extent, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6 , dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.title(str(h).zfill(2)+'00UTC | '+str(np.max(dic['cnt']))+' cores, LSTA day-1, ERA5$_{noon}$ 925hPa T anomaly', fontsize=9)


    ax1 = f.add_subplot(232)
    plt.contourf(((dic['lsta'])/ dic['cnt'])-((dic2['lsta'])/ dic2['cnt']), extend='both',  cmap='RdBu_r', levels=np.linspace(-5,5,12)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='K', format='%1.2f')
    contours = plt.contour((dic['probc']/ dic['cntc'])*100 - (dic2['probc']/ dic2['cntc'])*100, extend='both', cmap='jet', linewidths=2)  #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.1f')



    #qu = ax1.quiver(xquiv, yquiv, u, v, scale=50)
    plt.plot(extent, extent, 'bo')
    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('LSTA day0, Contours: CMORPH rainP>5mm [6am|day-1 to 10am|day0]', fontsize=9)

    ax1 = f.add_subplot(233)
    plt.contourf(((dic['q'])*1000/ dic['cntp'])-((dic2['q'])*1000/ dic2['cntp']), extend='both',  cmap='RdBu', levels=np.linspace(-0.8,0.8,10)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='g kg-1')
    contours = plt.contour((dic['shear'] / dic['cntp'])-(dic2['shear'] / dic2['cntp']), extend='both', cmap='viridis_r', linewidths=2) #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    plt.plot(extent, extent, 'bo')
    #qu = ax1.quiver(xquiv, yquiv, u, v, scale=50)
    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('925hPa q anomaly, contours: 650hPa-925hPa zonal shear', fontsize=9)

    ax1 = f.add_subplot(234)
 #   plt.contourf(((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r', vmin=-1.5, vmax=1.5) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.contourf((((dic['div'])/ dic['cntp'])-((dic2['div'])/ dic2['cntp']))*100, extend='both',  cmap='PuOr', levels=np.linspace(-0.7, 0.7, 10)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='10$^{-2}$ s$^{-1}$', format='%1.3f')
    plt.plot(extent, extent, 'bo')
    # contours = plt.contour((dic['v650_orig'] / dic['cntp']), extend='both', cmap='viridis') #np.arange(-15,-10,0.5)
    # plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    # qu = ax1.quiver(xquiv, yquiv, u, v, scale=15)
    #qk = plt.quiverkey(qu, 0.2, 0.02,1, '1 m s$^{-1}$',
    #                   labelpos='E', coordinates='figure')

    ax1.streamplot(xv, yv, (dic['u925']/ dic['cntp'])-(dic2['u925']/ dic2['cntp']), (dic['v925']/ dic['cntp'])-(dic2['v925']/ dic2['cntp']), density=[0.5, 1])

    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('925hPa divergence, vectors: 925hPa wind anomaly', fontsize=9)


    ax1 = f.add_subplot(235)
    #   plt.contourf(((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r', vmin=-1.5, vmax=1.5) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
     # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.contourf((dic['tciwmid'] / dic['cntp'])-(dic2['tciwmid'] / dic2['cntp']), extend='both', cmap='PuOr', levels=np.linspace(-0.02, 0.02, 10))
    plt.colorbar(label=r'Pa s$^{-1}$', format='%1.3f')
    plt.plot(extent, extent, 'bo')
    # contours = plt.contour(((dic['v925_orig']) / dic['cntp']) , extend='both', cmap='RdBu', levels=np.linspace(-2,2,11), linewidths=2) #np.arange(-15,-10,0.5)
    # plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    ax1.streamplot(xv, yv, dic['u925_orig']/ dic['cntp']-dic2['u925_orig']/ dic2['cntp'], dic['v925_orig']/ dic['cntp']-dic2['v925_orig']/ dic2['cntp'], density=[0.5, 1])

    # qu = ax1.quiver(xquiv, yquiv, u_orig, v_orig, scale=30)
    # qk = plt.quiverkey(qu, 0.55, 0.02, 1, '1 m s$^{-1}$',
    #                    labelpos='E', coordinates='figure')

    ax1.set_xticklabels(np.array((np.linspace(0, extent * 2, 9) - 100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent * 2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title(r'500hpa omega, vectors: 925hPa wind', fontsize=9)

    ax1 = f.add_subplot(236)
    #   plt.contourf(((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r', vmin=-1.5, vmax=1.5) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.contourf(((dic['theta_e']) / dic['cntp'])-((dic2['theta_e']) / dic2['cntp']), extend='both', cmap='RdBu', levels=np.linspace(-2.5,2.5,14)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label=r'K')
    plt.plot(extent, extent, 'bo')
    contours = plt.contour((dic['u650'] / dic['cntp'])-(dic2['u650'] / dic2['cntp']), extend='both', cmap='PuOr', levels=[-2,-1.5,-1,-0.5,-0.2,-0.1, 0,0.1,0.2,0.5,1,1.5,2], linewidths=2) #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    # qu = ax1.quiver(xquiv, yquiv, u600, v600, scale=20)
    # qk = plt.quiverkey(qu, 0.9, 0.02, 1, '1 m s$^{-1}$',
    #                    labelpos='E', coordinates='figure')

    ax1.set_xticklabels(np.array((np.linspace(0, extent * 2, 9) - 100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent * 2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title(r'$\Delta \theta_{e}$ anomaly, contours: 650hPa u-wind anomaly', fontsize=9)

    plt.tight_layout()
    plt.show()
    plt.savefig(cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/plots/"+name+str(h).zfill(2)+'_'+str(eh).zfill(2)+'_DIFF_9010.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png)
    #plt.close()



def plot_doug_all(h, eh):

    dic = {}
    dic2 = {}

    name = "ERA5_composite_cores_LSTA_500w04_15k_"#"ERA5_composite_cores_AMSRE_w1_15k_minusMean"


    def coll(dic, h, eh, year):
        print(h)
        core = pkl.load(open(
            cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/"+name+str(eh) + "UTCERA"+str(h).zfill(2)+'_'+str(year)+"_small_cores.p", "rb"))
        for id, k in enumerate(core.keys()):
            try:
                dic[k] = dic[k] + core[k]
            except KeyError:
                dic[k] = core[k]



    if abs(eh) <= h:
        time = h+eh
    else:
        time = (h+24)+eh


    for y in range(2006,2011):
        coll(dic, h, eh, y)

    extent = (dic['lsta'].shape[1]-1)/2
    xlen = dic['lsta'].shape[1]
    ylen = dic['lsta'].shape[0]

    xv, yv = np.meshgrid(np.arange(ylen), np.arange(xlen))
    st=30
    xquiv = xv[4::st, 4::st]
    yquiv = yv[4::st, 4::st]

    u = (dic['u925']/ dic['cntp'])[4::st, 4::st]
    v = (dic['v925']/ dic['cntp'])[4::st, 4::st]

    u600 = (dic['u650_orig']/ dic['cntp'])[4::st, 4::st]
    v600 = (dic['v650_orig']/ dic['cntp'])[4::st, 4::st]

    u_orig = (dic['v925_orig']/ dic['cntp'])[4::st, 4::st]
    v_orig = (dic['v925_orig']/ dic['cntp'])[4::st, 4::st]


    f = plt.figure(figsize=(10,7))
    ax = f.add_subplot(221)


    plt.contourf((dic['plsta'] / dic['plcnt'])-np.mean((dic['plsta'] / dic['plcnt'])), cmap='RdBu_r', levels=np.linspace(-1.5,1.5,16), extend='both') #-(rkernel2_sum / rcnt_sum)
    plt.plot(extent, extent, 'bo')
    plt.colorbar(label='K', format='%1.2f')
    plt.text(0.02,0.08, 'ITD 0-line', color='turquoise', fontsize=12, transform=ax.transAxes)
    plt.text(0.02, 0.03, 'V650 0-line', color='k', fontsize=12, transform=ax.transAxes)

    # plt.annotate('ITD 0-line', xy=(0.04, 0.1), xytext=(0, 4), size=15, color='turquoise', xycoords=('figure fraction', 'figure fraction'))
    #              #             textcoords='offset points')   #transform=ax.transAxes
    #pdb.set_trace()

    contours = plt.contour((dic['t'] / dic['cntp']), extend='both',levels=[ -0.8,-0.7, -0.6,-0.5,-0.4,-0.2,0, 0.2,0.4,0.5,0.6, 0.7, 0.8], cmap='PuOr_r', linewidths=2) # #, levels=np.arange(1,5, 0.5)
    plt.clabel(contours, inline=True, fontsize=11, fmt='%1.1f')
    #contours2 = plt.contour((dic['v925']) / dic['cntp'], extend='both', cmap='RdBu', levels=np.linspace(-1, 1,9))  # , levels=np.linspace(-1,1,10)#(dic['probc']/ dic['cntc'])*100, extend='both', levels=np.arange(15,70,12), cmap='jet') # #, levels=np.arange(1,5, 0.5)
    #plt.clabel(contours2, inline=True, fontsize=11, fmt='%1.0f')
    qu = ax.quiver(xquiv, yquiv, u600 * 0, v600, scale=15)
    contours = plt.contour((dic['v650_orig']/ dic['cntc']), extend='both',colors='k', linewidths=4, levels=[-50,0,50])  #np.arange(-15,-10,0.5)
    #plt.clabel(contours, inline=True, fontsize=9, fmt='%1.1f')

    contours = plt.contour((dic['v925_orig']/ dic['cntc']), extend='both',colors='turquoise', linewidths=4, levels=[-50,0,50])  #np.arange(-15,-10,0.5)
    #plt.clabel(contours, inline=True, fontsize=9, fmt='%1.1f')

    plt.plot(extent, extent, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6 , dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.title(str(h).zfill(2)+'00UTC | '+str(np.max(dic['cnt']))+' cores, LSTA day-1, ERA5$_{noon}$ 925hPa T anomaly', fontsize=9)


    ax1 = f.add_subplot(222)
    plt.contourf(((dic['lsta'])/ dic['cnt'])-np.mean((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r', levels=np.linspace(-1.5,1.5,16)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='K', format='%1.2f')
    contours = plt.contour((dic['probc']/ dic['cntc'])*100, extend='both', levels=np.arange(15,70,12), cmap='jet', linewidths=2)  #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.1f')



    #qu = ax1.quiver(xquiv, yquiv, u, v, scale=50)
    plt.plot(extent, extent, 'bo')
    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('LSTA day0, Contours: CMORPH rainP>5mm [6am|day-1 to 10am|day0]', fontsize=9)

    ax1 = f.add_subplot(223)
    plt.contourf(((dic['q'])*1000/ dic['cntp']), extend='both',  cmap='RdBu',levels=np.linspace(-0.9,0.9,16)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='g kg-1')
    contours = plt.contour((dic['shear'] / dic['cntp']), extend='both',levels=np.arange(-17,-12,0.5), cmap='viridis_r', linewidths=2) #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    plt.plot(extent, extent, 'bo')
    #qu = ax1.quiver(xquiv, yquiv, u, v, scale=50)
    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('925hPa q anomaly, contours: 650hPa-925hPa zonal shear', fontsize=9)

    ax1 = f.add_subplot(224)
 #   plt.contourf(((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r', vmin=-1.5, vmax=1.5) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.contourf(((dic['div'])/ dic['cntp'])*100, extend='both',  cmap='PuOr', levels=np.linspace(-0.7,0.7,10)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='10$^{-2}$ s$^{-1}$', format='%1.3f')
    plt.plot(extent, extent, 'bo')
    # contours = plt.contour((dic['v650_orig'] / dic['cntp']), extend='both', cmap='viridis') #np.arange(-15,-10,0.5)
    # plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    # qu = ax1.quiver(xquiv, yquiv, u, v, scale=15)
    # qk = plt.quiverkey(qu, 0.2, 0.02,1, '1 m s$^{-1}$',
    #                   labelpos='E', coordinates='figure')

    ax1.streamplot(xv, yv, (dic['u925']/ dic['cntp']), (dic['v925']/ dic['cntp']), density=[0.5, 1])

    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('925hPa divergence, vectors: 925hPa wind anomaly', fontsize=9)


    plt.tight_layout()
    plt.show()
    plt.savefig(cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/plots/"+name+str(h).zfill(2)+'_'+str(eh).zfill(2)+'_'+str(time).zfill(2)+'UTC.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png)
    #plt.close()


def plot_doug_CLIM(h, eh):

    dic = {}
    dic2 = {}

    name = "ERA5_composite_cores_LSTA_500w04_15k_-80_clim"#"ERA5_composite_cores_AMSRE_w1_15k_minusMean"


    def coll(dic, h, eh, year):
        print(h)
        core = pkl.load(open(
            cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/"+name+str(eh) + "UTCERA"+str(h).zfill(2)+'_'+str(year)+"_small_cores.p", "rb"))
        for id, k in enumerate(core.keys()):
            try:
                dic[k] = dic[k] + core[k]
            except KeyError:
                dic[k] = core[k]


    for y in range(2006,2011):
        coll(dic, h, eh, y)

    extent = (dic['lsta'].shape[1]-1)/2
    xlen = dic['lsta'].shape[1]
    ylen = dic['lsta'].shape[0]

    xv, yv = np.meshgrid(np.arange(ylen), np.arange(xlen))
    st=30
    xquiv = xv[4::st, 4::st]
    yquiv = yv[4::st, 4::st]

    u = (dic['u925']/ dic['cntp'])[4::st, 4::st]
    v = (dic['v925']/ dic['cntp'])[4::st, 4::st]

    u600 = (dic['u650_orig']/ dic['cntp'])[4::st, 4::st]
    v600 = (dic['v650_orig']/ dic['cntp'])[4::st, 4::st]

    u_orig = (dic['v925_orig']/ dic['cntp'])[4::st, 4::st]
    v_orig = (dic['v925_orig']/ dic['cntp'])[4::st, 4::st]


    f = plt.figure(figsize=(15,7))
    ax = f.add_subplot(231)


    plt.contourf((dic['plsta'] / dic['plcnt'])-np.mean((dic['plsta'] / dic['plcnt'])), cmap='RdBu_r', levels=np.linspace(-1.5,1.5,12), extend='both') #-(rkernel2_sum / rcnt_sum)
    plt.plot(extent, extent, 'bo')
    plt.colorbar(label='K', format='%1.2f')
    plt.text(0.02,0.08, 'ITD 0-line', color='turquoise', fontsize=12, transform=ax.transAxes)
    plt.text(0.02, 0.03, 'ITD anomaly 0-line', color='k', fontsize=12, transform=ax.transAxes)

    # plt.annotate('ITD 0-line', xy=(0.04, 0.1), xytext=(0, 4), size=15, color='turquoise', xycoords=('figure fraction', 'figure fraction'))
    #              #             textcoords='offset points')   #transform=ax.transAxes
    #pdb.set_trace()

    contours = plt.contour((dic['t'] / dic['cntp']), extend='both',levels=[ -0.8,-0.7, -0.6,-0.5,-0.4,-0.2,0, 0.2,0.4,0.5,0.6, 0.7, 0.8], cmap='PuOr_r', linewidths=2) # #, levels=np.arange(1,5, 0.5)
    plt.clabel(contours, inline=True, fontsize=11, fmt='%1.1f')
    #contours2 = plt.contour((dic['v925']) / dic['cntp'], extend='both', cmap='RdBu', levels=np.linspace(-1, 1,9))  # , levels=np.linspace(-1,1,10)#(dic['probc']/ dic['cntc'])*100, extend='both', levels=np.arange(15,70,12), cmap='jet') # #, levels=np.arange(1,5, 0.5)
    #plt.clabel(contours2, inline=True, fontsize=11, fmt='%1.0f')

    contours = plt.contour((dic['v925']/ dic['cntc']), extend='both',colors='k', linewidths=4, levels=[-50,0,50])  #np.arange(-15,-10,0.5)
    #plt.clabel(contours, inline=True, fontsize=9, fmt='%1.1f')

    contours = plt.contour((dic['v925_orig']/ dic['cntc']), extend='both',colors='turquoise', linewidths=4, levels=[-50,0,50])  #np.arange(-15,-10,0.5)
    #plt.clabel(contours, inline=True, fontsize=9, fmt='%1.1f')

    plt.plot(extent, extent, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6 , dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.title(str(h).zfill(2)+'00UTC | '+str(np.max(dic['cnt']))+' cores, LSTA day-1, ERA5$_{noon}$ 925hPa T anomaly', fontsize=9)


    ax1 = f.add_subplot(232)
    plt.contourf(((dic['lsta'])/ dic['cnt'])-np.mean((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r', levels=np.linspace(-1.5,1.5,12)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='K', format='%1.2f')
    contours = plt.contour((dic['probc']/ dic['cntc'])*100, extend='both', levels=np.arange(15,70,12), cmap='jet', linewidths=2)  #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.1f')



    #qu = ax1.quiver(xquiv, yquiv, u, v, scale=50)
    plt.plot(extent, extent, 'bo')
    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('LSTA day0, Contours: CMORPH rainP>5mm [6am|day-1 to 10am|day0]', fontsize=9)

    ax1 = f.add_subplot(233)
    plt.contourf(((dic['q'])*1000/ dic['cntp']), extend='both',  cmap='jet', levels=np.arange(11,17,0.5)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='g kg-1')
    contours = plt.contour((dic['shear'] / dic['cntp']), extend='both',levels=np.arange(-17,-12,0.5), cmap='viridis_r', linewidths=2) #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    plt.plot(extent, extent, 'bo')
    #qu = ax1.quiver(xquiv, yquiv, u, v, scale=50)
    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('925hPa q anomaly, contours: 650hPa-925hPa zonal shear', fontsize=9)

    ax1 = f.add_subplot(234)
 #   plt.contourf(((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r', vmin=-1.5, vmax=1.5) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.contourf(((dic['div'])/ dic['cntp'])*100, extend='both',  cmap='PuOr', levels=np.linspace(-0.7,0.7,10)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='10$^{-2}$ s$^{-1}$', format='%1.3f')
    plt.plot(extent, extent, 'bo')
    # contours = plt.contour((dic['v650_orig'] / dic['cntp']), extend='both', cmap='viridis') #np.arange(-15,-10,0.5)
    # plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    # qu = ax1.quiver(xquiv, yquiv, u, v, scale=15)
    #qk = plt.quiverkey(qu, 0.2, 0.02,1, '1 m s$^{-1}$',
    #                   labelpos='E', coordinates='figure')

    ax1.streamplot(xv, yv, (dic['u925']/ dic['cntp']), (dic['v925']/ dic['cntp']), density=[0.5, 1])

    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('925hPa divergence, vectors: 925hPa wind anomaly', fontsize=9)


    ax1 = f.add_subplot(235)
    #   plt.contourf(((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r', vmin=-1.5, vmax=1.5) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
     # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.contourf((dic['tciwmid'] / dic['cntp']), extend='both', cmap='PuOr', levels=np.linspace(-0.02, 0.02, 10))
    plt.colorbar(label=r'Pa s$^{-1}$', format='%1.3f')
    plt.plot(extent, extent, 'bo')
    # contours = plt.contour(((dic['v925_orig']) / dic['cntp']) , extend='both', cmap='RdBu', levels=np.linspace(-2,2,11), linewidths=2) #np.arange(-15,-10,0.5)
    # plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    ax1.streamplot(xv, yv, dic['u925_orig']/ dic['cntp'], dic['v925_orig']/ dic['cntp'], density=[0.5, 1])

    # qu = ax1.quiver(xquiv, yquiv, u_orig, v_orig, scale=30)
    # qk = plt.quiverkey(qu, 0.55, 0.02, 1, '1 m s$^{-1}$',
    #                    labelpos='E', coordinates='figure')

    ax1.set_xticklabels(np.array((np.linspace(0, extent * 2, 9) - 100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent * 2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title(r'500hpa omega, vectors: 925hPa wind', fontsize=9)

    ax1 = f.add_subplot(236)
    #   plt.contourf(((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r', vmin=-1.5, vmax=1.5) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.contourf(((dic['theta_e']) / dic['cntp']), extend='both', cmap='RdBu', levels=np.arange(5,13,0.5)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label=r'K')
    plt.plot(extent, extent, 'bo')
    contours = plt.contour((dic['u650'] / dic['cntp']), extend='both', cmap='PuOr', levels=[-2,-1.5,-1,-0.5,-0.2,-0.1, 0,0.1,0.2,0.5,1,1.5,2], linewidths=2) #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    # qu = ax1.quiver(xquiv, yquiv, u600, v600, scale=20)
    # qk = plt.quiverkey(qu, 0.9, 0.02, 1, '1 m s$^{-1}$',
    #                    labelpos='E', coordinates='figure')

    ax1.set_xticklabels(np.array((np.linspace(0, extent * 2, 9) - 100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent * 2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title(r'$\Delta \theta_{e}$ anomaly, contours: 650hPa u-wind anomaly', fontsize=9)

    plt.tight_layout()
    plt.show()
    plt.savefig(cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/plots/"+name+str(h).zfill(2)+'_'+str(eh).zfill(2)+'_CLIM.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png)
    #plt.close()


def plot_doug_paper(h, eh):

    dic = {}
    name = "ERA5_composite_cores_LSTA_500w04_15k_"


    def coll(dic, h, eh, year):
        print(h)
        core = pkl.load(open(
            cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/"+name+str(eh) + "UTCERA"+str(h).zfill(2)+'_'+str(year)+"_small_cores.p", "rb"))
        for id, k in enumerate(core.keys()):
            try:
                dic[k] = dic[k] + core[k]
            except KeyError:
                dic[k] = core[k]


    for y in range(2006,2011):
        coll(dic, h, eh, y)

    extent = (dic['lsta'].shape[1]-1)/2
    xlen = dic['lsta'].shape[1]
    ylen = dic['lsta'].shape[0]

    xv, yv = np.meshgrid(np.arange(ylen), np.arange(xlen))
    st=30
    xquiv = xv[4::st, 4::st]
    yquiv = yv[4::st, 4::st]

    u = (dic['u925']/ dic['cntp'])[4::st, 4::st]
    v = (dic['v925']/ dic['cntp'])[4::st, 4::st]

    u600 = (dic['u650_orig']/ dic['cntp'])[4::st, 4::st]
    v600 = (dic['v650_orig']/ dic['cntp'])[4::st, 4::st]

    u_orig = (dic['v925_orig']/ dic['cntp'])[4::st, 4::st]
    v_orig = (dic['v925_orig']/ dic['cntp'])[4::st, 4::st]


    f = plt.figure(figsize=(15,7))


    f = plt.figure(figsize=(10,8))
    ax = f.add_subplot(221)

    plt.contourf((dic['plsta'] / dic['plcnt']), cmap='RdBu_r', levels=[-1.5,-1.25,-1 ,-0.75,-0.5,-0.25,0.25,0.5,0.75,1,1.25, 1.5], extend='both') #-(rkernel2_sum / rcnt_sum)
    plt.plot(extent, extent, 'bo')
    plt.colorbar(label='K')
    #pdb.set_trace()

    contours = plt.contour((dic['probc']/ dic['cntc'])*100, extend='both', levels=np.arange(15,70,12), cmap='jet') # #, levels=np.arange(1,5, 0.5)
    plt.clabel(contours, inline=True, fontsize=11, fmt='%1.0f')
    plt.plot(extent, extent, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6 , dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.title('1800UTC | '+str(np.max(dic['cnt']))+' cores, LSTA', fontsize=9)


    ax1 = f.add_subplot(222)
    plt.contourf(((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r',levels=[-1.5,-1.25,-1 ,-0.75,-0.5,-0.25,0.25,0.5,0.75,1,1.25, 1.5]) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='K')
    contours = plt.contour((dic['t'] / dic['cntp']), extend='both',levels=[ -0.8,-0.7, -0.6,-0.5,-0.4,-0.2,-0.1, 0.1,0.2,0.4,0.5,0.6, 0.7, 0.8], cmap='PuOr_r') #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.1f')
    #qu = ax1.quiver(xquiv, yquiv, u, v, scale=50)
    plt.plot(extent, extent, 'bo')
    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('Shading: MSG LSTA, Contours: ERA5 925hPa T', fontsize=9)

    ax1 = f.add_subplot(223)
    plt.contourf(((dic['q'])*1000/ dic['cntp']), extend='both',  cmap='RdBu',levels=np.arange(-1,1.1,0.1)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='g kg-1')
    contours = plt.contour((dic['shear'] / dic['cntp']), extend='both',levels=np.arange(-17,-12,0.25), cmap='viridis') #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    plt.plot(extent, extent, 'bo')
    #qu = ax1.quiver(xquiv, yquiv, u, v, scale=50)
    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('Shading: 925hPa q anomaly, Contours: 650hPa-925hPa wind shear ', fontsize=9)

    ax1 = f.add_subplot(224)
 #   plt.contourf(((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r', vmin=-1.5, vmax=1.5) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.contourf(((dic['div'])/ dic['cntp'])*100, extend='both',  cmap='RdBu', levels=np.linspace(-0.5,0.5,10)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='K')
    plt.plot(extent, extent, 'bo')
    # contours = plt.contour((dic['shear'] / dic['cntp']), extend='both',levels=np.arange(-17,-12,0.5), cmap='viridis') #np.arange(-15,-10,0.5)
    # plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    qu = ax1.quiver(xquiv, yquiv, u, v, scale=15)
    qk = plt.quiverkey(qu, 0.9, 0.02,1, '1 m s$^{-1}$',
                       labelpos='E', coordinates='figure')

    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('Shading: Divergence, vectors: 925hPa wind anomaly', fontsize=9)


    plt.tight_layout()
    plt.show()
    plt.savefig(cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/plots/ERA5_"+str(h).zfill(2)+'_'+str(eh).zfill(2)+'.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png)
    plt.close()

def plot_doug_small(h, eh):

    dic = {}
    dic2 = {}

    name = "ERA5_composite_cores_LSTA_500w04_15k_"#"ERA5_composite_cores_AMSRE_w1_15k_minusMean"


    def coll(dic, h, eh, year):
        print(h)
        core = pkl.load(open(
            cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/"+name+str(eh) + "UTCERA"+str(h).zfill(2)+'_'+str(year)+"_small_cores.p", "rb"))
        for id, k in enumerate(core.keys()):
            try:
                dic[k] = dic[k] + core[k]
            except KeyError:
                dic[k] = core[k]


    for y in range(2006,2011):
        coll(dic, h, eh, y)

    extent = (dic['lsta'].shape[1]-1)/2
    xlen = dic['lsta'].shape[1]
    ylen = dic['lsta'].shape[0]

    xv, yv = np.meshgrid(np.arange(ylen), np.arange(xlen))
    st=30
    xquiv = xv[4::st, 4::st]
    yquiv = yv[4::st, 4::st]

    u = (dic['u925']/ dic['cntp'])[4::st, 4::st]
    v = (dic['v925']/ dic['cntp'])[4::st, 4::st]

    u600 = (dic['u650_orig']/ dic['cntp'])[4::st, 4::st]
    v600 = (dic['v650_orig']/ dic['cntp'])[4::st, 4::st]

    u_orig = (dic['v925_orig']/ dic['cntp'])[4::st, 4::st]
    v_orig = (dic['v925_orig']/ dic['cntp'])[4::st, 4::st]

    f = plt.figure(figsize=(10, 4))
    ax = f.add_subplot(121)

    plt.contourf((dic['lsta'] / dic['cnt']), cmap='RdBu_r', levels = np.linspace(-1.5, 1.5, 12), extend = 'both')
    # plt.plot(extent, extent, 'bo')
    plt.colorbar(label='K')
    # pdb.set_trace()

    contours = plt.contour(((dic['div'])/ dic['cntp'])*100, extend='both',  cmap='PuOr_r', levels=np.linspace(-0.7,0.7,10)) # #, levels=np.arange(1,5, 0.5)
    qu = ax.quiver(xquiv, yquiv, u, v, scale=15)
    qk = plt.quiverkey(qu, 0.9, 0.02,1, '1 m s$^{-1}$',
                       labelpos='E', coordinates='figure')


    plt.clabel(contours, inline=True, fontsize=11, fmt='%1.0f')
    plt.plot(extent, extent, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, extent * 2, 9) - 100) * 6, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent * 2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.title('23-01UTC | ' + str(np.max(dic['cnt'])) + ' cores, LSTA & 06-06UTC antecedent rain', fontsize=9)


    ax1 = f.add_subplot(122)
    plt.contourf(((dic['q'])*1000/ dic['cntp']), extend='both',  cmap='RdBu',levels=np.arange(-1,1.1,0.1)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='g kg-1')
    contours = plt.contour((dic['shear'] / dic['cntp']), extend='both',levels=np.arange(-15.5,-12,0.5), cmap='viridis_r') #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    plt.plot(extent, extent, 'bo')
    #qu = ax1.quiver(xquiv, yquiv, u, v, scale=50)
    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('Shading: 950hPa q anomaly, Contours: 600hPa-925hPa wind shear ', fontsize=9)

    plt.tight_layout()
    plt.show()
    #plt.savefig('/users/global/cornkle/figs/LSTA-bullshit/corrected_LSTA/system_scale/doug/large_scale/'+str(hour).zfill(2)+'_JUN.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png')
    # plt.close()
