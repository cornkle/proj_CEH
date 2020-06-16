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
from scipy import ndimage
from utils import u_statistics as ustats
import salem


import pickle as pkl


matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)




def plot_amsr_lsta_paper(h, eh):

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


    dic2 = {}

    name = "ERA5_composite_cores_AMSRE_500w04_15k_minusMean"# "ERA5_cores_WET_SM0GT0.01-1GT0.01_noMeteosatFilter_AMSRE"  # "ERA5_composite_cores_AMSRE_w1_15k_minusMean"

    def coll(dic2, h, eh, year):
        print(h)
        core = pkl.load(open(
            cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/" + name + str(eh) + "UTCERA" + str(h).zfill(
                2) + '_' + str(year) + "_small_cores.p", "rb"))
        for id, k in enumerate(core.keys()):
            try:
                dic2[k] = dic2[k] + core[k]
            except KeyError:
                dic2[k] = core[k]

    for y in range(2006, 2011):
        coll(dic2, h, eh, y)

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


    f = plt.figure(figsize=(10,8))
    ax = f.add_subplot(221)
    levels = list(np.arange(-2.25,0,0.25))+list(np.arange(0.25,2.5,0.25))
    plt.contourf((dic['plsta'] / dic['plcnt'])-0.5, cmap='RdBu_r', levels=levels, extend='both') #-(rkernel2_sum / rcnt_sum)
    plt.plot(extent, extent, 'ro')
    plt.colorbar(label='K')
    #pdb.set_trace()


    contours = plt.contour(((dic2['lsta-2']) / (dic2['cnt-2'])), extend='both', cmap='PuOr', levels=[-1.5,-1,-0.5,0,0.5,1,1.5])
    # plt.clabel(contours, inline=True, fontsize=11, fmt='%1.0f')


    contours = plt.contour((dic['v925']/ dic['cntp']), extend='both',colors='k', linewidths=5, levels=[-50,0.06,50])  #np.arange(-15,-10,0.5)
    contours = plt.contour((dic['v925_orig']/ dic['cntp']), extend='both',colors='k', linewidths=5, levels=[-50,0.06,50])

    contours = plt.contour((dic['v925']/ dic['cntp']), extend='both',colors='k', linewidths=3, levels=[-50,0.06,50])  #np.arange(-15,-10,0.5)
    contours = plt.contour((dic['v925_orig']/ dic['cntp']), extend='both',colors='r', linewidths=3, levels=[-50,0.06,50])
    plt.axvline(x=200, linestyle='dashed', color='k', linewidth=1.2)
    plt.axhline(y=200, linestyle='dashed', color='k', linewidth=1.2)
    plt.plot(extent, extent, 'ro')
    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6 , dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.title('1800UTC | '+str(np.max(dic['cnt']))+' cores, Day-1: LSTA (shading) & SM (contour)', fontsize=9)


    ax1 = f.add_subplot(222)
    plt.contourf(((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r',levels=levels) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='K')

    contours = plt.contour(((dic2['lsta0']) / (dic2['cnt0'])), extend='both', cmap='PuOr',
                     levels=[-1.2,-0.8,-0.4, 0,0.4,0.8,1.2], linewidths=2) #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.1f')
    #qu = ax1.quiver(xquiv, yquiv, u, v, scale=50)
    plt.axvline(x=200, linestyle='dashed', color='k', linewidth=1.2)
    plt.axhline(y=200, linestyle='dashed', color='k', linewidth=1.2)
    plt.plot(extent, extent, 'ro')
    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('Day0: LSTA (shading) & SM (contour)', fontsize=9)

    ax1 = f.add_subplot(223)
    plt.contourf(((dic['q']-dic['qclim'])*1000/ dic['cntp']), extend='both',  cmap='RdBu',levels=np.arange(-1,1.1,0.1)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='g kg-1')
    #contours = plt.contour((dic['u650'] / dic['cntp']), extend='both',levels=np.arange(-16.5,-11.5,1), colors='k', linestyles='solid', linewidths=1) #np.arange(-15,-10,0.5)
    contours = plt.contour((dic['u650'] / dic['cntp']), extend='both', levels=np.arange(-1.5, -0.5, 0.1), colors='k',
                           linestyles='solid', linewidths=1)  # np.arange(-15,-10,0.5)

    ax1.streamplot(xv, yv, dic['u650_orig'] * 0.01, dic['v650_orig'], density=[0.5, 1], linewidth=0.5, color='k')

    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    plt.axvline(x=200, linestyle='dashed', color='k', linewidth=1.2)
    plt.axhline(y=200, linestyle='dashed', color='k', linewidth=1.2)
    plt.plot(extent, extent, 'ro')
    #qu = ax1.quiver(xquiv, yquiv, u, v, scale=50)
    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('Shading: 925hPa q anomaly, Contours: 650hPa wind anomaly ', fontsize=9)

    ax1 = f.add_subplot(224)
 #   plt.contourf(((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r', vmin=-1.5, vmax=1.5) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.contourf(((dic['div'])/ dic['cntp'])*100, extend='both',  cmap='RdBu', levels=np.linspace(-0.7,0.7,12)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='K')
    plt.axvline(x=200, linestyle='dashed', color='k', linewidth=1.2)
    plt.axhline(y=200, linestyle='dashed', color='k', linewidth=1.2)
    plt.plot(extent, extent, 'ro')
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
    plt.savefig(cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/plots/ERA5_NEW_"+str(h).zfill(2)+'_'+str(eh).zfill(2)+'_u650_streamline.eps')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png)
    plt.close()

def plot_amsr_paper(h, eh):

    dic = {}

    name = "ERA5_cores_2hOverlap_AMSRE_SMALL_AMSL_new" #"ERA5_composite_cores_AMSRE_500w04_15k_minusMean"  # "ERA5_composite_cores_AMSRE_w1_15k_minusMean"
    #name = "ERA5_cores_DRY_SM0LT-3-1LT-1.5_noMeteosatFilter_AMSRE"
    #name='ERA5_cores_WET_TAG_noMeteosatFilter_AMSRE'
    def coll(dic2, h, eh, year):
        print(h)
        core = pkl.load(open(
            cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/core_txt/" + name + str(eh) + "UTCERA" + str(h).zfill(
                2) + '_' + str(year) + ".p", "rb"))

        # core = pkl.load(open(
        #     cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/" + name + str(eh) + "UTCERA" + str(h).zfill(
        #         2) + '_' + str(year) + "_small_cores.p", "rb"))

        for id, k in enumerate(core.keys()):
            try:
                dic2[k] = dic2[k] + core[k]
            except KeyError:
                dic2[k] = core[k]

    for y in range(2006, 2011):
        coll(dic, h, eh, y)

    extent = (dic['lsta0'].shape[1]-1)/2
    xlen = dic['lsta0'].shape[1]
    ylen = dic['lsta0'].shape[0]

    xv, yv = np.meshgrid(np.arange(ylen), np.arange(xlen))
    st=30
    xquiv = xv[4::st, 4::st]
    yquiv = yv[4::st, 4::st]

    u = (dic['u925']/ dic['cnte'])[4::st, 4::st]
    v = (dic['v925']/ dic['cnte'])[4::st, 4::st]
    levels = list(np.arange(-1.2,0,0.3)) + list(np.arange(0.3,1.5,0.3))

    f = plt.figure(figsize=(9.5,7), dpi=200)
    ax = f.add_subplot(221)
    plt.contourf(((dic['lsta-2']) / (dic['cnt-2'])), extend='both', cmap='RdBu',
                 levels=levels)  # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='%') #-(rkernel2_sum / rcnt_sum)
    plt.plot(extent, extent, 'ro')

    #pdb.set_trace()


    #contours = plt.contour(((dic['lsta-2']) / (dic['cnt-2'])), extend='both', cmap='PuOr', levels=[-1.5,-1,-0.5,0,0.5,1,1.5])
    # plt.clabel(contours, inline=True, fontsize=11, fmt='%1.0f')


    contours = plt.contour((dic['v925']/ dic['cnte']), extend='both',colors='k', linewidths=5, levels=[-50,0.06,50])  #np.arange(-15,-10,0.5)
    contours = plt.contour((dic['v925_orig']/ dic['cnte']), extend='both',colors='k', linewidths=5, levels=[-50,0.06,50])

    contours = plt.contour((dic['v925']/ dic['cnte']), extend='both',colors='k', linewidths=3, levels=[-50,0.06,50])  #np.arange(-15,-10,0.5)
    contours = plt.contour((dic['v925_orig']/ dic['cnte']), extend='both',colors='r', linewidths=3, levels=[-50,0.06,50])
    plt.axvline(x=200, linestyle='dashed', color='k', linewidth=1.2)
    plt.axhline(y=200, linestyle='dashed', color='k', linewidth=1.2)
    plt.plot(extent, extent, 'ro')
    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6 , dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.title('1700UTC | '+str(np.max(dic['cnte']))+' cores, AMSR-E Day-1', fontsize=9)


    ax1 = f.add_subplot(222)
    plt.contourf(((dic['lsta0']) / (dic['cnt0'])), extend='both', cmap='RdBu',
                 levels=levels)  # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='%')
    plt.plot(extent, extent, 'bo')
    contours = plt.contour((dic['t'] / dic['cnte']), extend='both',
                           levels=[-0.8, -0.7, -0.6, -0.5, -0.4, -0.2, 0, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8],
                           colors='k',linestyles = 'solid', linewidths = 1)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.1f')
    #qu = ax1.quiver(xquiv, yquiv, u, v, scale=50)
    plt.axvline(x=200, linestyle='dashed', color='k', linewidth=1.2)
    plt.axhline(y=200, linestyle='dashed', color='k', linewidth=1.2)
    plt.plot(extent, extent, 'ro')
    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('AMSR-E Day0, contours: 925hPa temperature', fontsize=9)

    ax1 = f.add_subplot(223)
    plt.contourf(((dic['q'])*1000/ dic['cnte']), extend='both',  cmap='RdBu',levels=np.arange(-1,1.1,0.1)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='g kg-1')
    contours = plt.contour((dic['shear'] / dic['cnte']), extend='both',levels=np.arange(-16.5,-11.5,1), colors='k', linestyles='solid', linewidths=1) #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    plt.axvline(x=200, linestyle='dashed', color='k', linewidth=1.2)
    plt.axhline(y=200, linestyle='dashed', color='k', linewidth=1.2)
    plt.plot(extent, extent, 'ro')
    #qu = ax1.quiver(xquiv, yquiv, u, v, scale=50)
    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('shading: 925hPa q anomaly, contours: 650hPa-925hPa wind shear ', fontsize=9)

    ax1 = f.add_subplot(224)
 #   plt.contourf(((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r', vmin=-1.5, vmax=1.5) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.contourf(((dic['div'])/ dic['cnte'])*100, extend='both',  cmap='RdBu', levels=np.linspace(-0.7,0.7,12)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='K')
    plt.axvline(x=200, linestyle='dashed', color='k', linewidth=1.2)
    plt.axhline(y=200, linestyle='dashed', color='k', linewidth=1.2)
    plt.plot(extent, extent, 'ro')
    # contours = plt.contour((dic['shear'] / dic['cntp']), extend='both',levels=np.arange(-17,-12,0.5), cmap='viridis') #np.arange(-15,-10,0.5)
    # plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    qu = ax1.quiver(xquiv, yquiv, u, v, scale=15)
    qk = plt.quiverkey(qu, 0.9, 0.02,1, '1 m s$^{-1}$',
                       labelpos='E', coordinates='figure')

    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('shading: divergence, vectors: 925hPa wind anomaly', fontsize=9)


    plt.tight_layout()
    plt.show()
    #plt.savefig(cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/plots/ERA5_SMALL_"+name+"_"+str(h).zfill(2)+'_'+str(eh).zfill(2)+'.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png)
    #plt.close()

def plot_amsr_paper_trio(h, eh):

    dic = {}
    dic2 = {}
    dic3={}
    dic4 = {}


    name2='ERA5_cores_2hOverlap_AMSRE_SMALL_MCSfilter_minusMean_smallDomain_'
    #
    name3 = 'ERA5_cores_2hOverlap_AMSRE_SMALL_AMSL_MCSfilter'
    name="ERA5_cores_2hOverlap_AMSRE_SMALL_MCSfilter_minusMean_bigDomain_init400km_"

    # name ='ERA5_cores_NEWTRACKING_AMSRE_DRY2_'
    # name2 ='ERA5_cores_NEWTRACKING_AMSRE_DRY2_'
    # name3 ='ERA5_cores_NEWTRACKING_AMSRE_DRY2_'

    # name ='ERA5_cores_NEWTRACKING_AMSRE_WET2f_'
    # name2 ='ERA5_cores_NEWTRACKING_AMSRE_WET2f_'
    # name3 ='ERA5_cores_NEWTRACKING_AMSRE_WET2f_'

    #name ='ERA5_cores_NEWTRACKING_AMSRE_ALL222_'
    #name2 ='ERA5_cores_NEWTRACKING_AMSRE_ALL222_'
    name4 ='ERA5_cores_NEWTRACKING_AMSRE_ALL222_'



    def coll(dic2, h, eh, year, name):
        print(h)
        core = pkl.load(open(
            cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/core_txt/" + name + str(eh) + "UTCERA" + str(h).zfill(
                2) + '_' + str(year) + ".p", "rb"))

        # core = pkl.load(open(
        #     cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/" + name + str(eh) + "UTCERA" + str(h).zfill(
        #         2) + '_' + str(year) + "_small_cores.p", "rb"))

        for id, k in enumerate(core.keys()):
            try:
                dic2[k] = dic2[k] + core[k]
            except KeyError:
                dic2[k] = core[k]

    for y in range(2006, 2010):
        coll(dic, h, eh, y, name)

    for y in range(2006, 2010):
        coll(dic2, h, eh, y, name2)

    for y in range(2006, 2010):
        coll(dic3, h, eh, y, name3)

    for y in range(2006, 2010):
        coll(dic4, h, eh, y, name4)


    print('NUMBER CORES ', np.max(dic4['cnte']))


    extent = (dic['lsta0'].shape[1]-1)/2
    xlen = dic['lsta0'].shape[1]
    ylen = dic['lsta0'].shape[0]

    xv, yv = np.meshgrid(np.arange(ylen), np.arange(xlen))
    st=30
    xquiv = xv[4::st, 4::st]
    yquiv = yv[4::st, 4::st]

    u = (dic2['u925']/ dic2['cnte'])[4::st, 4::st]
    v = (dic2['v925']/ dic2['cnte'])[4::st, 4::st]
    levels = list(np.arange(-2.7,-0.1,0.3)) + list(np.arange(0.3,2.71,0.3))
    qlevels = list(np.round(np.arange(-0.55, -0.01, 0.05),2)) + list(np.round(np.arange(0.05, 0.56, 0.05),2))
    divlevels = list(np.arange(-0.8, -0.01, 0.1)) + list(np.arange(0.1, 0.805, 0.1))
    f = plt.figure(figsize=(15.8,4), dpi=250)


    ax1 = f.add_subplot(131)
    amsr = ndimage.gaussian_filter(((dic2['lsta0']) / (dic2['cnt0'])), 2, mode='nearest')
    vals = np.percentile(amsr,[8,92])
    mask = (amsr <= vals[0]) | (amsr >= vals[1])

    plt.contourf(amsr, extend='both', cmap='RdBu',
                 levels=levels)  # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=12)
    cbar.ax.set_ylabel(r'%', size=12)

    # plt.contour(amsr, extend='both', colors='k',
    #              levels=levels, linewidths=0.05)
    cb=plt.contourf(mask, colors='none', hatches='.',
                 levels=[0.5, 1])
    for i, collection in enumerate(cb.collections):
        collection.set_edgecolor('peachpuff')
    for collection in cb.collections:
        collection.set_linewidth(0.)

    # stpos = np.where(mask)
    # plt.scatter(stpos[1][::60], stpos[0][::60], color='white', s=0.5)

    lev = [-0.6, -0.5, -0.4,-0.3, -0.2, 0, 0.2,0.3, 0.4, 0.5,0.6]
    #lev = list(np.arange(-2.7,-0.1,0.3)) + list(np.arange(0.3,2.71,0.3))
    t = ndimage.gaussian_filter((dic['t']-dic['tclim'])/ dic['cnte'], 6, mode='nearest')
    contours = plt.contour(t, extend='both', #(dic['t']-dic['tclim'])
                           levels=lev,
                           colors='k',linestyles = 'solid', linewidths = 1)


    # t = ndimage.gaussian_filter((dic['skt'])/ dic['cnte'], 6, mode='nearest')
    # contours = plt.contour(t, extend='both', levels=np.linspace(-1.5,1.6,20),
    #                        colors='r',linestyles = 'solid', linewidths = 1)


    plt.clabel(contours, inline=True, fontsize=11, fmt='%1.1f')
    #qu = ax1.quiver(xquiv, yquiv, u, v, scale=50)
    plt.axvline(x=200, linestyle='dashed', color='dimgrey', linewidth=1.2)
    plt.axhline(y=200, linestyle='dashed', color='dimgrey', linewidth=1.2)
    plt.plot(extent, extent, marker='o', color='dimgrey')
    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -200) * 3, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) -200) * 3, dtype=int))
    ax1.set_xlabel('km', fontsize=12)
    ax1.set_ylabel('km', fontsize=12)
    #plt.title('AMSR-E Day0, contours: 925hPa temperature anomaly', fontsize=9)

    ax1 = f.add_subplot(133)
    qq = ((dic['q']-dic['qclim'])*1000/ dic['cnte'])
    qq = ndimage.gaussian_filter(qq, 3, mode='nearest')
    vals = np.percentile(qq,[9,95])
    #vals = dic3['qclim']*1000 / dic3['cnte']

    #tval, pval = ustats.welch_t_test(dic['q']/dic['cnte'], dic3['q']/dic3['cnte'], 5000, dic['qclim']/dic3['cnte'], dic3['qclim']/dic3['cnte'],50)
    mask =  qq>=vals[1]


    plt.contourf(qq, extend='both',  cmap='RdBu',levels=qlevels) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=12)
    cbar.ax.set_ylabel(r'g kg$^{-1}$', size=12)

    # plt.contour(qq, extend='both', colors='k',
    #              levels=qlevels, linewidths=0.05)

    cb=plt.contourf(mask, colors='none', hatches='.',
                 levels=[0.5, 1])
    for i, collection in enumerate(cb.collections):
        collection.set_edgecolor('peachpuff')
    for collection in cb.collections:
        collection.set_linewidth(0.)

    shear = ndimage.gaussian_filter((dic3['shear'] / dic3['cnte']), 3, mode='nearest')
    contours = plt.contour(shear, extend='both',levels=np.arange(-15.25,-12.75,0.5), colors='k', linestyles='solid', linewidths=1) #np.arange(-15,-10,0.5)

    # contours = plt.contour((dic2['u650']-dic2['u650_clim']) / dic2['cnte'], extend='both', levels=np.arange(-2, 2, 0.3), colors='k',
    #                        linestyles='solid', linewidths=1)  # np.arange(-15,-10,0.5)

    plt.clabel(contours, inline=True, fontsize=11, fmt='%1.2f')
    plt.axvline(x=200, linestyle='dashed', color='dimgrey', linewidth=1.2)
    plt.axhline(y=200, linestyle='dashed', color='dimgrey', linewidth=1.2)
    plt.plot(extent, extent, marker='o', color='dimgrey')
    #qu = ax1.quiver(xquiv, yquiv, u, v, scale=50)
    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -200) * 3, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) -200) * 3, dtype=int))
    ax1.set_xlabel('km', fontsize=12)
    #ax1.set_ylabel('km')
    #plt.title('shading: 925hPa q anomaly, contours: 650hPa-925hPa wind shear ', fontsize=9)

    ax1 = f.add_subplot(132)
 #   plt.contourf(((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r', vmin=-1.5, vmax=1.5) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    div = ndimage.gaussian_filter(((dic2['div'])/ dic2['cnte'])*100, 3, mode='nearest')
    vals = np.percentile(div,[9,91])
    #ipdb.set_trace()
    plt.contourf(div, extend='both',  cmap='RdBu', levels=divlevels) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    mask =  (div <= vals[0])
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=12)
    cbar.ax.set_ylabel(r'10$^{-3}$s$^{-1}$', size=11)

    # plt.contour(div, extend='both', colors='k',
    #              levels=divlevels, linewidths=0.05)

    cb = plt.contourf(np.arange(mask.shape[1])[::15], np.arange(mask.shape[0])[::15], mask[::15,::15], colors='none', hatches='.',
                 levels=[0.5, 1])
    for i, collection in enumerate(cb.collections):
        collection.set_edgecolor('peachpuff')
    for collection in cb.collections:
        collection.set_linewidth(0.)

    v925 = ndimage.gaussian_filter(dic3['v925_orig'] / dic3['cnte'], 3, mode='nearest')
    plt.contour(v925, extend='both', colors='k', linewidths=5,
                           levels=[-50, 0.06, 50])

    #itd = ndimage.gaussian_filter(dic3['itd'] / dic3['cnte'], 8, mode='nearest')

    #
    #itd = xr.DataArray(np.array((dic4['itd']/dic4['cnte'])[1::,1::]).astype(float))
    itd = xr.DataArray(np.array((dic4['itd'] / dic4['cnte'])).astype(float))
    # cnt = xr.DataArray(np.array(dic3['cnte'][1::, 1::]).astype(float))
    #
    #itd = salem.reduce(itd, factor=8,how=np.sum)
    # cnt = salem.reduce(cnt, factor=16,how=np.max)

    itd = ndimage.gaussian_filter((itd)*100, 28, mode='nearest')

    #ipdb.set_trace()
    #cs = plt.contour(np.arange(0,400, 8),np.arange(0,400,8),itd, extend='both', colors='r', linewidths=2) #
    cs = plt.contour( itd, extend='both', colors='k', linewidths=0.5, levels=[24.9,30,35], linestyles='dashed')
    cs = plt.contourf(itd, colors='lightgrey',levels=[24.9, 31, 35], alpha=0.27)
    plt.clabel(cs, inline=True, fontsize=11, fmt='%1.0f')

    plt.axvline(x=200, linestyle='dashed', color='dimgrey', linewidth=1.2)
    plt.axhline(y=200, linestyle='dashed', color='dimgrey', linewidth=1.2)
    plt.plot(extent, extent, marker='o', color='dimgrey')
    # contours = plt.contour((dic['shear'] / dic['cntp']), extend='both',levels=np.arange(-17,-12,0.5), cmap='viridis') #np.arange(-15,-10,0.5)
    # plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    qu = ax1.quiver(xquiv, yquiv, u, v, scale=15, headwidth=5)
    qk = plt.quiverkey(qu, 0.9, 0.03,1, r'1 m s$^{-1}$',
                       labelpos='E', coordinates='figure')

    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -200) * 3, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) -200) * 3, dtype=int))
    ax1.set_xlabel('km', fontsize=12)
    #ax1.set_ylabel('km')
    #plt.title('shading: divergence, vectors: 925hPa wind anomaly', fontsize=9)

    plt.tight_layout()
    text = ['a', 'b', 'c']
    plt.annotate(text[0], xy=(0.006, 0.92), xytext=(0, 4),xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points', fontweight='bold', fontname='Ubuntu', fontsize=16)
    plt.annotate(text[1], xy=(0.33, 0.92), xytext=(0, 4), xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points', fontweight='bold', fontname='Ubuntu', fontsize=16)
    plt.annotate(text[2], xy=(0.66, 0.92), xytext=(0, 4),  xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points', fontweight='bold', fontname='Ubuntu', fontsize=16)

    plt.show()
    plt.savefig(cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/plots/"+name+"_"+str(h).zfill(2)+'_'+str(eh).zfill(2)+'-NEW2.pdf')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png)
    plt.close()

    # plt.figure()
    # plt.imshow(itd*100)

def plot_amsr_paper_quatro(h, eh):

    dic = {}

    name = "ERA5_cores_2hOverlap_AMSRE_SMALL_AMSL_MCSfilter" #"ERA5_composite_cores_AMSRE_500w04_15k_minusMean"  # "ERA5_composite_cores_AMSRE_w1_15k_minusMean"
    #name = "ERA5_cores_DRY_SM0LT-3-1LT-1.5_noMeteosatFilter_AMSRE"
    #name='ERA5_cores_WET_TAG_noMeteosatFilter_AMSRE'
    def coll(dic2, h, eh, year):
        print(h)
        core = pkl.load(open(
            cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/core_txt/" + name + str(eh) + "UTCERA" + str(h).zfill(
                2) + '_' + str(year) + ".p", "rb"))

        # core = pkl.load(open(
        #     cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/" + name + str(eh) + "UTCERA" + str(h).zfill(
        #         2) + '_' + str(year) + "_small_cores.p", "rb"))

        for id, k in enumerate(core.keys()):
            try:
                dic2[k] = dic2[k] + core[k]
            except KeyError:
                dic2[k] = core[k]

    for y in range(2006, 2011):
        coll(dic, h, eh, y)

    extent = (dic['lsta0'].shape[1]-1)/2
    xlen = dic['lsta0'].shape[1]
    ylen = dic['lsta0'].shape[0]

    xv, yv = np.meshgrid(np.arange(ylen), np.arange(xlen))
    st=30
    xquiv = xv[4::st, 4::st]
    yquiv = yv[4::st, 4::st]

    u = (dic['u925']/ dic['cnte'])[4::st, 4::st]
    v = (dic['v925']/ dic['cnte'])[4::st, 4::st]
    levels = list(np.arange(-2.1,0,0.3)) + list(np.arange(0.3,2.4,0.3))

    f = plt.figure(figsize=(9.5,7), dpi=200)
    ax = f.add_subplot(221)
    plt.contourf(((dic['lsta-2']) / (dic['cnt-2'])), extend='both', cmap='RdBu',
                 levels=levels)  # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='%') #-(rkernel2_sum / rcnt_sum)
    plt.plot(extent, extent, 'ro')
    contours = plt.contour((dic['skt'] / dic['cnte']), extend='both', cmap='RdBu_r') #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    ax.streamplot(xv, yv, (dic['u925'] / dic['cnte']), (dic['v925'] / dic['cnte']), density=[0.5, 1])

    #contours = plt.contour(((dic['lsta-2']) / (dic['cnt-2'])), extend='both', cmap='PuOr', levels=[-1.5,-1,-0.5,0,0.5,1,1.5])
    # plt.clabel(contours, inline=True, fontsize=11, fmt='%1.0f')


    #contours = plt.contour((dic['v925']/ dic['cnte']), extend='both',colors='k', linewidths=5, levels=[-50,0.06,50])  #np.arange(-15,-10,0.5)

    plt.axvline(x=200, linestyle='dashed', color='k', linewidth=1.2)
    plt.axhline(y=200, linestyle='dashed', color='k', linewidth=1.2)
    plt.plot(extent, extent, 'ro')
    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6 , dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.title('1700UTC | '+str(np.max(dic['cnte']))+' cores, AMSR-E Day-1', fontsize=9)


    ax1 = f.add_subplot(222)
    plt.contourf(((dic['lsta0']) / (dic['cnt0'])), extend='both', cmap='RdBu',
                 levels=levels)  # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='%')
    plt.plot(extent, extent, 'bo')
    contours = plt.contour((dic['t'] / dic['cnte']), extend='both',
                           levels=[-0.8, -0.7, -0.6, -0.5, -0.4, -0.2, 0, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8],
                           colors='k',linestyles = 'solid', linewidths = 1)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.1f')
    #qu = ax1.quiver(xquiv, yquiv, u, v, scale=50)
    plt.axvline(x=200, linestyle='dashed', color='k', linewidth=1.2)
    plt.axhline(y=200, linestyle='dashed', color='k', linewidth=1.2)

    contours = plt.contour((dic['v925_orig']/ dic['cnte']), extend='both',colors='k', linewidths=5, levels=[-50,0.06,50])

    #contours = plt.contour((dic['v925']/ dic['cnte']), extend='both',colors='k', linewidths=3, levels=[-50,0.06,50])  #np.arange(-15,-10,0.5)
    contours = plt.contour((dic['v925_orig']/ dic['cnte']), extend='both',colors='r', linewidths=3, levels=[-50,0.06,50])
    ax1.streamplot(xv, yv, (dic['u925_orig'] / dic['cnte']), (dic['v925_orig'] / dic['cnte']), density=[0.5, 1])
    plt.plot(extent, extent, 'ro')
    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('AMSR-E Day0, contours: 925hPa temperature', fontsize=9)

    ax1 = f.add_subplot(223)
    plt.contourf(((dic['q']-dic['qclim'])*1000/ dic['cnte']), extend='both',  cmap='RdBu',levels=np.arange(-1,1.1,0.1)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='g kg-1')
    contours = plt.contour((dic['shear'] / dic['cnte']), extend='both',levels=np.arange(-16.5,-11.5,1), colors='k', linestyles='solid', linewidths=1) #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    plt.axvline(x=200, linestyle='dashed', color='k', linewidth=1.2)
    plt.axhline(y=200, linestyle='dashed', color='k', linewidth=1.2)
    ax1.streamplot(xv, yv, (dic['u650'] / dic['cnte']), (dic['v650'] / dic['cnte']), density=[0.5, 1])
    plt.plot(extent, extent, 'ro')
    #qu = ax1.quiver(xquiv, yquiv, u, v, scale=50)
    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('shading: 925hPa q anomaly, contours: 650hPa-925hPa wind shear ', fontsize=9)

    ax1 = f.add_subplot(224)
 #   plt.contourf(((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r', vmin=-1.5, vmax=1.5) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.contourf(((dic['div'])/ dic['cnte'])*100, extend='both',  cmap='RdBu', levels=np.linspace(-0.7,0.7,12)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='K')
    plt.axvline(x=200, linestyle='dashed', color='k', linewidth=1.2)
    plt.axhline(y=200, linestyle='dashed', color='k', linewidth=1.2)
    plt.plot(extent, extent, 'ro')
    # contours = plt.contour((dic['shear'] / dic['cntp']), extend='both',levels=np.arange(-17,-12,0.5), cmap='viridis') #np.arange(-15,-10,0.5)
    # plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    ax1.streamplot(xv, yv, (dic['u650_orig'] / dic['cnte']), (dic['v650_orig']*10 / dic['cnte']), density=[0.5, 1])

    contours = plt.contour((dic['probmsg'] / dic['cntm']) * 100, extend='both', cmap='PuOr_r',
                           levels=[-30, -20, -10, 0, 10, 20, 30])
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')

    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('shading: divergence, vectors: 925hPa wind anomaly', fontsize=9)


    plt.tight_layout()
    plt.show()
    plt.savefig(cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/plots/ERA5_SMALL_quatro_"+name+"_"+str(h).zfill(2)+'_'+str(eh).zfill(2)+'.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png)
    plt.close()



def plot_all_quatro():

    for eh in np.arange(-50, 31, 3):
        plot_amsr_paper_quatro(17, eh)


#
#
#
def plot_amsr_paper_diff(h, eh):

    name = "ERA5_cores_DRY_TAG_noMeteosatFilter_AMSRE"#"ERA5_composite_cores_AMSRE_w1_15k_minusMean"
    name2 = 'ERA5_cores_WET_TAG_noMeteosatFilter_AMSRE'#"ERA5_cores_WET_SM0GT0.18_SM-1GT0.1_noMeteosatFilter_AMSRE"

    dic={}
    dic2 ={}

    def coll(dic, h, eh, year,name):
        print(h)
        core = pkl.load(open(
            cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/core_txt/"+name+str(eh) + "UTCERA"+str(h).zfill(2)+'_'+str(year)+".p", "rb"))
        for id, k in enumerate(core.keys()):
            try:
                dic[k] = dic[k] + core[k]
            except KeyError:
                dic[k] = core[k]


    for y in range(2006,2011):
        coll(dic, h, eh, y,name)

    for y in range(2006,2011):
        coll(dic2, h, eh, y,name2)


    extent = (dic['lsta0'].shape[1]-1)/2
    xlen = dic['lsta0'].shape[1]
    ylen = dic['lsta0'].shape[0]

    xv, yv = np.meshgrid(np.arange(ylen), np.arange(xlen))
    st=30
    xquiv = xv[4::st, 4::st]
    yquiv = yv[4::st, 4::st]

    u = ((dic['u925']/dic['cnte'])-(dic2['u925']/dic2['cnte']))[4::st, 4::st]
    v = ((dic['v925'] / dic['cnte']) - (dic2['v925'] / dic2['cnte']))[4::st, 4::st]
    levels = list(np.arange(-3,0,0.5)) + list(np.arange(0.5,3,0.5))

    f = plt.figure(figsize=(9.5,7), dpi=200)
    ax = f.add_subplot(221)

    plt.contourf((dic['lsta-2']/dic['cnt-2'])-(dic2['lsta-2']/dic2['cnt-2']), extend='both', cmap='RdBu', levels=levels)  # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='%') #-(rkernel2_sum / rcnt_sum)
    plt.plot(extent, extent, 'ro')

    #pdb.set_trace()


    #contours = plt.contour(((dic['lsta-2']) / (dic['cnt-2'])), extend='both', cmap='PuOr', levels=[-1.5,-1,-0.5,0,0.5,1,1.5])
    # plt.clabel(contours, inline=True, fontsize=11, fmt='%1.0f')


    contours = plt.contour((dic['v925']/ dic['cnte'])-(dic2['v925']/dic2['cnte']), extend='both',colors='k', linewidths=5, levels=[-50,0.06,50])  #np.arange(-15,-10,0.5)
    contours = plt.contour((dic['v925_orig']/dic['cnte'])-(dic2['v925_orig']/ dic2['cnte']), extend='both',colors='k', linewidths=5, levels=[-50,0.06,50])

    contours = plt.contour((dic['v925']/dic['cnte'])-(dic2['v925']/dic2['cnte']), extend='both',colors='k', linewidths=3, levels=[-50,0.06,50])  #np.arange(-15,-10,0.5)
    contours = plt.contour((dic['v925_orig']/dic['cnte'])-(dic2['v925_orig']/dic2['cnte']), extend='both',colors='r', linewidths=3, levels=[-50,0.06,50])
    plt.axvline(x=200, linestyle='dashed', color='k', linewidth=1.2)
    plt.axhline(y=200, linestyle='dashed', color='k', linewidth=1.2)
    plt.plot(extent, extent, 'ro')
    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6 , dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.title('1700UTC | '+str(np.max(dic['cnte']))+' cores, AMSR-E Day-1', fontsize=9)


    ax1 = f.add_subplot(222)
    plt.contourf((dic['lsta0']/(dic['cnt0'])-(dic2['lsta0']/dic2['cnt0'])), extend='both', cmap='RdBu', levels=levels)  # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='%')
    plt.plot(extent, extent, 'bo')
    contours = plt.contour((dic['t']/dic['cnte'])-(dic2['t'] / dic2['cnte']), extend='both',
                           levels=[-0.8, -0.7, -0.6, -0.5, -0.4, -0.2, 0, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8],
                           colors='k',linestyles = 'solid', linewidths = 1)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.1f')
    #qu = ax1.quiver(xquiv, yquiv, u, v, scale=50)
    plt.axvline(x=200, linestyle='dashed', color='k', linewidth=1.2)
    plt.axhline(y=200, linestyle='dashed', color='k', linewidth=1.2)
    plt.plot(extent, extent, 'ro')
    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('AMSR-E Day0, contours: 925hPa temperature', fontsize=9)

    ax1 = f.add_subplot(223)
    plt.contourf((dic['q']*1000/ dic['cnte'])-(dic['q']*1000/dic2['cnte']), extend='both',  cmap='RdBu',levels=np.arange(-0.5,0.5,0.05)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='g kg-1')
    contours = plt.contour((dic['shear']/dic['cnte'])-(dic2['shear']/ dic2['cnte']), extend='both',levels=np.arange(-1,1,0.1), colors='k', linestyles='solid', linewidths=1) #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    plt.axvline(x=200, linestyle='dashed', color='k', linewidth=1.2)
    plt.axhline(y=200, linestyle='dashed', color='k', linewidth=1.2)
    plt.plot(extent, extent, 'ro')
    #qu = ax1.quiver(xquiv, yquiv, u, v, scale=50)
    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('shading: 925hPa q anomaly, contours: 650hPa-925hPa wind shear ', fontsize=9)

    ax1 = f.add_subplot(224)
 #   plt.contourf(((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r', vmin=-1.5, vmax=1.5) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.contourf(((dic['div']/dic['cnte'])-(dic2['div']/ dic2['cnte']))*100, extend='both',  cmap='RdBu', levels=np.linspace(-0.5,0.5,12)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='K')
    plt.axvline(x=200, linestyle='dashed', color='k', linewidth=1.2)
    plt.axhline(y=200, linestyle='dashed', color='k', linewidth=1.2)
    plt.plot(extent, extent, 'ro')
    # contours = plt.contour((dic['shear'] / dic['cntp']), extend='both',levels=np.arange(-17,-12,0.5), cmap='viridis') #np.arange(-15,-10,0.5)
    # plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    qu = ax1.quiver(xquiv, yquiv, u, v, scale=15)
    qk = plt.quiverkey(qu, 0.9, 0.02,1, '1 m s$^{-1}$',
                       labelpos='E', coordinates='figure')

    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('shading: divergence, vectors: 925hPa wind anomaly', fontsize=9)


    plt.tight_layout()
    plt.show()
    plt.savefig(cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/plots/ERA5_DIFF_"+name+"_"+str(h).zfill(2)+'_'+str(eh).zfill(2)+'.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png)
    plt.close()

