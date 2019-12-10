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

    name = "ERA5_cores_WET_SM0GT0.01-1GT0.01_noMeteosatFilter_AMSRE"  # "ERA5_composite_cores_AMSRE_w1_15k_minusMean"

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
    plt.contourf(((dic['q'])*1000/ dic['cntp']), extend='both',  cmap='RdBu',levels=np.arange(-1,1.1,0.1)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='g kg-1')
    contours = plt.contour((dic['shear'] / dic['cntp']), extend='both',levels=np.arange(-16.5,-11.5,1), colors='k', linestyles='solid', linewidths=1) #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    plt.axvline(x=200, linestyle='dashed', color='k', linewidth=1.2)
    plt.axhline(y=200, linestyle='dashed', color='k', linewidth=1.2)
    plt.plot(extent, extent, 'ro')
    #qu = ax1.quiver(xquiv, yquiv, u, v, scale=50)
    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('Shading: 925hPa q anomaly, Contours: 650hPa-925hPa wind shear ', fontsize=9)

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
    plt.savefig(cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/plots/ERA5_LS_"+str(h).zfill(2)+'_'+str(eh).zfill(2)+'.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png)
    plt.close()

def plot_amsr_paper(h, eh):

    dic = {}

    name = "ERA5_composite_cores_AMSRE_500w04_15k_minusMean"  # "ERA5_composite_cores_AMSRE_w1_15k_minusMean"
    #name = "ERA5_cores_DRY_SM0LT-3-1LT-1.5_noMeteosatFilter_AMSRE"
    #name='ERA5_cores_WET_TAG_noMeteosatFilter_AMSRE'
    def coll(dic2, h, eh, year):
        print(h)
        core = pkl.load(open(
            cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/" + name + str(eh) + "UTCERA" + str(h).zfill(
                2) + '_' + str(year) + "_small_cores.p", "rb"))

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
    plt.savefig(cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/plots/ERA5_"+name+"_"+str(h).zfill(2)+'_'+str(eh).zfill(2)+'.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png)
    plt.close()

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

