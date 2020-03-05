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
import glob
import os
from collections import OrderedDict
import salem
from utils import u_met, u_parallelise, u_gis, u_arrays as ua, constants as cnst, u_grid, u_darrays
from scipy.interpolate import griddata
import multiprocessing
import metpy
from metpy import calc
from metpy.units import units

import pickle as pkl


matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)



def plot_timeseries_coarse():

    xx = len(list(range(-38, 4, 3)))
    yy = 401

    ranges = np.arange(-38, 4, 3)

    h = 17

    def coll(dic, h, eh, year, ff):
        print(h)
        core = pkl.load(open(
            cnst.network_data + ff + str(eh) + "UTCERA" + str(h).zfill(2) + '_' + str(year) + ".p", "rb"))
        for id, k in enumerate(core.keys()):
            try:
                dic[k] = dic[k] + core[k]
            except KeyError:
                dic[k] = core[k]

    filewet = "figs/LSTA/corrected_LSTA/new/ERA5/core_txt/ERA5_cores_2hOverlap_AMSRE_SMWET_AMSL_mini"

    filedry = "figs/LSTA/corrected_LSTA/new/ERA5/core_txt/ERA5_cores_2hOverlap_AMSRE_SMDRY_AMSL_mini"

    file = "figs/LSTA/corrected_LSTA/new/ERA5/core_txt/ERA5_cores_2hOverlap_AMSRE_SMALL_MCSfilter_minusMean_smallDomain_"#ERA5_cores_2hOverlap_AMSRE_SMALL_AMSL_new"

    #file = "figs/LSTA/corrected_LSTA/new/ERA5/core_txt/ERA5_cores_DRY_SM0LT3-1LT1.5_noMeteosatFilter_AMSRE"
    #file = "figs/LSTA/corrected_LSTA/new/ERA5/core_txt/ERA5_cores_WET_SM0GT0.01-1GT0.01_noMeteosatFilter_AMSRE"
    #####file = "figs/LSTA/corrected_LSTA/new/ERA5/core_txt/ERA5_cores_WET_TAG_noMeteosatFilter_AMSRE"
    fs = [file,filedry,filewet]

    f = plt.figure(figsize=(14,7), dpi=300)


    for fids, ff in enumerate(fs):

        outdic = {}

        dummyfile = glob.glob(
            cnst.network_data + ff+"*.p")
        dummy = pkl.load(open(dummyfile[0], "rb"))

        for k in dummy.keys():
            outdic[k] = np.zeros((yy, xx))

        for ids, eh in enumerate(range(-38, 4, 3)):
            dic = {}

            for y in range(2006, 2011):
                coll(dic, h, eh, y, ff)
            #if fids == 1:
                #ipdb.set_trace()
            for k in dic.keys():
                if (k=='lsta0') & (eh<-11):
                    outdic[k][:, ids] = dic['lsta-1'][:, 190:211].mean(axis=1)
                else:
                    outdic[k][:, ids] = dic[k][:, 190:211].mean(axis=1) # [190:211]

        outdic['lsta0'][0:25, :] = 0

        #

        def groupedAvg(myArray, N=2):
            result = np.cumsum(myArray, 0)[N - 1::N, :] / float(N)
            result[1:,:] = result[1:, :] - result[:-1, :]
            return result

        diff = {}
        for k in outdic.keys():

            if 'cnt' in k:
                continue

            if 'lsta0' in k:
                diff[k] = groupedAvg((outdic[k] / outdic['cnt0']) , N=4)
            elif 'lsta-1' in k:
                diff[k] = groupedAvg((outdic[k] / outdic['cnt-1']) , N=4)
            elif 'lsta-2' in k:
                diff[k] = groupedAvg((outdic[k] / outdic['cnt-2']), N=4)
            elif 'lsta-3' in k:
                diff[k] = groupedAvg((outdic[k] / outdic['cnt-3']) , N=4)
            elif 'msg' in k:
                diff[k] = groupedAvg((outdic[k] / outdic['cntm']) , N=4)
            elif 'probc' in k:
                diff[k] = groupedAvg((outdic[k] / outdic['cntc']) , N=4)
            else:
                diff[k] = groupedAvg((outdic[k] / outdic['cnte']) , N=4)


        print(diff.keys())

        titles = ['ALL', 'DRY', 'WET']

        ax = f.add_subplot(2,3,fids+1)

        yax = np.arange(100)

        levels = list(np.arange(-0.9, 0, 0.1)) + list(np.arange(0.1, 0.95, 0.1))
        plt.contourf(ranges, yax, (diff['t'])-diff['tclim'], extend='both', cmap='RdBu_r',
                    levels=levels)
        plt.colorbar(label=r'K')
        plt.contour(ranges, yax, (diff['t'])-diff['tclim'], extend='both', colors='k', linewidths=0.1,
                     levels=levels)
        # if (fids == 1):
        #     cont = plt.contour(ranges, yax, diff['lsta0']*-1, extend='both', colors='r', linewidths=2, levels=[1.5,5], linestyle='solid')
        #
        # if (fids == 2):
        #     diff['lsta0'][75:100, :] = 0
        #     #ipdb.set_trace()
        #     cont = plt.contour(ranges, yax, diff['lsta0'], extend='both', colors='b', linewidths=2,
        #                        levels=[0.4, 1.5,10], linestyle='solid')

        #plt.clabel(cont, inline=True, fontsize=9, fmt='%1.2f')

        contours = plt.contour(ranges, yax, (diff['u650']), extend='both', colors='k',
                               levels=np.arange(-4,0,0.5), linewidths=1)


        plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')

        contours = plt.contour(ranges, yax,(diff['v925_orig']), extend='both',colors='k', linewidths=5, levels=[-50,0.06,50])

        plt.axvline(x=-5, color='slategrey')
        plt.axvline(x=-29, color='slategrey')
        plt.axhline(y=50, color='slategrey')
        plt.plot(-5, 50, 'ko')
        plt.plot(0, 50, 'ro')

        ax.set_yticklabels(np.array((np.linspace(0, 100, 6) - 50) * 12, dtype=int))
        if (fids == 0):
            plt.ylabel('North-South distance from core (km)')

        plt.title(titles[fids], fontweight='bold', fontname='Ubuntu', fontsize=12, ha='left', loc='left')


        ax = f.add_subplot(2,3,fids+4)

        plt.contourf(ranges, yax, (diff['q']-diff['qclim']) * 1000, levels=[-0.8, -0.7, -0.6, -0.5, -0.4, -0.3,-0.2, -0.1, 0.1, 0.2, 0.3,0.4, 0.5, 0.6, 0.7, 0.8],
                     cmap='RdBu', extend='both')
        plt.colorbar(label=r'g kg$^{-1}$')

        plt.contour(ranges, yax, (diff['q']-diff['qclim']) * 1000, levels=[-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3,0.4, 0.5, 0.6, 0.7, 0.8],
         colors = 'k', linewidths = 0.1)

        contours = plt.contour(ranges, yax, (diff['theta_e']-diff['theta_e_clim']), colors='white',
                                linewidths=2, levels=[0.8,1.2,1.6,2,2.4]) #levels=np.arange(-3, 3, 0.5),
        plt.clabel(contours, inline=True, fontsize=8, fmt='%1.1f')


        plt.axvline(x=-5, color='slategrey')
        plt.axvline(x=-29, color='slategrey')
        plt.axhline(y=50, color='slategrey')
        plt.plot(-5, 50, 'ko')
        plt.plot(0, 50, 'ro')

        ax.set_yticklabels(np.array((np.linspace(0, 100, 6) - 50) *12 , dtype=int))


        plt.xlabel('Hour relative to t0')

        if (fids == 0):
            plt.ylabel('North-South distance from core (km)')

    plt.tight_layout()

    text = ['a', 'b', 'c', 'd', 'e', 'f']
    fs = 14
    plt.annotate(text[0], xy=(0.006, 0.96), xytext=(0, 4),xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points', fontweight='bold', fontname='Ubuntu', fontsize=fs)
    plt.annotate(text[1], xy=(0.34, 0.96), xytext=(0, 4), xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points', fontweight='bold', fontname='Ubuntu', fontsize=fs)
    plt.annotate(text[2], xy=(0.67, 0.96), xytext=(0, 4),  xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points', fontweight='bold', fontname='Ubuntu', fontsize=fs)
    plt.annotate(text[3], xy=(0.006, 0.48), xytext=(0, 4),xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points', fontweight='bold', fontname='Ubuntu', fontsize=fs)
    plt.annotate(text[4], xy=(0.34, 0.48), xytext=(0, 4), xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points', fontweight='bold', fontname='Ubuntu', fontsize=fs)
    plt.annotate(text[5], xy=(0.67, 0.48), xytext=(0, 4),  xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points', fontweight='bold', fontname='Ubuntu', fontsize=fs)

    plt.savefig(cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/plots/ERA5_amsre_"+str(h).zfill(2)+'_timeseries_ALL2.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png)
    plt.close()

def plot_timeseries_coarse_lon():

    yy = len(list(range(-38, 24, 3)))
    xx = 401

    ranges = np.arange(-38, 24, 3)

    h = 17

    def coll(dic, h, eh, year, ff):
        print(h)
        core = pkl.load(open(
            cnst.network_data + ff + str(eh) + "UTCERA" + str(h).zfill(2) + '_' + str(year) + ".p", "rb"))
        for id, k in enumerate(core.keys()):
            try:
                dic[k] = dic[k] + core[k]
            except KeyError:
                dic[k] = core[k]

    filewet = "figs/LSTA/corrected_LSTA/new/ERA5/core_txt/ERA5_cores_2hOverlap_AMSRE_SMWET_AMSL_mini"

    filedry = "figs/LSTA/corrected_LSTA/new/ERA5/core_txt/ERA5_cores_2hOverlap_AMSRE_SMDRY_AMSL_mini"

    file = "figs/LSTA/corrected_LSTA/new/ERA5/core_txt/ERA5_cores_2hOverlap_AMSRE_SMALL_MCSfilter_minusMean_smallDomain_"  # ERA5_cores_2hOverlap_AMSRE_SMALL_AMSL_new"

    # file = "figs/LSTA/corrected_LSTA/new/ERA5/core_txt/ERA5_cores_DRY_SM0LT3-1LT1.5_noMeteosatFilter_AMSRE"
    # file = "figs/LSTA/corrected_LSTA/new/ERA5/core_txt/ERA5_cores_WET_SM0GT0.01-1GT0.01_noMeteosatFilter_AMSRE"
    #####file = "figs/LSTA/corrected_LSTA/new/ERA5/core_txt/ERA5_cores_WET_TAG_noMeteosatFilter_AMSRE"
    fs = [filedry, filewet]

    f = plt.figure(figsize=(14, 7), dpi=200)

    for fids, ff in enumerate(fs):

        outdic = {}

        dummyfile = glob.glob(
            cnst.network_data + ff + "*.p")
        dummy = pkl.load(open(dummyfile[0], "rb"))

        for k in dummy.keys():
            outdic[k] = np.zeros((yy, xx))

        for ids, eh in enumerate(range(-38, 24, 3)):
            dic = {}

            for y in range(2006, 2011):
                coll(dic, h, eh, y, ff)
            # if fids == 1:
            # ipdb.set_trace()
            for k in dic.keys():
                if (k == 'lsta0') & (eh < -11):
                    outdic[k][ids, :] = dic['lsta-1'][190:211, :].mean(axis=0)
                else:
                    outdic[k][ids, :] = dic[k][190:211, : ].mean(axis=0)  # [190:211]

        #outdic['lsta0'][0:25, :] = 0

        #

        def groupedAvg(myArray, N=2):
            myArray = myArray.T
            result = np.cumsum(myArray, 0)[N - 1::N, :] / float(N)
            result[1:,:] = result[1:, :] - result[:-1, :]
            return result.T

        diff = {}
        for k in outdic.keys():

            if 'cnt' in k:
                continue

            if 'lsta0' in k:
                diff[k] = groupedAvg((outdic[k] / outdic['cnt0']) , N=4)
            elif 'lsta-1' in k:
                diff[k] = groupedAvg((outdic[k] / outdic['cnt-1']) , N=4)
            elif 'lsta-2' in k:
                diff[k] = groupedAvg((outdic[k] / outdic['cnt-2']), N=4)
            elif 'lsta-3' in k:
                diff[k] = groupedAvg((outdic[k] / outdic['cnt-3']) , N=4)
            elif 'msg' in k:
                diff[k] = groupedAvg((outdic[k] / outdic['cntm']) , N=4)
            elif 'probc' in k:
                diff[k] = groupedAvg((outdic[k] / outdic['cntc']) , N=4)
            else:
                diff[k] = groupedAvg((outdic[k] / outdic['cnte']) , N=4)

        print(diff.keys())

        titles = ['DRY', 'WET']

        ax = f.add_subplot(2, 2, fids + 1)

        yax = np.arange(100) #np.arange(401)

        levels = list(np.arange(-0.9, 0, 0.1)) + list(np.arange(0.1, 0.95, 0.1))
        plt.contourf(yax, ranges, (diff['t']) - diff['tclim'], extend='both', cmap='RdBu_r',
                     levels=levels)
        plt.colorbar(label=r'K')
        plt.contour(yax, ranges,  (diff['t']) - diff['tclim'], extend='both', colors='k', linewidths=0.1,
                    levels=levels)
        # if (fids == 1):
        #     cont = plt.contour(ranges, yax, diff['lsta0']*-1, extend='both', colors='r', linewidths=2, levels=[1.5,5], linestyle='solid')
        #
        # if (fids == 2):
        #     diff['lsta0'][75:100, :] = 0
        #     #ipdb.set_trace()
        #     cont = plt.contour(ranges, yax, diff['lsta0'], extend='both', colors='b', linewidths=2,
        #                        levels=[0.4, 1.5,10], linestyle='solid')

        # plt.clabel(cont, inline=True, fontsize=9, fmt='%1.2f')

        contours = plt.contour(yax, ranges, (diff['u650']), extend='both', colors='k',
                               levels=np.arange(-4, 0, 0.5), linewidths=1)

        plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')

        # contours = plt.contour(yax, ranges, (diff['v925_orig']), extend='both', colors='k', linewidths=5,
        #                        levels=[-50, 0.06, 50])

        # plt.axvline(x=-5, color='slategrey')
        plt.axvline(x=50, color='slategrey')
        plt.axhline(y=0, color='slategrey')
        # plt.plot(-5, 50, 'ko')
        # plt.plot(0, 50, 'ro')

        ax.set_xticklabels(np.array((np.linspace(0, 100, 6) - 50) * 12, dtype=int))
        if (fids == 0):
            plt.ylabel('Hour relative to t0')

        plt.title(titles[fids], fontweight='bold', fontname='Ubuntu', fontsize=12, ha='left', loc='left')

        ax = f.add_subplot(2, 2, fids + 3)

        # plt.contourf(yax, ranges, (diff['q'] - diff['qclim']) * 1000,
        #              levels=[-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7,
        #                      0.8],
        #              cmap='RdBu', extend='both')
        # plt.colorbar(label=r'g kg$^{-1}$')
        #
        # plt.contour(yax, ranges, (diff['q'] - diff['qclim']) * 1000,
        #             levels=[-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8],
        #             colors='k', linewidths=0.1)

        plt.contourf(yax, ranges, (diff['v650'] ),
                     levels=np.arange(-2.25,2.26,0.25),
                     cmap='RdBu', extend='both')
        plt.colorbar(label=r'm s$^{-1}$')

        plt.contour(yax, ranges, (diff['v650']),
                    levels=np.arange(-2.25, 2.26, 0.25),
                    colors='k', linewidths=0.1)

        # contours = plt.contour(yax, ranges, (diff['theta_e'] - diff['theta_e_clim']), colors='white',
        #                        linewidths=2, levels=[0.8, 1.2, 1.6, 2, 2.4])  # levels=np.arange(-3, 3, 0.5),
        # plt.clabel(contours, inline=True, fontsize=8, fmt='%1.1f')

        contours=plt.contour(yax, ranges, (diff['q'] - diff['qclim']) * 1000,
                    levels=[-0.8,  -0.6,  -0.4, -0.2, 0, 0.2,  0.4, 0.6,  0.8],
                    colors='white', linewidths=2)
        plt.clabel(contours, inline=True, fontsize=8, fmt='%1.1f')
        # plt.axvline(x=-5, color='slategrey')
        plt.axvline(x=50, color='slategrey')
        plt.axhline(y=0, color='slategrey')
        plt.plot(y=0, x=50, color='k', marker='o')
        # plt.plot(0, 50, 'ro')

        ax.set_xticklabels(np.array((np.linspace(0, 100, 6) - 50) * 12, dtype=int))

        plt.xlabel('East-West distance from core (km)')

        if (fids == 0):
            plt.ylabel('Hour relative to t0')

    plt.tight_layout()

    text = ['a', 'b', 'c', 'd', 'e', 'f']
    fs = 14
    plt.annotate(text[0], xy=(0.006, 0.96), xytext=(0, 4), xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points', fontweight='bold', fontname='Ubuntu', fontsize=fs)
    plt.annotate(text[1], xy=(0.52, 0.96), xytext=(0, 4), xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points', fontweight='bold', fontname='Ubuntu', fontsize=fs)
    plt.annotate(text[2], xy=(0.006, 0.48), xytext=(0, 4), xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points', fontweight='bold', fontname='Ubuntu', fontsize=fs)
    plt.annotate(text[3], xy=(0.52, 0.48), xytext=(0, 4), xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points', fontweight='bold', fontname='Ubuntu', fontsize=fs)


    plt.savefig(cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/plots/ERA5_amsre_LON_" + str(h).zfill(
        2) + '_timeseries_ALL2.png')  # str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png)
    plt.close()
