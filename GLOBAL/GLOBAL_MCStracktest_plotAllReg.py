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
from utils import constants as cnst
import matplotlib.patches as patches
import pickle as pkl


matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)


MREGIONS = {'WAf' : [[-18,25,4,25], 'spac', 0, (1,7), (8,12), (1,12)], # last is hourly offset to UCT # 12    # [-18,25,4,25]
 'SAf' : [[20,35, -35,-15], 'spac', 2, (9,12), (1,5), (1,12)], # 10
 'india' : [[70,90, 5,30], 'asia', 5, (1,7), (8,12), (1,12)], # 7
 'china' : [[105,115,25,40], 'asia', 8 , (1,7), (8,12), (1,12)], # 4
 'australia' : [[120,140,-23, -11], 'asia', 9, (10,12), (1,5), (1,12)], # 3
 'sub_SA' : [[-68,-47, -40, -20.5], 'spac', -4, (9,12), (1,5), (1,12)] , # 16
 'trop_SA' : [[-75, -50, -20, -5], 'spac', -5, (1,12), (1,12), (1,12)], # 17
 'GPlains' : [[-100,-90,32,47], 'nam', -6, (1,7), (8,12), (1,12)] # # 18

}


REGIONS = ['GPlains', 'sub_SA', 'WAf', 'china', 'india', 'australia']
SENSOR = 'AMSR2'

OUT = '/home/ck/DIR/cornkle/figs/GLOBAL_MCS/'


MONTHS=(1,12)
rawhour=2


f = plt.figure(figsize=(23, 8))



for ids, regs in enumerate(REGIONS):

        h = rawhour - (MREGIONS[regs])[2]
        print('Hour: ', h)
        extag = '_day-1'
        pin = cnst.network_data + 'data/GLOBAL_MCS/save_composites/' + "/" + regs + "_SM_" + SENSOR + "_SWA_2012-2019_MCSTRACK_BOX-ANNUAL_PF_DAY_" + str(
            MONTHS[0]).zfill(2) + '-' + str(MONTHS[1]).zfill(2) + '_'
        y1 = 2013
        y2 = 2020

        dic = pkl.load(open(
            pin + str(y1) + '_h' + str(rawhour).zfill(2) + extag + ".p", "rb"))

        def coll(dic, h, year):
            print(h)
            core = pkl.load(open(
                pin + str(year) + '_h' + str(rawhour).zfill(2) + extag + ".p", "rb"))
            for id, k in enumerate(core.keys()):
                try:
                    dic[k] = dic[k] + core[k]
                except KeyError:
                    dic[k] = core[k]

        for y in range(y1 + 1, y2):
            print('Coll', y)
            coll(dic, h, y)

        extent = 11


        ax = f.add_subplot(3,6,ids+1)

        ano = dic['init']#(dic['ano'] / dic['cnt']) - np.mean((dic['ano']) / dic['cnt'])  # / dic['cnt']

        thresh = np.max(np.abs(np.percentile(ano, [5,95])))

        plt.contourf(ano, cmap='rainbow', levels=np.arange(1,20), #np.linspace(thresh * -1, thresh, 10),
                     extend='both')  # -(rkernel2_sum / rcnt_sum)
        plt.plot(extent, extent, 'bo')
        ax.set_xticks((np.linspace(0, 2 * extent, 9)))
        ax.set_xticklabels(((np.linspace(0, (2 * extent), 9) - extent)*0.25).round(1).astype(float))
        ax.set_yticks((np.linspace(0, 2 * extent, 9)))
        ax.set_yticklabels(((np.linspace(0, (2 * extent), 9) - extent)*0.25).round(1).astype(float))
        ax.set_xlabel('Degrees')
        ax.set_ylabel('Degrees')
        ax.axvline(extent, linestyle='dashed', color='k')
        ax.axhline(extent, linestyle='dashed', color='k')

        rect = patches.Rectangle((5.5, 5.5), 11, 11, linewidth=0.5, edgecolor='k', facecolor='none')
        ax.add_patch(rect)

        plt.colorbar(label='count')
        plt.title(regs+' day 0 initiations',  fontsize=10)

        ######################################################

        ax = f.add_subplot(3, 6, ids + 1 + 6)

        ano = (dic['ano'] / dic['cnt']) #- np.mean((dic['ano']) / dic['cnt']) # / dic['cnt']

        thresh = np.max(np.abs(np.percentile(ano, [5, 95])))

        plt.contourf(ano, cmap='RdBu', levels=np.linspace(thresh * -1, thresh, 10),
                     extend='both')
        plt.plot(extent, extent, 'bo')
        ax.set_xticks((np.linspace(0, 2 * extent, 9)))
        ax.set_xticklabels(((np.linspace(0, (2 * extent), 9) - extent)*0.25).round(1).astype(float))
        ax.set_yticks((np.linspace(0, 2 * extent, 9)))
        ax.set_yticklabels(((np.linspace(0, (2 * extent), 9) - extent)*0.25).round(1).astype(float))
        ax.set_xlabel('Degrees')
        ax.set_ylabel('Degrees')
        ax.axvline(extent, linestyle='dashed', color='k')
        ax.axhline(extent, linestyle='dashed', color='k')

        rect = patches.Rectangle((5.5, 5.5), 11, 11, linewidth=0.5, edgecolor='k', facecolor='none')
        ax.add_patch(rect)

        cb = plt.colorbar(label='%', format='%.2f')
        plt.title(regs+' day-1 ano  | n='+str(np.max(dic['allcnt'])),
                  fontsize=10)

        ###############################################


        extag = '_day0'
        pin = cnst.network_data + 'data/GLOBAL_MCS/save_composites/' + "/" + regs + "_SM_" + SENSOR + "_SWA_2012-2019_MCSTRACK_BOX-ANNUAL_PF_DAY_" + str(
            MONTHS[0]).zfill(2) + '-' + str(MONTHS[1]).zfill(2) + '_'

        dic = pkl.load(open(
            pin + str(y1) + '_h' + str(rawhour).zfill(2) + extag + ".p", "rb"))

        def coll(dic, h, year):
            print(h)
            core = pkl.load(open(
                pin + str(year) + '_h' + str(rawhour).zfill(2) + extag + ".p", "rb"))
            for id, k in enumerate(core.keys()):
                try:
                    dic[k] = dic[k] + core[k]
                except KeyError:
                    dic[k] = core[k]

        for y in range(y1 + 1, y2):
            print('Coll', y)
            coll(dic, h, y)

        ax = f.add_subplot(3, 6, ids + 1 + 12)

        ano = (dic['ano'] / dic['cnt'])# - np.mean((dic['ano']) / dic['cnt'])  # / dic['cnt']

        thresh = np.max(np.abs(np.percentile(ano, [5, 95])))

        plt.contourf(ano, cmap='RdBu', levels=np.linspace(thresh * -1, thresh, 10),
                     extend='both')  # -(rkernel2_sum / rcnt_sum)
        plt.plot(extent, extent, 'bo')
        ax.set_xticks((np.linspace(0, 2 * extent, 9)))
        ax.set_xticklabels(((np.linspace(0, (2 * extent), 9) - extent)*0.25).round(1).astype(float))
        ax.set_yticks((np.linspace(0, 2 * extent, 9)))
        ax.set_yticklabels(((np.linspace(0, (2 * extent), 9)  - extent)*0.25).round(1).astype(float))
        ax.set_xlabel('Degrees')
        ax.set_ylabel('Degrees')
        ax.axvline(extent, linestyle='dashed', color='k')
        ax.axhline(extent, linestyle='dashed', color='k')

        rect = patches.Rectangle((5.5, 5.5), 11, 11, linewidth=0.5, edgecolor='k', facecolor='none')
        ax.add_patch(rect)

        plt.colorbar(label='%', format='%.2f')
        plt.title(regs+' day0 ano  | n='+str(np.max(dic['allcnt'])), fontsize=10)

plt.tight_layout()

plt.savefig(
    OUT + 'all_reg/' + str(rawhour) + '/AllRegions_' + SENSOR + '_SM_TS_' + str(MONTHS[0]).zfill(
        2) + '-' + str(
        MONTHS[1]).zfill(2) + '_2012-2019_allPF_' + str(rawhour).zfill(2) + 'h_pureAnom.jpg')
plt.close('all')

