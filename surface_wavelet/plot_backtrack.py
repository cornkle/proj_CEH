# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 10:15:40 2016

@author: cornkle
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib
import pdb
import pandas as pd
from utils import u_met, u_parallelise, u_gis, u_arrays, constants, u_grid
import salem
from scipy.interpolate import griddata

import pickle as pkl


matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)

def plot_gewex(lag):
    hour=21
    path = '/users/global/cornkle/figs/LSTA-bullshit/corrected_LSTA/lsta_sameHour/passage_plots/'
    dic = pkl.load(open(path+"composite_backtrack_"+str(lag)+"_"+str(hour).zfill(2)+".p", "rb"))

    extent = (dic['ano'].shape[1]-1)/2

    f = plt.figure(figsize=(7, 5))
    ax = f.add_subplot(111)
    pdb.set_trace()
    plt.contourf((dic['ano'] / dic['cnt']), cmap='RdBu_r',  levels=[-0.5,-0.4,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5], extend='both') #-(rkernel2_sum / rcnt_sum)
    #plt.plot(extent, extent, 'bo')
    plt.colorbar(label='K')

    contours = plt.contour((dic['prob']/ dic['cnt']) * 100, extend='both', levels=np.arange(25,101,10)) # #, levels=np.arange(1,5, 0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.0f')

    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')

    plt.title(str(hour).zfill(2)+'00 UTC | '+str(np.max(dic['cnt']))+' cores', fontsize=17)


    plt.tight_layout()
    # plt.savefig('/users/global/cornkle/figs/LSTA-bullshit/GEWEX/'+str(hour).zfill(2)+'_single.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png')
    # plt.close()

def plot_gewex_double(h):
    hour=h
    dic = pkl.load(open("/users/global/cornkle/figs/LSTA-bullshit/corrected_LSTA/system_scale/composite_backtrack_"+str(hour).zfill(2)+".p", "rb"))

    extent = (dic['ano'].shape[1]-1)/2

    f = plt.figure(figsize=(12, 5))
    ax = f.add_subplot(121)
#
    plt.contourf((dic['ano']/ dic['cnt']), cmap='RdBu_r',  levels=[-0.5,-0.4,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5], extend='both') #-(rkernel2_sum / rcnt_sum)
    plt.colorbar(label='K')

    contours = plt.contour((dic['prob']/ dic['pcnt']) * 100, extend='both', levels=np.arange(25,101,5)) # #, levels=np.arange(1,5, 0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.1f')

    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')

    plt.title('17-19UTC | '+str(np.max(dic['cnt']))+' cores, Monthly LSTA', fontsize=17)

    ax = f.add_subplot(122)

    #plt.pcolormesh( dic['cnt'], cmap='viridis') #-(rkernel2_sum / rcnt_sum)


    contours = plt.pcolormesh((dic['cnt']) ) # #, levels=np.arange(1,5, 0.5)
    #plt.clabel(contours, inline=True, fontsize=9, fmt='%1.1f')

    #plt.contourf((dic['prob']/ dic['pcnt'])*100)
    #plt.colorbar(label='%')
    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.title('17-19UTC | ' + str(np.max(dic['cnt'])) + ' cores, Regional LSTA', fontsize=17) #str(hour).zfill(2) + '

    plt.tight_layout()
    plt.savefig('/users/global/cornkle/figs/LSTA-bullshit/corrected_LSTA/system_scale/composite_backtrack_'+str(hour).zfill(2)+'.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png')
    plt.close()


def plot_all():

    hours = [0,1,2,3,4,5,6,7,8,9,17,18,19,20,21,22,23]

    for h in hours:

        plot_gewex(h)

