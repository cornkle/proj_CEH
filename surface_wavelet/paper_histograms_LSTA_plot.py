# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 10:15:40 2016

@author: cornkle
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib
import pickle as pkl
matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)
import pdb
import multiprocessing
from statsmodels.stats.proportion import proportion_confint
from utils import constants as cnst
import ipdb

def diurnal_loop():

    pool = multiprocessing.Pool(processes=4)
    h = [15,16,17,18,19,20,21,22,23,0,1,2,3,4,5,6,7,8,9,10,11,12,13]
    res = pool.map(plot, h)
    pool.close()



def plot(hour):
    path = cnst.network_data + 'figs/LSTA/corrected_LSTA/new/wavelet_coefficients'
    dic = pkl.load(open(path+"/LSTA_histograms_"+str(hour).zfill(2)+".p", "rb"))


    for k in dic.keys():
        coll = []
        for ll in dic[k]:
            coll.extend(ll)
        dic[k] = coll
    print(dic.keys())

    cinput = np.array(dic['e100'])
    rinput = np.array(dic['re100'])

    cinput = cinput[np.isfinite(cinput)]
    rinput = rinput[np.isfinite(rinput)]

    f = plt.figure()
    ax = f.add_subplot(121)

    nball, bins,v = plt.hist(cinput, bins=np.arange(-7,7,1), normed=True, edgecolor='k', color=None, alpha=0.3)
    nbpoint, bins, v = plt.hist(rinput, bins=np.arange(-7,7,1), normed=True, edgecolor='k', color=None, alpha=0.3)
    plt.xlabel('Local wavelet coefficient, 30km')
    plt.ylabel('Frequency')
    stri = (np.sum(cinput >= np.percentile(rinput, 90)) / cinput.size * 100).round(2)
    plt.title('18-19UTC:'+str(stri)+'% of Cells occur in warmest decile')

    stri = (np.sum(cinput >= np.percentile(rinput, 90)) / cinput.size * 100).round(2)
    print('18-19UTC:' + str(stri) + '% of Cells occur in warmest half')

    stri = (np.sum(cinput >= 0.03) / cinput.size * 100).round(2)
    print('18-19UTC:' + str(stri) + '% of Cells occur over significantly positive LSTA')

    # dic = pkl.load(open("/users/global/cornkle/figs/LSTA-bullshit/scales/new/histo0.p", "rb"))
    #
    # all = dic['all']
    # allgrad = dic['allgrad']
    # point = dic['point']
    # pointgrad = dic['pointgrad']
    #
    # perc = np.percentile(all, np.arange(0, 101, 10))
    # percg = np.percentile(allgrad, np.arange(0, 101, 10))
    # ax = f.add_subplot(122)
    #
    # nball, bins, v = plt.hist(all, bins=np.arange(-10, 10, 1), normed=True, edgecolor='k', color=None, alpha=0.3)
    # nbpoint, bins, v = plt.hist(point, bins=np.arange(-10, 10, 1), normed=True, edgecolor='k', color=None, alpha=0.3)
    # plt.xlabel('Max. LSTA[9x9km box] | 120km East of cell / random point in same LAT ')
    # plt.ylabel('Frequency')
    # stri = (np.sum(point >= np.percentile(all, 90)) / point.size * 100).round(2)
    # plt.title('0-4UTC:' + str(stri) + '% of Cells occur in warmest quartile')
    #
    #
    # plt.figure()
    # ngall, gbins, v = plt.hist(allgrad, bins=np.arange(-10,10,1), normed=True, edgecolor='k', color=None, alpha=0.3)
    # ngpoint, gbins, v = plt.hist(pointgrad, bins=np.arange(-10,10,1), normed=True, edgecolor='k', color=None, alpha=0.3)
    # p = 90
    # prob = np.sum(cinput > np.percentile(rinput, p)) / cinput.size
    # print(prob)
    # print( (prob - (100-p)*0.01) / ((100-p)*0.01)) # percentage of cells in warmest 25% of LSTA
    #
    # p = 10
    # prob = np.sum(cinput < np.percentile(rinput, p)) / cinput.size
    # print(prob)
    # print( (prob - p*0.01) / (p*0.01)) # percentage of cells in warmest 25% of LSTA

    # ((np.sum(point >= np.percentile(all, 75)) / point.size) -0.25) / 0.25



    return nball, nbpoint, bins


def plot_diurn():

    percmmax = []
    percmmin = []
    nbmax = []
    nbmin = []
    err90_up = []
    err90_low = []
    err10_up = []
    err10_low = []


    rrange = [15,16,17,18,19,20,21,22,23,0,1,2,3,4,5, 6,7,8,9,10]

    for h in rrange:
        path = cnst.network_data + 'figs/LSTA/corrected_LSTA/new/wavelet_coefficients'
        dic = pkl.load(open(path + "/LSTA_histograms_" + str(h).zfill(2) + ".p", "rb"))
        print('Open')

        for k in dic.keys():
            coll = []
            for ll in dic[k]:
                coll.extend(ll)
            dic[k] = coll

        cinput = np.array(dic['e100'])
        rinput = np.array(dic['re100'])
        point = cinput[np.isfinite(cinput)]
        all = rinput[np.isfinite(rinput)]


        p = 90
        pprob = np.sum(point > np.percentile(all, p))
        prob = pprob / point.size

        percmax = (prob - (100-p)*0.01) / ((100-p)*0.01) *100 # percentage of cells in warmest 25% of LSTA
        percmmax.append(percmax)
        nbmax.append(point > np.percentile(all, p))

        low90, upp90 = proportion_confint(pprob, point.size)

        err90_up.append( ((upp90 - (100-p)*0.01) / ((100-p)*0.01) *100) - percmax)
        err90_low.append(percmax -((low90 - (100 - p) * 0.01) / ((100 - p) * 0.01) * 100) )

        p = 10
        pprob = np.sum(point < np.percentile(all, p))
        prob = pprob / point.size
        percmin =  (prob - p*0.01) / (p*0.01) *100 # percentage of cells in warmest 25% of LSTA

        percmmin.append(percmin)
        nbmin.append(len(point))
        low10, upp10 = proportion_confint(pprob, point.size)

        err10_up.append((upp10 - p*0.01) / (p*0.01) *100 - percmin)
        err10_low.append( percmin - (low10 - p*0.01) / (p*0.01) *100 )

    f = plt.figure(figsize=(9,5), dpi=200)
    ax = f.add_subplot(111)
    ax.bar(np.arange(0,len(rrange)), percmmin,  label='10th centile',yerr=np.vstack((err10_up, err10_low)), edgecolor='k', alpha=0.5) #
    ax.bar(np.arange(0, len(rrange)), percmmax, label='90th centile', yerr=np.vstack((err90_up, err90_low)), edgecolor='k')
    ax.set_xticks(np.arange(0, len(rrange)))
    ax.set_xticklabels(rrange)

    ax.set_xlabel('Hour')
    #plt.axvline(16, color='k')

    #ax.set_xlim(-2,25)
    #plt.ylabel('NbCells|LSTA gt/lt perc. / NbCells')
    plt.ylabel('Difference in probability (%)')
    plt.legend()

    ax1 = ax.twiny()
    ax1.bar(np.arange(0, len(rrange)), percmmin, label='10th centile', yerr=np.vstack((err10_up, err10_low)), edgecolor='k', alpha=0.5)
    ax1.bar(np.arange(0, len(rrange)), percmmax, label='90th centile', yerr=np.vstack((err90_up, err90_low)), edgecolor='k')
    ax1.set_xticks(np.arange(0,len(rrange)))
    ax1.set_xticklabels(nbmin, rotation=45)
    ax1.set_xlabel('Number of convective cores')

    plt.tight_layout()
    #plt.annotate('a)', xy=(0.04, 0.94), xytext=(0, 4), size=15, xycoords=('figure fraction', 'figure fraction'),
    #             textcoords='offset points')  # transform=ax.transAxes,
    plt.savefig(path + '/paper/LSTA_core_probability_100kmScale_upstream.png')


def plot_diurn_triple():


    loop = [('c30', 'r30','Co-located, 30km length scale'), ('e100', 're100', '150km upstream, 100km length scale'), ('s100','rs100', '150km south, 100km length scale')]
    f = plt.figure(figsize=(5,10), dpi=200)
    for ids, input in enumerate(loop):
        rrange = [15,16,17,18,19,20,21,22,23,0,1,2,3,4,5]#, 6,7,8,9,10]
        percmmax = []
        percmmin = []
        nbmax = []
        nbmin = []
        err90_up = []
        err90_low = []
        err10_up = []
        err10_low = []

        for h in rrange:
            path = cnst.network_data + 'figs/LSTA/corrected_LSTA/new/wavelet_coefficients'
            dic = pkl.load(open(path + "/LSTA_histograms_" + str(h).zfill(2) + ".p", "rb"))
            print('Open')

            for k in dic.keys():
                coll = []
                for ll in dic[k]:
                    coll.extend(ll)
                dic[k] = coll

            cinput = np.array(dic[input[0]])
            rinput = np.array(dic[input[1]])
            point = cinput[np.isfinite(cinput)]
            all = rinput[np.isfinite(rinput)]


            p = 90
            pprob = np.sum(point > np.percentile(all, p))
            prob = pprob / point.size

            percmax = (prob - (100-p)*0.01) / ((100-p)*0.01) *100 # percentage of cells in warmest 25% of LSTA
            percmmax.append(percmax)
            nbmax.append(point > np.percentile(all, p))

            low90, upp90 = proportion_confint(pprob, point.size)

            err90_up.append( ((upp90 - (100-p)*0.01) / ((100-p)*0.01) *100) - percmax)
            err90_low.append(percmax -((low90 - (100 - p) * 0.01) / ((100 - p) * 0.01) * 100) )

            p = 10
            pprob = np.sum(point < np.percentile(all, p))
            prob = pprob / point.size
            percmin =  (prob - p*0.01) / (p*0.01) *100 # percentage of cells in warmest 25% of LSTA

            percmmin.append(percmin)
            nbmin.append(len(point))
            low10, upp10 = proportion_confint(pprob, point.size)

            err10_up.append((upp10 - p*0.01) / (p*0.01) *100 - percmin)
            err10_low.append( percmin - (low10 - p*0.01) / (p*0.01) *100 )

        ax = f.add_subplot(3,1,ids+1)
        ax.bar(np.arange(0,len(rrange)), percmmin,  label='10th centile',yerr=np.vstack((err10_up, err10_low)), edgecolor='k', alpha=0.5) #
        ax.bar(np.arange(0, len(rrange)), percmmax, label='90th centile', yerr=np.vstack((err90_up, err90_low)), edgecolor='k')
        ax.set_xticks(np.arange(0, len(rrange)))
        ax.set_xticklabels(rrange)

        ax.set_xlabel('Hour')

        plt.ylabel('Difference in probability (%)')
        plt.legend()

        ax1 = ax.twiny()
        ax1.bar(np.arange(0, len(rrange)), percmmin, label='10th centile', yerr=np.vstack((err10_up, err10_low)), edgecolor='k', alpha=0.5)
        ax1.bar(np.arange(0, len(rrange)), percmmax, label='90th centile', yerr=np.vstack((err90_up, err90_low)), edgecolor='k')
        ax1.set_xticks(np.arange(0,len(rrange)))
        ax1.set_xticklabels(nbmin, rotation=45)
        ax1.set_xlabel('Number of convective cores')

        plt.title(input[2])

    plt.tight_layout()
    #plt.annotate('a)', xy=(0.04, 0.94), xytext=(0, 4), size=15, xycoords=('figure fraction', 'figure fraction'),
    #             textcoords='offset points')  # transform=ax.transAxes,
    plt.savefig(path + '/paper/LSTA_core_probability_triple.png')