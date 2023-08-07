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

    pool = multiprocessing.Pool(processes=1)
    h = [15,16,17,18,19,20,21,22,23,0,1,2,3,4,5]
    res = pool.map(composite, h)
    pool.close()



def composite(hour):
    #hour = 18
    path = cnst.network_data + 'figs/LSTA/corrected_LSTA/new/wavelet_coefficients'
    dic = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_-60.p", "rb"))

    print('Doing with hour ', hour)

    all = []
    point = []

    scales = dic['scales']
    nbcores = dic['nbcores']

    del dic['scales']
    del dic['nbcores']
    del dic['nbrcores']
    del dic['kernel']
    del dic['lsta']
    del dic['cnt']

    dicout = {'30kmCore' : [],
              '30kmRand' : [],
              '100kmSouthCore' : [],
              '100kmSouthRand' : [],
              '100kmEastCore' : [],
              '100kmEastRand' : [],
              'scales' : [],
              'sumCore30': [],
              'sumRand30' : [],
              }

    pos = np.where(np.isclose(scales, 40, atol=4))
    spos = int(pos[0])
    dicout['scales'].append(scales[spos])
    cnt = []
    dist = 100
    for dat0, dat1 in zip(dic['SN-pos'][0], dic['WE-pos'][0]):
        print('30 loop', scales[spos])

        #ipdb.set_trace()
        arr = dat0[spos, dist-5:dist+5] #(dat0[spos, 96:105]).tolist()+(dat1[spos, 96:105]).tolist()
        means = max(np.nanmax(arr), np.nanmin(arr), key=abs)
        if np.isfinite(means):
            dicout['30kmCore'].append(means) # 24km across
            cnt.append(1)
        else:
            cnt.append(0)
        # if l == 'SN-pos':
        #     f = plt.figure()
        #     pplot = dat[spos,:]
        #     pplot[np.isnan(pplot)] = 0
        #     plt.plot(np.arange(0,201,1),pplot)
        #     plt.title(str(l))

    for dat0, dat1 in zip(dic['SN-pos'][1], dic['WE-pos'][1]):
        arr = dat0[spos, dist-5:dist+5] #(dat0[spos, 96:105]).tolist()+(dat1[spos, 96:105]).tolist()
        meansr = max(np.nanmax(arr), np.nanmin(arr), key=abs) #max(np.nanmax(dat[spos, 96:105]), np.nanmin(dat[spos, 96:105]), key=abs)
        if np.isfinite(meansr):
            dicout['30kmRand'].append(meansr)  # 24km across
    #return

    pos = np.where(np.isclose(scales, 109, atol=2))
    spos = int(pos[0])
    dicout['scales'].append(scales[spos])

    centre = dist-67
    for dat, cn in zip(dic['SN-pos'][0],cnt):
            if cn == 0:
                continue
            means = max(np.nanmax(dat[spos, centre-20:centre+20]), np.nanmin(dat[spos, centre-20:centre+20]), key=abs)
            if np.isfinite(means):
                dicout['100kmSouthCore'].append(means)  #24km across
    for dat, cn  in zip(dic['SN-pos'][1],cnt):
            if cn == 0:
                continue
            meansr = max(np.nanmax(dat[spos, centre-20:centre+20]), np.nanmin(dat[spos, centre-20:centre+20]), key=abs)
            if np.isfinite(meansr):
                dicout['100kmSouthRand'].append(meansr)  #24km across
    centre = dist+51
    for dat, cn  in zip(dic['WE-pos'][0],cnt):
            if cn == 0:
                continue
            means = max(np.nanmax(dat[spos, centre-20:centre+20]), np.nanmin(dat[spos, centre-20:centre+20]), key=abs)
            if np.isfinite(means):
                dicout['100kmEastCore'].append(means) #24km across
    for dat, cn  in zip(dic['WE-pos'][1],cnt):
            if cn == 0:
                continue
            meansr = max(np.nanmax(dat[spos, centre-20:centre+20]), np.nanmin(dat[spos, centre-20:centre+20]), key=abs) #at 450km, -150km upstream
            if np.isfinite(meansr):
                 dicout['100kmEastRand'].append(meansr)  #24km across

    pkl.dump(dicout, open(path+"wcoeff_histograms_"+str(hour).zfill(2)+"_randomTest.p", "wb"))


def plot(hour):
    path = cnst.network_data + 'figs/LSTA/corrected_LSTA/new/wavelet_coefficients'
    dic = pkl.load(open(path+"wcoeff_histograms_"+str(hour).zfill(2)+"_randomTest.p", "rb"))

    c30 = np.array(dic['30kmCore'])
    r30 = np.array(dic['30kmRand'])

    w100 = np.array(dic['100kmSouthCore'])
    wr100 = np.array(dic['100kmSouthRand'])

    d100 = np.array(dic['100kmEastCore'])
    dr100 = np.array(dic['100kmEastRand'])

    cinput = d100
    rinput = dr100


    f = plt.figure()
    ax = f.add_subplot(121)

    nball, bins,v = plt.hist(cinput, bins=np.arange(-7,7,0.25), normed=True, edgecolor='k', color=None, alpha=0.3)
    nbpoint, bins, v = plt.hist(rinput, bins=np.arange(-7,7,0.25), normed=True, edgecolor='k', color=None, alpha=0.3)
    plt.xlabel('Local wavelet coefficient, 30km')
    plt.ylabel('Frequency')
    stri = (np.sum(cinput >= np.percentile(rinput, 99)) / cinput.size * 100).round(2)
    plt.title('18-19UTC:'+str(stri)+'% of Cells occur in warmest decile')

    stri = (np.sum(cinput >= np.percentile(rinput, 50)) / cinput.size * 100).round(2)
    print('18-19UTC:' + str(stri) + '% of Cells occur in warmest half')


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
    p = 90
    prob = np.sum(c30 > np.percentile(r30, p)) / r30.size
    print(prob)
    print( (prob - (100-p)*0.01) / ((100-p)*0.01)) # percentage of cells in warmest 25% of LSTA

    p = 10
    prob = np.sum(c30 < np.percentile(r30, p)) / r30.size
    print(prob)
    print( (prob - p*0.01) / (p*0.01)) # percentage of cells in warmest 25% of LSTA

    # ((np.sum(point >= np.percentile(all, 75)) / point.size) -0.25) / 0.25



    #return nball, nbpoint, bins


def plot_diurn():

    percmmax = []
    percmmin = []
    nbmax = []
    nbmin = []
    err90_up = []
    err90_low = []
    err10_up = []
    err10_low = []


    rrange = [15,16,17,18,19,20,21,22,23,0,1,2,3,4,5]

    for h in rrange:
        path = cnst.network_data + 'figs/LSTA/corrected_LSTA/new/wavelet_coefficients'
        dic = pkl.load(open(path+"wcoeff_histograms_"+str(h).zfill(2)+"_randomTest.p", "rb"))
        print('Open')

        # point = np.array(dic['100kmEastCore'])
        # all = np.array(dic['100kmEastRand'])

        # point = np.array(dic['30kmCore'])
        # all = np.array(dic['30kmRand'])

        point = np.array(dic['30kmCore'])
        all = np.array(dic['30kmRand'])

        p = 75
        pprob = np.sum(point > np.percentile(all, p))
        prob = pprob / point.size

        percmax = (prob - (100-p)*0.01) / ((100-p)*0.01) *100 # percentage of cells in warmest 25% of LSTA
        percmmax.append(percmax)
        nbmax.append(point > np.percentile(all, p))

        low90, upp90 = proportion_confint(pprob, point.size)

        err90_up.append( ((upp90 - (100-p)*0.01) / ((100-p)*0.01) *100) - percmax)
        err90_low.append(percmax -((low90 - (100 - p) * 0.01) / ((100 - p) * 0.01) * 100) )

        p = 25
        pprob = np.sum(point < np.percentile(all, p))
        prob = pprob / point.size
        percmin =  (prob - p*0.01) / (p*0.01) *100 # percentage of cells in warmest 25% of LSTA

        percmmin.append(percmin)
        nbmin.append(len(point))
        low10, upp10 = proportion_confint(pprob, point.size)

        err10_up.append((upp10 - p*0.01) / (p*0.01) *100 - percmin)
        err10_low.append( percmin - (low10 - p*0.01) / (p*0.01) *100 )

    f = plt.figure(figsize=(9,5))
    ax = f.add_subplot(111)
    ax.bar(np.arange(0,len(rrange)), percmmin,  label='10th centile',yerr=np.vstack((err10_up, err10_low)), edgecolor='k') #
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
    ax1.bar(np.arange(0, len(rrange)), percmmin, label='10th centile', yerr=np.vstack((err10_up, err10_low)), edgecolor='k')
    ax1.bar(np.arange(0, len(rrange)), percmmax, label='90th centile', yerr=np.vstack((err90_up, err90_low)), edgecolor='k')
    ax1.set_xticks(np.arange(0,len(rrange)))
    ax1.set_xticklabels(nbmin, rotation=45)
    ax1.set_xlabel('Number of convective cores')

    plt.tight_layout()
    plt.annotate('a)', xy=(0.04, 0.94), xytext=(0, 4), size=15, xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points')  # transform=ax.transAxes,
