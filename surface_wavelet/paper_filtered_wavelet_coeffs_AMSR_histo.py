# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 10:15:40 2016

@author: cornkle
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib
import multiprocessing
import ipdb
import pandas as pd
from utils import u_arrays, constants as cnst
import pickle as pkl
import collections
from utils import u_statistics as u_stat
from statsmodels.stats.proportion import proportion_confint
import salem
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import norm
from scipy.stats import ttest_ind

plt.style.use('seaborn-ticks')

matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)
matplotlib.rcParams['hatch.linewidth'] = 0.1

def run_hours():

    l = [15,16,17,18,19,20,21,22,23,0,1,2,3,4,5,6,7] #15,16,17,18,19
    for ll in l:
        composite(ll)

def composite(h):

    key = 'NEWTRACKING'

    msgopen = pd.read_csv(
        cnst.network_data + 'figs/LSTA/corrected_LSTA/new/ERA5/core_txt/init_merged/cores_gt15000km2_table_AMSRE_tracking_' + str(
            h) + '_init.csv', na_values=[-999, -99])


    msg = pd.DataFrame.from_dict(msgopen)# &  &
    msg['date'] = pd.to_datetime(msg[['year','month','day']])
    print('Start core number ', len(msg))

    msgopen = msg

    #basic filter
    msgopen = msgopen[(msgopen['lat']>9.5) & (msgopen['lat']<20.5) & (msgopen['topo']<=450) & (msgopen['dtime']<=2)]
    #propagation filter
    msgopen = msgopen[(msgopen['xdiff']>=100) | (msgopen['initTime'] <= 2.5)]

    #SM_filter
    msgopen = msgopen[np.isfinite(msgopen['SMmean0'])]# & np.isfinite(msgopen['SMmean-1'])]
    #wetness_filter
    # msgopen = msgopen[np.isfinite(msgopen['SMmean0']) & np.isfinite(msgopen['SMmean-1'])]
    # #eraq_filter
    # msgopen = msgopen[(msgopen['ERAqmean'] >= 13.5)]#& (msgopen['ERAqmean'] <= 15.7)]
    # #dry_filter
    # #msgopen = msgopen[(msgopen['SMmean0']<=-6.7)&(msgopen['SMmean-1']<=-0.01)] #294 cases, with q 312
    # # #wet_filter
    # msgopen = msgopen[(msgopen['SMmean0']>=0.3) & (msgopen['SMmean-1']>=-0.01)] #295 cases, with q 318

    msgin = msgopen


    print('Number of cores', len(msgin))
    msgin.sort_values(by='date')


    msgy = msgin

    chunk, chunk_ind, chunk_count = np.unique(msgy.date, return_index=True, return_counts=True)

    chunks = [msgy.loc[msgy.index[ci:ci + cc]] for ci, cc in zip(chunk_ind, chunk_count)] # daily chunks

    # res = []
    # for m in chunks[0:2]:
    #     print(m['SMmean0'])
    #     out = file_loop(m)
    #     res.append(out)

    #return
    pool = multiprocessing.Pool(processes=5)

    res = pool.map(file_loop, chunks)
    pool.close()

    print('Returned from parallel')

    res = [x for x in res if x is not None]

    amsrk30 = []
    amsre100 = []

    ramsrk30 = []
    ramsre100 = []

    cores = 0
    for r in res:

        amsrk30.extend(r[0])
        amsre100.extend(r[1])


        ramsrk30.extend(r[2])
        ramsre100.extend(r[3])

        cores += r[4]

    dic = collections.OrderedDict([
                                   ('amsr' , [amsrk30, amsre100]),
                                    ('ramsr', [ramsrk30, ramsre100]),
                                   ('cores', cores)])


    outpath = cnst.network_data + '/figs/LSTA/corrected_LSTA/new/wavelet_coefficients/'
    pkl.dump(dic, open(outpath + "/AMSR_histograms_" + str(h).zfill(2) + "_SMFINITE_150kcircle_minusMean_bigDomain_topo_"+key+".p", "wb"))
    print('Save file written!')




def cut_kernel_lsta(xpos, ypos, arr):

    dist=200

    kernel = u_arrays.cut_kernel(arr,xpos, ypos,dist)

    if (np.sum(np.isfinite(kernel)) < 0.01 * kernel.size):
        return

    kmean = kernel - np.nanmean(kernel)

    if kernel.shape != (2*dist+1, 2*dist+1):
        print('WRONG SHAPE')
        return
    ycirc30, xcirc30 = u_arrays.draw_circle(dist+1, dist+1,6) # 15km radius
    k30 = np.nanmean(kmean[ycirc30, xcirc30])

    ycirc100e, xcirc100e = u_arrays.draw_circle(dist+51, dist+1, 17)  # at - 150km, draw 50km radius circle
    data = kmean[ycirc100e,xcirc100e]
   # data = kernel[dist - 10:dist + 10, dist+6:dist + 67]

    e100 = np.nanmean(data)

    if np.isnan(e100):
        return None

    if (np.sum(np.isfinite(data)) / data.size) < 0.01:
        print('Too small valid area')
        return

    return k30, e100



def file_loop(df):

    date = df['date'].iloc[0]
    hour = df['hour'].iloc[0]
    print('Doing day: ', date)

    storm_date = date

    dayd = pd.Timedelta('1 days')

    if (hour) <= 13:
        print('Nighttime')
        daybefore = storm_date - dayd
    else:
        print('Daytime')
        daybefore = storm_date

    fdate = str(daybefore.year) + str(daybefore.month).zfill(2) + str(daybefore.day).zfill(2)

    amsre = xr.open_dataset(cnst.AMSRE_ANO_DAY + 'sma_' + fdate + '.nc')
    amsre = amsre.sel(time=str(daybefore.year)+'-'+str(daybefore.month)+'-'+str(daybefore.day))
    amsre = amsre.sel(lon=slice(-12.5, 12.5), lat=slice(7, 23))
    print('Doing '+ 'AMSR_' + str(daybefore.year) + str(daybefore.month).zfill(2) + str(
        daybefore.day).zfill(2) + '.nc')

    amsr_da = amsre['SM'].squeeze()


    #topo = xr.open_dataset(cnst.LSTA_TOPO)
    topo = xr.open_dataset(cnst.WA_TOPO_3KM)
    topo = topo.sel(lon=slice(-12, 12), lat=slice(7, 23))

    ttopo = topo['h']
    grad = np.gradient(ttopo.values)
    gradsum = abs(grad[0]) + abs(grad[1])

    # if (np.sum(np.isfinite(amsr_da)) / amsr_da.size) < 0.01:
    #     print('Not enough valid')
    #     return None

    try:
        amsr_da = topo.salem.transform(amsr_da)
    except RuntimeError:
        print('amsr_da on LSTA interpolation problem')
        return None

    amsr_da.values[ttopo.values >= 450] = np.nan
    amsr_da.values[gradsum > 30] = np.nan

    del topo

    amsrk30 = []
    amsre100 = []

    ramsrk30 = []
    ramsre100 = []


    ###############################Blob loop
    cores = 0

    for dids, dit in df.iterrows():

        lat = dit.lat
        lon = dit.lon

        try:
            point = amsr_da.sel(lat=lat, lon=lon, method='nearest', tolerance=0.03)
        except KeyError:
            print('Nearest point finding error')
            continue

        plat = point['lat'].values
        plon = point['lon'].values

        xpos = np.where(amsr_da['lon'].values == plon)
        xpos = int(xpos[0])
        ypos = np.where(amsr_da['lat'].values == plat)
        ypos = int(ypos[0])

        try:
            ak30, ae100 = cut_kernel_lsta(xpos, ypos, amsr_da.values)
        except TypeError:
            print('AMSR kernel error')
            continue


        amsrk30.append(ak30)
        amsre100.append(ae100)

        cores += 1

        ##### random

        y = ypos
        x = xpos

        rdist = 50
        randy50 = [y - rdist, y - rdist, y - rdist, y, y, y + rdist, y + rdist, y + rdist]
        randx50 = [x - rdist, x, x + rdist, x - rdist, x + rdist, x - rdist, x, x + rdist]
        randy50_100 = [y - rdist, y - rdist, y, y, y + rdist, y + rdist]

        rdist = 100
        randx100 = [x - rdist, x + rdist, x - rdist, x + rdist, x - rdist, x + rdist]

        rdist = 150
        randx150 = [x - rdist, x + rdist, x - rdist, x + rdist, x - rdist, x + rdist]

        randy = np.array(randy50 + randy50_100 + randy50_100)
        randx = np.array(randx50 + randx100 + randx150)

        for ry, rx in zip(randy,randx):

            if ry < 0:
                continue
            if ry > amsr_da.shape[0] - 1:
                continue

            if rx < 0:
                continue
            if rx > amsr_da.shape[1] - 1:
                continue

            try:
                lat = amsr_da['lat'][ry]
            except IndexError:
                ipdb.set_trace()
            try:
                lon = amsr_da['lon'][rx]
            except IndexError:
                ipdb.set_trace()

            point = amsr_da.sel(lat=lat, lon=lon, method='nearest')
            plat = point['lat'].values
            plon = point['lon'].values

            xpos = np.where(amsr_da['lon'].values == plon)
            xpos = int(xpos[0])
            ypos = np.where(amsr_da['lat'].values == plat)
            ypos = int(ypos[0])

            try:
                arc30, arce100 = cut_kernel_lsta(xpos, ypos, amsr_da.values)
            except TypeError:
                continue

            ramsrk30.append(arc30)
            ramsre100.append(arce100)

    del amsr_da

    if (len(amsrk30) == 0) :
        return None


    print('Returning with kernel, success!!')

    return (amsrk30, amsre100,
            ramsrk30, ramsre100, cores)

if __name__ == "__main__":
    run_hours()



def plot_diurn_double_relative():


    loop = [('c30', 'r30','Co-located, 30km length scale'), ('e100', 're100', '100km upstream, 100km length scale')]

    f = plt.figure(figsize=(5,8), dpi=200)
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
            path = cnst.network_data + 'figs/LSTA/corrected_LSTA/new/wavelet_coefficients/'
            #dic = pkl.load(open(path + "LSTA_histograms_AMSRE_"+str(h).zfill(2)+"SlotFilter_+150km_validCheck.p", "rb"))
            dic = pkl.load(
                open(path + "AMSR_histograms_" + str(h).zfill(2) + "_SMFINITE.p", "rb"))
            print('Open')


            cinput = np.array(dic['amsr'][ids])
            rinput = np.array(dic['ramsr'][ids])


            point = cinput[np.isfinite(cinput)]
            all = rinput[np.isfinite(rinput)]

            p = 75
            pprob = np.sum(point > np.percentile(all, p))
            prob = pprob / point.size

            percmax = prob  *100 # percentage of cells in warmest 25% of LSTA
            percmmax.append(percmax)
            nbmax.append(point > np.percentile(all, p))

            low90, upp90 = proportion_confint(pprob, point.size)

            err90_up.append((upp90 *100) - percmax)
            err90_low.append(percmax -(low90 * 100) )

            p = 25
            pprob = np.sum(point < np.percentile(all, p))
            prob = pprob / point.size
            #ipdb.set_trace()
            percmin =  prob  *100 # percentage of cells in warmest 25% of LSTA

            print(h, '10prob', prob)
            print(h, 'percent increase', (prob-0.1)/0.1*100)

            percmmin.append(percmin)
            nbmin.append(len(point))
            low10, upp10 = proportion_confint(pprob, point.size)

            err10_up.append((upp10 *100) - percmin)
            err10_low.append( percmin - (low10 *100 ))

        ax = f.add_subplot(2,1,ids+1)
        ax.bar(np.arange(0,len(rrange)), np.array(percmmin)-25,  label='25th centile',yerr=np.vstack((err10_up, err10_low)), edgecolor='k', color='darkorange') #
        ax.bar(np.arange(0, len(rrange)), np.array(percmmax)-25, label='75th centile', yerr=np.vstack((err90_up, err90_low)), edgecolor='k',color='powderblue')
        ax.set_xticks(np.arange(0, len(rrange)))
        ax.set_xticklabels(rrange)

        lw = 0.5


        ax.set_xlabel('Hour')

        plt.ylabel('Probability (%)')
        plt.legend()

        ax1 = ax.twiny()
        ax1.bar(np.arange(0, len(rrange)), np.array(percmmin)-25, label='25th centile', yerr=np.vstack((err10_up, err10_low)), edgecolor='k', alpha=0.8, color='darkorange')
        ax1.bar(np.arange(0, len(rrange)), np.array(percmmax)-25, label='75th centile', yerr=np.vstack((err90_up, err90_low)), edgecolor='k',color='powderblue')
        ax1.set_xticks(np.arange(0,len(rrange)))
        ax1.set_xticklabels(nbmin, rotation=45)
        ax1.set_xlabel('Number of convective cores')

        ax.set_ylim(0 - 25, 60 - 25)
        locs, ylabels = plt.yticks()
        # print(ids, locs)
        ax.set_yticklabels(locs+25)


        plt.title(input[2])

        # plt.axhline(y=-50, linewidth=lw, color='k', linestyle='dashed')
        #plt.axhline(y=25, linewidth=lw, color='k', linestyle='dashed')
        plt.axhline(y=0, linewidth=3, color='k', linestyle='solid')
        #plt.axhline(y=-25, linewidth=lw, color='k', linestyle='dashed')
    plt.tight_layout()

    #plt.annotate('a)', xy=(0.04, 0.94), xytext=(0, 4), size=15, xycoords=('figure fraction', 'figure fraction'),
    #             textcoords='offset points')  # transform=ax.transAxes,
    plt.show()
    plt.savefig(path + '/initVSprop/AMSRE_core_probability_RELATIVE_100km_SMFINITE.png')

def plot_diurn_double_relative_double():


    loop = [('e100', 're100', '100km upstream, 100km length scale')] #('c30', 'r30','Co-located, 30km length scale'),

    f = plt.figure(figsize=(11,4), dpi=200)

    ax = plt.axes([0.47,0.12,0.5,0.7])

    for ids, input in enumerate(loop):
        ids = ids+1
        rrange = [15,16,17,18,19,20,21,22,23,0,1,2,3,4,5]#, 6,7,8,9,10]
        percmmax = []
        percmmin = []
        nbmax = []
        nbmin = []
        err90_up = []
        err90_low = []
        err10_up = []
        err10_low = []

        name = '_SMFINITE_150kcircle_minusMean_bigDomain_topo_NEWTRACKING'

        for h in rrange:
            path = cnst.network_data + 'figs/LSTA/corrected_LSTA/new/wavelet_coefficients/'
            #dic = pkl.load(open(path + "LSTA_histograms_AMSRE_"+str(h).zfill(2)+"SlotFilter_+150km_validCheck.p", "rb"))
            # dic = pkl.load(
            #     open(path + "AMSR_histograms_" + str(h).zfill(2) + "_SMFINITE.p", "rb"))
            #AMSR_histograms_03_SMFINITE_150kcircle_minusMean_bigDomain_topo
            # dic = pkl.load(
            #     open(path + "AMSR_histograms_" + str(h).zfill(2) + "_SMFINITE.p", "rb"))
            dic = pkl.load(
                open(path + "AMSR_histograms_" + str(h).zfill(2) + name + ".p", "rb"))
            print('Open')


            cinput = np.array(dic['amsr'][ids])
            rinput = np.array(dic['ramsr'][ids])
            point = cinput[np.isfinite(cinput)]
            all = rinput[np.isfinite(rinput)]

            p = 75
            pprob = np.sum(point > np.percentile(all, p))
            prob = pprob / point.size

            percmax = prob  *100 # percentage of cells in warmest 25% of LSTA
            percmmax.append(percmax)
            nbmax.append(point > np.percentile(all, p))

            low90, upp90 = proportion_confint(pprob, point.size)

            err90_up.append((upp90 *100) - percmax)
            err90_low.append(percmax -(low90 * 100) )

            p = 25
            pprob = np.sum(point < np.percentile(all, p))
            prob = pprob / point.size
            #ipdb.set_trace()
            percmin =  prob  *100 # percentage of cells in warmest 25% of LSTA

            print(h, '10prob', prob)
            print(h, 'percent increase', (prob-0.1)/0.1*100)

            percmmin.append(percmin)
            nbmin.append(len(point))
            low10, upp10 = proportion_confint(pprob, point.size)

            err10_up.append((upp10 *100) - percmin)
            err10_low.append( percmin - (low10 *100 ))


        ax.bar(np.arange(0,len(rrange)), np.array(percmmin)-25,  label='25th centile',yerr=np.vstack((err10_up, err10_low)), edgecolor='k', color='coral') #
        ax.bar(np.arange(0, len(rrange)), np.array(percmmax)-25, label='75th centile', yerr=np.vstack((err90_up, err90_low)), edgecolor='k',color='powderblue')
        ax.set_xticks(np.arange(0, len(rrange)))
        ax.set_xticklabels(rrange)

        lw = 0.5


        ax.set_xlabel('Hour')

        plt.ylabel('Probability (%)')
        plt.legend()

        ax1 = ax.twiny()
        ax1.bar(np.arange(0, len(rrange)), np.array(percmmin)-25, label='25th centile', yerr=np.vstack((err10_up, err10_low)), edgecolor='k', alpha=0.6, color='coral')
        ax1.bar(np.arange(0, len(rrange)), np.array(percmmax)-25, label='75th centile', yerr=np.vstack((err90_up, err90_low)), edgecolor='k',color='powderblue')
        ax1.set_xticks(np.arange(0,len(rrange)))
        ax1.set_xticklabels(nbmin, rotation=45)
        ax1.set_xlabel('Number of convective cores')

        ax.set_ylim(0 - 25, 60 - 25)
        locs, ylabels = plt.yticks()
        # print(ids, locs)
        ax.set_yticklabels(locs+25)


        #plt.title(input[2])

        # plt.axhline(y=-50, linewidth=lw, color='k', linestyle='dashed')
        #plt.axhline(y=25, linewidth=lw, color='k', linestyle='dashed')
        plt.axhline(y=0, linewidth=3, color='k', linestyle='solid')

    ax1 = plt.axes([0.07, 0.12, 0.30, 0.8])
    hour=17
    dic = pkl.load(
        open(path + "AMSR_histograms_" + str(hour).zfill(2) + "_SMFINITE.p", "rb"))

    cinput = np.array(dic['amsr'][1])
    rinput = np.array(dic['ramsr'][1])
    point = cinput[np.isfinite(cinput)]
    all = rinput[np.isfinite(rinput)]

    print('NUMBEROFCORE', dic['cores'])

    nbpoint, pointcount, bins = u_stat.histo_frequency(point, bins=np.arange(-15, 15, 1))
    nball, allcount, bins = u_stat.histo_frequency(all, bins=np.arange(-15, 15, 1))
    print(bins)
    bin_centre = bins[0:-1] + ((bins[1::] - bins[0:-1]) / 2)
    bin_edge = bins[0:-1]
    width = bins[1::] - bins[0:-1]


    randomp = np.zeros_like(bin_centre)
    loop = 1000
    for num in range(loop):
        print('Loop', num)
        rand = np.random.choice(all, size=len(point), replace=False)
        #ipdb.set_trace()
        (mu_random, sigma_random) = norm.fit(rand)
        randoms = norm.pdf(bin_centre, mu_random, sigma_random)

        randomp += randoms
    randomp = randomp/loop

    # best fit of data
    (mu_core, sigma_core) = norm.fit(point)
    #(mu_random, sigma_random) = norm.fit(all)

    corep = norm.pdf(bin_centre,mu_core,sigma_core)
    #randomp = norm.pdf(bins,mu_random,sigma_random)

    tstat, pvalue = ttest_ind(point, all)
    ipdb.set_trace()


    ax1.bar(bin_edge, nbpoint, label='Cores', edgecolor='k', alpha=0.95, align='edge', width=width, facecolor='coral', zorder=9)
    ax1.bar(bin_edge, nball, label='Random', edgecolor='k', alpha=0.7, align='edge', width=width, facecolor='grey', zorder=10)

    ax1.plot(bin_centre,corep, color='red',zorder=12)
    ax1.plot(bin_centre, randomp, color='k',zorder=11)
    plt.ylabel('Frequency')
    plt.xlabel('Soil moisture anomaly (%)')
    plt.axvline(0, color='k', linestyle='solid', linewidth=2,zorder=13)
    #ipdb.set_trace()
    # plt.axvline(norm.ppf(0.25,mu_random,sigma_random), color='k', linestyle='dashed')
    # plt.axvline(norm.ppf(0.75,mu_random,sigma_random), color='k', linestyle='dashed')
    # ax1.axvspan(norm.ppf(0.00001,mu_random,sigma_random), norm.ppf(0.25,mu_random,sigma_random), alpha=0.5, color='grey', zorder=0)
    ax1.axvspan(norm.ppf(0.25, mu_random, sigma_random), norm.ppf(0.75, mu_random, sigma_random), alpha=0.5,
                color='grey', zorder=0)

    ax1.set_ylim(0,0.175)

    plt.axvline(norm.ppf(0.25, mu_random, sigma_random), color='k', linestyle='solid', linewidth=0.2, zorder=0)
    plt.axvline(norm.ppf(0.75, mu_random, sigma_random), color='k', linestyle='solid', linewidth=0.2, zorder=0)

    stri = (np.sum(cinput <= np.percentile(rinput, 25)) / cinput.size * 100).round(2)
    print(str(hour) + 'UTC:' + str(stri) + '% of Cells occur in warmest decile')

    stri = (np.sum(cinput <= np.percentile(rinput, 50)) / cinput.size * 100).round(2)
    print(str(hour) + 'UTC:' + str(stri) + '% of Cells occur in warmest half')
    plt.legend()


        #plt.axhline(y=-25, linewidth=lw, color='k', linestyle='dashed')
    plt.tight_layout()
    text = ['a', 'b']
    fs = 14
    plt.annotate(text[0], xy=(0.05, 0.91), xytext=(0, 4),xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points', fontweight='bold', fontname='Ubuntu', fontsize=fs)
    plt.annotate(text[1], xy=(0.44, 0.91), xytext=(0, 4), xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points', fontweight='bold', fontname='Ubuntu', fontsize=fs)
    plt.show()
    plt.savefig(path + '/initVSprop/AMSRE_core_probability_RELATIVE_' + name + '.png')

    plt.figure()
    plt.plot(np.arange(-15,15), norm.cdf(np.arange(-15,15),mu_random,sigma_random))

    print(norm.ppf(0.25,mu_random,sigma_random))


def plot(hour):
    path = cnst.network_data + 'figs/LSTA/corrected_LSTA/new/wavelet_coefficients'
    dic = pkl.load(open(path+"/AMSR_histograms_" + str(hour).zfill(2) + "_SMFINITE.p", "rb"))
    #dic = pkl.load(open(path + "/LSTA_histograms_AMSRE_" + str(hour).zfill(2) + "SlotFilter_+150km_validCheck.p", "rb"))

    cinput = np.array(dic['amsr'][1])
    rinput = np.array(dic['ramsr'][1])
    point = cinput[np.isfinite(cinput)]
    all = rinput[np.isfinite(rinput)]

    print('NUMBEROFCORE', dic['cores'])

    nbpoint, pointcount, bins = u_stat.histo_frequency(point, bins=np.arange(-15,15,1))
    nball, allcount, bins = u_stat.histo_frequency(all, bins=np.arange(-15, 15, 1))
    ipdb.set_trace()
    print(bins)
    bin_centre = bins[0:-1] + ((bins[1::] - bins[0:-1]) / 2)
    bin_edge = bins[0:-1]
    width = bins[1::] - bins[0:-1]

    f = plt.figure()
    ax = f.add_subplot(111)

    ax.bar(bin_edge, nbpoint, label='core', edgecolor='k', alpha=0.5, align='edge', width=width)
    ax.bar(bin_edge, nball, label='core', edgecolor='k', alpha=0.5, align='edge', width=width)
    plt.ylabel('Frequency')
    stri = (np.sum(cinput <= np.percentile(rinput, 25)) / cinput.size * 100).round(2)
    plt.title(str(hour)+'UTC:'+str(stri)+'% of Cells occur in warmest decile')

    stri = (np.sum(cinput <= np.percentile(rinput, 50)) / cinput.size * 100).round(2)
    print(str(hour)+'UTC:' + str(stri) + '% of Cells occur in warmest half')
