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
from utils import u_statistics as u_stat

def diurnal_loop():

    pool = multiprocessing.Pool(processes=4)
    h = [15,16,17,18,19,20,21,22,23,0,1,2,3,4,5,6,7,8,9,10,11,12,13]
    res = pool.map(plot, h)
    pool.close()



def plot(hour):
    path = cnst.network_data + 'figs/LSTA/corrected_LSTA/new/wavelet_coefficients'
    dic = pkl.load(open(path+"/AMSL_histograms_" + str(hour).zfill(2) + "_SMFINITE.p", "rb"))
    #dic = pkl.load(open(path + "/LSTA_histograms_AMSRE_" + str(hour).zfill(2) + "SlotFilter_+150km_validCheck.p", "rb"))

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

    print('NUMBEROFCORE', cinput.shape)

    nbpoint, pointcount, bins = u_stat.histo_frequency(cinput, bins=np.arange(-15,15,1))
    nball, allcount, bins = u_stat.histo_frequency(rinput, bins=np.arange(-15, 15, 1))
    print(bins)
    bin_centre = bins[0:-1] + ((bins[1::] - bins[0:-1]) / 2)
    bin_edge = bins[0:-1]
    width = bins[1::] - bins[0:-1]

    f = plt.figure()
    ax = f.add_subplot(121)

    # nball, bins,v = plt.hist(cinput, bins=np.arange(-7,7,0.5), normed=True, edgecolor='k', color=None, alpha=0.3)
    # nbpoint, bins, v = plt.hist(rinput, bins=np.arange(-7,7,0.5), normed=True, edgecolor='k', color=None, alpha=0.3)
    # plt.xlabel('Local wavelet coefficient, 30km')

    ax.bar(bin_edge, nbpoint, label='core', edgecolor='k', alpha=0.5, align='edge', width=width)
    ax.bar(bin_edge, nball, label='core', edgecolor='k', alpha=0.5, align='edge', width=width)
    plt.ylabel('Frequency')
    stri = (np.sum(cinput <= np.percentile(rinput, 10)) / cinput.size * 100).round(2)
    plt.title(str(hour)+'UTC:'+str(stri)+'% of Cells occur in warmest decile')

    stri = (np.sum(cinput <= np.percentile(rinput, 50)) / cinput.size * 100).round(2)
    print(str(hour)+'UTC:' + str(stri) + '% of Cells occur in warmest half')

    cumulative_random = np.cumsum(nball)
    cumulative = np.cumsum(nbpoint)

    ax = f.add_subplot(122)
    plt.title(str(hour))
    plt.plot(bin_centre,cumulative, label='cores')
    plt.plot(bin_centre, cumulative_random, label='random')
    plt.axvline(0,ymin=0, ymax=1, linestyle='dashed', color='k')
    plt.legend()


    plt.figure()
    plt.hist(rinput, range=(-15,15), bins=30)
    plt.axvline(0)
    plt.show()


def plot_double(hour):

    tags = [('c30', 'r30', 'Local - 30km scale'), ('e100', 're100', '100km upstream - 100km scale')]
    f = plt.figure(figsize=(9,6), dpi=200)
    for id,h in enumerate(hour):
        path = cnst.network_data + 'figs/LSTA/corrected_LSTA/new/wavelet_coefficients'
        dic = pkl.load(open(path+"/LSTA_histograms_AMSRE_"+str(h).zfill(2)+"_corrected_SouthBox.p", "rb"))
        #dic = pkl.load(open(path + "/LSTA_histograms_AMSRE_" + str(hour).zfill(2) + "SlotFilter_+150km_validCheck.p", "rb"))

        for k in dic.keys():
            coll = []
            for ll in dic[k]:
                coll.extend(ll)
            dic[k] = coll
        print(dic.keys())

        cinput = np.array(dic[tags[id][0]])
        rinput = np.array(dic[tags[id][1]])

        cinput = cinput[np.isfinite(cinput)]
        rinput = rinput[np.isfinite(rinput)]

        nbpoint, pointcount, bins = u_stat.histo_frequency(cinput, bins=np.arange(-15,15,1))
        nball, allcount, bins = u_stat.histo_frequency(rinput, bins=np.arange(-15, 15, 1))
        print(bins)
        bin_centre = bins[0:-1] + ((bins[1::] - bins[0:-1]) / 2)
        bin_edge = bins[0:-1]
        width = bins[1::] - bins[0:-1]


        ax = f.add_subplot(1,2,id+1)

        # nball, bins,v = plt.hist(cinput, bins=np.arange(-7,7,0.5), normed=True, edgecolor='k', color=None, alpha=0.3)
        # nbpoint, bins, v = plt.hist(rinput, bins=np.arange(-7,7,0.5), normed=True, edgecolor='k', color=None, alpha=0.3)
        # plt.xlabel('Local wavelet coefficient, 30km')

        ax.bar(bin_edge, nbpoint, label='core', edgecolor='k', alpha=0.5, align='edge', width=width)
        ax.bar(bin_edge, nball, label='core', edgecolor='k', alpha=0.5, align='edge', width=width)
        plt.ylabel('Frequency')
        stri = (np.sum(cinput <= np.percentile(rinput, 10)) / cinput.size * 100).round(2)
        plt.title(str(hour[0])+'UTC: '+tags[id][2])
        plt.axvline(x=0, color='k', linewidth=2)
        stri = (np.sum(cinput <= np.percentile(rinput, 50)) / cinput.size * 100).round(2)
        print(tags[id][2])

    plt.tight_layout()
    plt.savefig(path + '/initVSprop/AMSRE_histo_perHOUR_'+str(hour[0])+'.png')


def plot_double_relative(hour):

    tags = [('c30', 'r30', 'Local - 30km scale'), ('e100', 're100', '100km upstream - 100km scale')]
    f = plt.figure(figsize=(9,6), dpi=200)
    for id,h in enumerate(hour):
        path = cnst.network_data + 'figs/LSTA/corrected_LSTA/new/wavelet_coefficients'
        dic = pkl.load(open(path+"/LSTA_histograms_AMSRE_"+str(h).zfill(2)+"_corrected_SouthBox.p", "rb"))
        #dic = pkl.load(open(path + "/LSTA_histograms_AMSRE_" + str(hour).zfill(2) + "SlotFilter_+150km_validCheck.p", "rb"))

        for k in dic.keys():
            coll = []
            for ll in dic[k]:
                coll.extend(ll)
            dic[k] = coll
        print(dic.keys())

        cinput = np.array(dic[tags[id][0]])
        rinput = np.array(dic[tags[id][1]])

        cinput = cinput[np.isfinite(cinput)]
        rinput = rinput[np.isfinite(rinput)]

        nbpoint, pointcount, bins = u_stat.histo_frequency(cinput, bins=np.arange(-15,15,1))
        nball, allcount, bins = u_stat.histo_frequency(rinput, bins=np.arange(-15, 15, 1))
        print(bins)
        bin_centre = bins[0:-1] + ((bins[1::] - bins[0:-1]) / 2)
        bin_edge = bins[0:-1]
        width = bins[1::] - bins[0:-1]


        ax = f.add_subplot(1,2,id+1)

        # nball, bins,v = plt.hist(cinput, bins=np.arange(-7,7,0.5), normed=True, edgecolor='k', color=None, alpha=0.3)
        # nbpoint, bins, v = plt.hist(rinput, bins=np.arange(-7,7,0.5), normed=True, edgecolor='k', color=None, alpha=0.3)
        # plt.xlabel('Local wavelet coefficient, 30km')

        ax.bar(bin_edge, nbpoint, label='core', edgecolor='k', alpha=0.5, align='edge', width=width)
        ax.bar(bin_edge, nball, label='core', edgecolor='k', alpha=0.5, align='edge', width=width)
        plt.ylabel('Frequency')
        stri = (np.sum(cinput <= np.percentile(rinput, 10)) / cinput.size * 100).round(2)
        plt.title(str(hour[0])+'UTC: '+tags[id][2])
        plt.axvline(x=0, color='k', linewidth=2)
        stri = (np.sum(cinput <= np.percentile(rinput, 50)) / cinput.size * 100).round(2)
        print(tags[id][2])

    plt.tight_layout()
    plt.savefig(path + '/initVSprop/AMSRE_histo_perHOUR_'+str(hour[0])+'.png')



def plot_rand():

    f = plt.figure()
    ax = f.add_subplot(121)
    cinput = []
    rinput = []
    for hour in [16,17,18]:
        path = cnst.network_data + 'figs/LSTA/corrected_LSTA/new/wavelet_coefficients'
        dic = pkl.load(open(path+"/LSTA_histograms_AMSRE_"+str(hour).zfill(2)+"_corrected_SouthBox.p", "rb"))

        for k in dic.keys():
            coll = []
            for ll in dic[k]:
                coll.extend(ll)
            dic[k] = coll
        print(dic.keys())

        ccinput = dic['c30']
        rrinput = dic['r30']

        cinput.extend(ccinput)
        rinput.extend(rrinput)

    cinput = np.array(cinput)
    rinput = np.array(rinput)

    cinput = cinput[np.isfinite(cinput)]
    rinput = rinput[np.isfinite(rinput)]

    nbpoint, pointcount, bins = u_stat.histo_frequency(cinput, bins=np.arange(-15,15,1))
    nball, allcount, bins = u_stat.histo_frequency(rinput, bins=np.arange(-15, 15, 1))
    print(bins)
    bin_centre = bins[0:-1] + ((bins[1::] - bins[0:-1]) / 2)
    bin_edge = bins[0:-1]
    width = bins[1::] - bins[0:-1]



    # nball, bins,v = plt.hist(cinput, bins=np.arange(-7,7,0.5), normed=True, edgecolor='k', color=None, alpha=0.3)
    # nbpoint, bins, v = plt.hist(rinput, bins=np.arange(-7,7,0.5), normed=True, edgecolor='k', color=None, alpha=0.3)
    # plt.xlabel('Local wavelet coefficient, 30km')

    ax.bar(bin_edge, nbpoint, label='core', edgecolor='k', alpha=0.5, align='edge', width=width)
    ax.bar(bin_edge, nball, label='core', edgecolor='k', alpha=0.5, align='edge', width=width)
    plt.ylabel('Frequency')
    stri = (np.sum(cinput <= np.percentile(rinput, 10)) / cinput.size * 100).round(2)
    plt.title(str(hour)+'UTC:'+str(stri)+'% of Cells occur in warmest decile')

    stri = (np.sum(cinput <= np.percentile(rinput, 50)) / cinput.size * 100).round(2)
    print(str(hour)+'UTC:' + str(stri) + '% of Cells occur in warmest half')

    cumulative_random = np.cumsum(nball)
    cumulative = np.cumsum(nbpoint)

    ax = f.add_subplot(122)
    plt.title(str(hour))
    plt.plot(bin_centre,cumulative, label='cores')
    plt.plot(bin_centre, cumulative_random, label='random')
    plt.axvline(0,ymin=0, ymax=1, linestyle='dashed', color='k')
    plt.legend()

def plot_rand_separate():

    f = plt.figure()
    ax = f.add_subplot(121)

    for hour in [17,21,3]:
        path = cnst.network_data + 'figs/LSTA/corrected_LSTA/new/wavelet_coefficients'
        dic = pkl.load(open(path+"/LSTA_histograms_AMSRE_"+str(hour).zfill(2)+"_corrected_SouthBox.p", "rb"))

        for k in dic.keys():
            coll = []
            for ll in dic[k]:
                coll.extend(ll)
            dic[k] = coll
        print(dic.keys())

        cinput = dic['c30']
        rinput = dic['r30']

        cinput = np.array(cinput)
        rinput = np.array(rinput)

        cinput = cinput[np.isfinite(cinput)]
        rinput = rinput[np.isfinite(rinput)]

        nbpoint, pointcount, bins = u_stat.histo_frequency(cinput, bins=np.arange(-15,15,1))
        nball, allcount, bins = u_stat.histo_frequency(rinput, bins=np.arange(-15, 15, 1))
        print(bins)
        bin_centre = bins[0:-1] + ((bins[1::] - bins[0:-1]) / 2)
        bin_edge = bins[0:-1]
        width = bins[1::] - bins[0:-1]



        # nball, bins,v = plt.hist(cinput, bins=np.arange(-7,7,0.5), normed=True, edgecolor='k', color=None, alpha=0.3)
        # nbpoint, bins, v = plt.hist(rinput, bins=np.arange(-7,7,0.5), normed=True, edgecolor='k', color=None, alpha=0.3)
        # plt.xlabel('Local wavelet coefficient, 30km')

        ax.bar(bin_edge, nbpoint, label='core', edgecolor='k', alpha=0.5, align='edge', width=width)
        #ax.bar(bin_edge, nball, label='core', edgecolor='k', alpha=0.4, align='edge', width=width)
        plt.ylabel('Frequency')
        stri = (np.sum(cinput <= np.percentile(rinput, 10)) / cinput.size * 100).round(2)
        plt.title(str(hour)+'UTC:'+str(stri)+'% of Cells occur in warmest decile')

        stri = (np.sum(cinput <= np.percentile(rinput, 50)) / cinput.size * 100).round(2)
        print(str(hour)+'UTC:' + str(stri) + '% of Cells occur in warmest half')

        cumulative_random = np.cumsum(nball)
        cumulative = np.cumsum(nbpoint)

        ax1 = f.add_subplot(122)
        plt.title(str(hour))
        ax1.plot(bin_centre,cumulative, label=str(hour))
        #ax1.plot(bin_centre, cumulative_random, label=str(hour))
    plt.axvline(0,ymin=0, ymax=1, linestyle='dashed', color='k')
    plt.legend()


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
        dic = pkl.load(open(path + "/LSTA_histograms_AMSRE_" + str(h).zfill(2) + "_corrected_SouthBox.p", "rb"))
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


        p = 75
        print('Random sample 75 centile', np.percentile(all, p))
        pprob = np.sum(point > np.percentile(all, p))
        prob = pprob / point.size
        ipdb.set_trace()
        percmax = (prob - (100-p)*0.01) / ((100-p)*0.01) *100 # percentage of cells in warmest 25% of LSTA
        percmmax.append(percmax)
        nbmax.append(point > np.percentile(all, p))

        low90, upp90 = proportion_confint(pprob, point.size)

        err90_up.append( ((upp90 - (100-p)*0.01) / ((100-p)*0.01) *100) - percmax)
        err90_low.append(percmax -((low90 - (100 - p) * 0.01) / ((100 - p) * 0.01) * 100) )

        p = 25
        print('Random sample 25 centile', np.percentile(all, p))
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
    plt.savefig(path + '/paper/AMSRE_core_probability_100kmScale_upstream.png')


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
            path = cnst.network_data + 'figs/LSTA/corrected_LSTA/new/wavelet_coefficients/'
            dic = pkl.load(open(path + "LSTA_histograms_AMSRE_"+str(h).zfill(2)+"_corrected_SouthBox.p", "rb"))
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

            print(h, '10prob', prob)
            print(h, 'percent increase', (prob-0.1)/0.1*100)

            percmmin.append(percmin)
            nbmin.append(len(point))
            low10, upp10 = proportion_confint(pprob, point.size)

            err10_up.append((upp10 - p*0.01) / (p*0.01) *100 - percmin)
            err10_low.append( percmin - (low10 - p*0.01) / (p*0.01) *100 )

        ax = f.add_subplot(3,1,ids+1)
        ax.bar(np.arange(0,len(rrange)), percmmin,  label='25th centile',yerr=np.vstack((err10_up, err10_low)), edgecolor='k', alpha=0.5) #
        ax.bar(np.arange(0, len(rrange)), percmmax, label='75th centile', yerr=np.vstack((err90_up, err90_low)), edgecolor='k')
        ax.set_xticks(np.arange(0, len(rrange)))
        ax.set_xticklabels(rrange)

        ax.set_xlabel('Hour')

        plt.ylabel('Difference in probability (%)')
        plt.legend()

        ax1 = ax.twiny()
        ax1.bar(np.arange(0, len(rrange)), percmmin, label='25th centile', yerr=np.vstack((err10_up, err10_low)), edgecolor='k', alpha=0.5)
        ax1.bar(np.arange(0, len(rrange)), percmmax, label='75th centile', yerr=np.vstack((err90_up, err90_low)), edgecolor='k')
        ax1.set_xticks(np.arange(0,len(rrange)))
        ax1.set_xticklabels(nbmin, rotation=45)
        ax1.set_xlabel('Number of convective cores')

        plt.title(input[2])

    plt.tight_layout()
    #plt.annotate('a)', xy=(0.04, 0.94), xytext=(0, 4), size=15, xycoords=('figure fraction', 'figure fraction'),
    #             textcoords='offset points')  # transform=ax.transAxes,
    plt.show()
    plt.savefig(path + '/paper/AMSRE_core_probability_triple.png')


def plot_diurn_double():


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
            dic = pkl.load(open(path + "LSTA_histograms_AMSRE_"+str(h).zfill(2)+"SlotFilter_+150km_validCheck.p", "rb"))
            dic = pkl.load(
                open(path + "LSTA_histograms_AMSRE_" + str(h).zfill(2) + "_corrected_SouthBox.p", "rb"))
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

            print(h, '10prob', prob)
            print(h, 'percent increase', (prob-0.1)/0.1*100)

            percmmin.append(percmin)
            nbmin.append(len(point))
            low10, upp10 = proportion_confint(pprob, point.size)

            err10_up.append((upp10 - p*0.01) / (p*0.01) *100 - percmin)
            err10_low.append( percmin - (low10 - p*0.01) / (p*0.01) *100 )

        ax = f.add_subplot(2,1,ids+1)
        ax.bar(np.arange(0,len(rrange)), percmmin,  label='25th centile',yerr=np.vstack((err10_up, err10_low)), edgecolor='k') #
        ax.bar(np.arange(0, len(rrange)), percmmax, label='75th centile', yerr=np.vstack((err90_up, err90_low)), edgecolor='k')
        ax.set_xticks(np.arange(0, len(rrange)))
        ax.set_xticklabels(rrange)
        plt.ylim(-75,150)
        lw = 0.5
        plt.axhline(y=-50, linewidth=lw, color='k', linestyle='dashed')
        plt.axhline(y=-25, linewidth=lw, color='k', linestyle='dashed')
        plt.axhline(y=25, linewidth=lw, color='k', linestyle='dashed', zorder=0)
        plt.axhline(y=50, linewidth=lw, color='k', linestyle='dashed', zorder=0)
        plt.axhline(y=75, linewidth=lw, color='k', linestyle='dashed', zorder=0)
        plt.axhline(y=100, linewidth=lw, color='k', linestyle='dashed', zorder=0)
        plt.axhline(y=125, linewidth=lw, color='k', linestyle='dashed', zorder=0)
        plt.axhline(y=0, linewidth=1, color='k', linestyle='solid', zorder=0)

        ax.set_xlabel('Hour')

        plt.ylabel('Difference in probability (%)')
        plt.legend()

        ax1 = ax.twiny()
        ax1.bar(np.arange(0, len(rrange)), percmmin, label='25th centile', yerr=np.vstack((err10_up, err10_low)), edgecolor='k')
        ax1.bar(np.arange(0, len(rrange)), percmmax, label='75th centile', yerr=np.vstack((err90_up, err90_low)), edgecolor='k')
        ax1.set_xticks(np.arange(0,len(rrange)))
        ax1.set_xticklabels(nbmin, rotation=45)
        ax1.set_xlabel('Number of convective cores')

        plt.title(input[2])

    plt.tight_layout()
    #plt.annotate('a)', xy=(0.04, 0.94), xytext=(0, 4), size=15, xycoords=('figure fraction', 'figure fraction'),
    #             textcoords='offset points')  # transform=ax.transAxes,
    plt.show()
    plt.savefig(path + '/initVSprop/AMSRE_core_probability_DOUBLE_100km.png')


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
                open(path + "LSTA_histograms_AMSRE_" + str(h).zfill(2) + "_corrected_SouthBox.p", "rb"))
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
    plt.savefig(path + '/initVSprop/AMSRE_core_probability_RELATIVE_100km.png')


def plot_scatter():

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