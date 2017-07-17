import seaborn as sns
pal = sns.color_palette('Blues')
sns.set_context("paper", font_scale=1.5)
sns.set_style("ticks")
import pdb
import matplotlib.cm as cm

from utils import u_statistics as ug

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl
import statsmodels.stats.proportion as stats
import sys
import scipy.stats as ss


def comp_lat():
    fpath = '/users/global/cornkle/C_paper/wavelet/figs/paper/'
  #  fpath = 'D://data/wavelet/saves/pandas/'
  #  path = 'D://data/wavelet/saves/pandas/'
    path = '/users/global/cornkle/C_paper/wavelet/saves/pandas/'
    dic = pkl.load(open(path+'3dmax_gt15000_TR.p', 'rb'))
    dic2 = pkl.load(open(path+'3dmax_gt15000_noR.p', 'rb'))

    thresh = 10

    hour = np.array(dic['hour'])
    hour2 = np.array(dic2['hour'])

    scales = np.array(dic['scale'])
    psum = np.array(dic['circle_p'])#[(hour>17) & (hour<=23)] #[(hour>15) & (hour<23)]
    psumm = np.array(dic['circle_pc'])
    tmin = np.array(dic['circle_t'])#[(hour>17) & (hour<=23)]
    lat = np.array(dic['clat'])#[(hour>17) & (hour<=23)]
    lon = np.array(dic['clon'])

    pmean = np.array(dic['bulk_pmean'])

    # psum = np.array(dic['circle_p'])[(hour>17) & (hour<=23)] #[(hour>15) & (hour<23)]
    # tmin = np.array(dic['circle_t'])[(hour>17) & (hour<=23)]
    # lat = np.array(dic['clat'])[(hour>17) & (hour<=23)]

    scales2 = np.array(dic2['scale'])
    psum2 = np.array(dic2['circle_p'])[(scales2<=35) ]
    psumm2 = np.array(dic2['circle_pc'])[(scales2 <= 35)]
    tmin2 = np.array(dic2['circle_t'])[(scales2<=35) ]
    lat2 = np.array(dic2['clat'])[(scales2<=35) ]
    lon2 = np.array(dic2['clon'])[(scales2 <= 35)]

    lat22 = np.array(dic2['clat'])

    bins = np.arange(5, 19, 1)  # compute probability per temperature range (1degC)
    print(bins)
    centre = bins[1::]-1.5
    print(centre)

    prob1 = []
    prob2 = []
    low =[]
    up =[]
    low2 =[]
    up2 =[]
    std_error = []


    areas = []

    ids = np.array(dic['id'])
    area = np.array(dic['area'])
    ids2 = np.array(dic2['id'])
    area2 = np.array(dic2['area'])
    ni2, uni2 = np.unique(ids2, return_index = True)
    ni, uni = np.unique(ids, return_index = True)

    print(np.mean(area[uni])*25)
    print(np.mean(area2[uni2])*25)


    for id, b in enumerate(bins):

        if id == 0:
            continue

        p_at_lat = np.concatenate(psum[(lat<b) & (lat>=bins[id-1]) & (lon<=30)])
        t_at_lat = np.concatenate(tmin[(lat < b) & (lat >= bins[id - 1])& (lon<=30)])
        p_t = p_at_lat[t_at_lat<=-78]


        p_at_lat2 = np.concatenate(psum2[(lat2<b) & (lat2>=bins[id-1])& (lon2<=30)])
        t_at_lat2 = np.concatenate(tmin2[(lat2 < b) & (lat2 >= bins[id - 1])& (lon2<=30)])
        p_t2 = p_at_lat2[t_at_lat2<=-78]


        lower, upper = stats.proportion_confint(np.sum(p_t>=thresh), np.sum(np.isfinite(p_t)))
        lower2, upper2 = stats.proportion_confint(np.sum(p_t2 >= thresh), np.sum(np.isfinite(p_t2)))

        print(b, bins[id-1])

        ar = area[(lat<b) & (lat>=bins[id-1])& (lon<=30)]
        ii = ids[(lat<b) & (lat>=bins[id-1])& (lon<=30)]
        iuni, induni = np.unique(ii, return_index = True)
        a_lat = np.mean(np.array(ar[induni]))
        std = ss.sem(np.array(ar[induni]))

        std_error.append(std)
        prob1.append(np.sum(p_t>=thresh) / np.sum(np.isfinite(p_t)))
        prob2.append(np.sum(p_t2 >= thresh) / np.sum(np.isfinite(p_t2)))
        low.append(lower)
        low2.append(lower2)
        up.append(upper)
        up2.append(upper2)
        areas.append(a_lat*25)

    f = plt.figure()
    ax1 = f.add_subplot(111)

    ax1.plot(centre, np.array(prob1)* 100-2,  linewidth=1.5 , marker='o', label='Temperature only')
    ax1.plot(centre, np.array(prob2)* 100+3,  linewidth=1.5 , marker='o', color='r', label='Scales < 35km')
    ax1.legend()
    ax1.set_title('Probability Precip>30mm')
    ax1.fill_between(centre, np.array(low) * 100-2, np.array(up) * 100-2, alpha=0.3)
    ax1.fill_between(centre, np.array(low2) * 100+3, np.array(up2) * 100+3, alpha=0.3, color='r')

    ax11 = ax1.twinx()
    ax11.plot(centre, areas, linestyle='', marker='o', color='grey')
    #ax11.fill_between(centre, areas-std, areas+std, alpha=0.3)


    #plt.tight_layout()
    #plt.savefig(fpath + 'lat_prob.png')
    # plt.savefig(path + 'wavelet_scale_p_T.pdf')
    #plt.close('all')
    return centre, prob1, prob2, low, up, low2, up2, areas
