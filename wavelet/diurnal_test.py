import seaborn as sns
import pickle as pkl
pal = sns.color_palette('Blues')
sns.set_context("paper", font_scale=1.5)
sns.set_style("ticks")
import ipdb
import math
from mpl_toolkits.mplot3d import Axes3D
from utils import u_graphics as ug

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def run():
    df = pd.read_pickle('/users/global/cornkle/C_paper/wavelet/saves/pandas/3dmax_gt15000.pkl')

    scales = df['scale'].unique()

  #   scales = np.array([15,20,30,40,50,60,90,120,150,180,202])
  # #  scales = np.array([19,24,30,40, 50, 67, 85, 113, 151, 202])#[15,20,30,40,50,60,90,120,150,180,202])
  #   #scales = np.arange(15,202,20)
    scenter = scales[0:-1] + (scales[1::] - scales[0:-1]) / 2
  #   print(np.max(df['tmin']), np.min(df['tmin']) )
  #   print(len(df['tmin']))
  #
  #
   ############## PLOT ALL SYSTEMS, AVERAGE MAX SCALE Fig1
    ll1 = []
    ll2 = []
    ll3 = []
    ll4 = []

    for iind, s in enumerate(scales):
        if iind == 0:
            continue
        sums = np.sum(df['sum'][(df['scale']>=scales[iind-1]) & (df['scale']<s) ])
        sums30 = np.sum(df['sum30'][(df['scale']>=scales[iind-1]) & (df['scale']<s) ])
        nz = np.sum(df['sumnz'][(df['scale']>=scales[iind-1]) & (df['scale']<s) ])
        valid = np.sum(df['sumvalid'][(df['scale']>=scales[iind-1]) & (df['scale']<s) ])

        if (nz == 0) or (valid == 0):
            nz = np.nan
            valid = np.nan
        mean = sums / nz
        exmean = sums30 / nz
        meanval = sums / valid
        exmeanval = sums30 / valid

        ll1.append(mean)
        ll2.append(exmean)
        ll3.append(meanval)
        ll4.append(exmeanval)

    f = plt.figure()

    ax = f.add_subplot(221)
    plt.title('Precip / nz: all scales')
    plt.scatter(scenter, ll1, label='Sum/ nozero')
    plt.text(0.7, 0.3, scenter[np.argmax(ll1)], transform=ax.transAxes, fontsize=15)
    plt.legend()

    ax = f.add_subplot(222)
    plt.scatter(scenter, ll2, label='Extreme / nozero')
    plt.legend()
    plt.text(0.7, 0.3, scenter[np.argmax(ll2)], transform=ax.transAxes, fontsize=15)

    ax = f.add_subplot(223)
    plt.scatter(scenter, ll3, label='sum / valid')
    plt.text(0.7, 0.3, scenter[np.argmax(ll3)], transform=ax.transAxes, fontsize=15)
    plt.legend()

    ax = f.add_subplot(224)
    plt.scatter(scenter, ll4, label='Extreme / valid')
    plt.text(0.7, 0.3, scenter[np.argmax(ll4)], transform=ax.transAxes, fontsize=15)
    plt.legend()

  #
  #   ############## PLOT SYSTEM IN AREA BINS, MAX SCALE PER BIN, MEAN RAIN FOR MAX SCALE Fig2
  #
  #   # In[110]:
  #   sizes = np.array(
  #       [15000, 25000, 35000, 50000, 75000, 100000, 150000, 250000, 300000, 400000,
  #        500000])
  #   center = sizes[0:-1] + (sizes[1::] - sizes[0:-1]) / 2
  #
  #   max_ll1 = []
  #   max_ll2 = []
  #   max_ll3 = []
  #   max_ll4 = []
  #   max_ll5 = []
  #   max_ll6 = []
  #   pmaxx1 = []
  #   pmaxx2 = []
  #   pmaxx3 = []
  #   pmaxx4 = []
  #   nbb=[]
  #   for ind, siz in enumerate(sizes):
  #       ll1 = []
  #       ll2 = []
  #       ll3 = []
  #       ll4 = []
  #       ll5 = []
  #       ll6 = []
  #       nnb=[]
  #
  #       if ind == 0:
  #           continue
  #       for iind, s in enumerate(scales):
  #           if iind == 0:
  #               continue
  #           sums = np.nansum(df['sum'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['area'] * 25 >= sizes[ind - 1]) & (df['area'] * 25 < siz) ])
  #           sums30 = np.nansum(
  #               df['sum30'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['area'] * 25 >= sizes[ind - 1]) & (df['area'] * 25 < siz)  ])
  #           nz = np.nansum(df['sumnz'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['area'] * 25 >= sizes[ind - 1]) & (df['area'] * 25 < siz) ])
  #           valid = np.nansum(
  #               df['sumvalid'][(df['scale']>=scales[iind-1]) & (df['scale']<s) &(df['area'] * 25 >= sizes[ind - 1]) & (df['area'] * 25 < siz)  ])
  #
  #           nb = len(df['scale'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['area'] * 25 >= sizes[ind - 1]) & (df['area'] * 25 < siz)  ])
  #
  #
  #           if (nz == 0) or (valid == 0):
  #               nz = np.nan
  #               valid = np.nan
  #
  #           mean = sums / nz
  #           exmean = sums30 / nz
  #           meanval = sums / valid
  #           exmeanval = sums30 / valid
  #           ll1.append(mean)
  #           ll2.append(exmean)
  #           ll3.append(meanval)
  #           ll4.append(exmeanval)
  #           ll5.append(sums)
  #           ll6.append(valid)
  #           nnb.append(nb)
  #
  #
  #       max_ll1.append(scenter[np.argmax(ll1)])
  #       max_ll2.append(scenter[np.argmax(ll2)])
  #       max_ll3.append(scenter[np.argmax(ll3)])
  #       max_ll4.append(scenter[np.argmax(ll4)])
  #       max_ll5.append(ll5[np.argmax(ll3)])
  #       max_ll6.append(ll6[np.argmax(ll3)])
  #       nbb.append(nnb[np.argmax(ll1)])
  #       pmaxx1.append(np.max(ll1))
  #       pmaxx2.append(np.max(ll2))
  #       pmaxx3.append(np.max(ll3))
  #       pmaxx4.append(np.max(ll4))
  #
  #
  #   ##############  MEAN RAIN FOR MAX SCALE Fig3
  #
  #   f = plt.figure()
  #
  #   ax = f.add_subplot(121)
  #   plt.scatter(center, max_ll1)
  #  # plt.scatter(center, np.sqrt(np.array(center)/math.pi)*2, color='r')
  #   plt.title('All points, Mean nozero')
  #
  #   ax = f.add_subplot(122)
  #   plt.scatter(center, max_ll2)
  #  # plt.scatter(center, np.sqrt(np.array(center) / math.pi) * 2, color='r')
  #   plt.title('Extreme nozero')
  #
  #
  #
  #   f = plt.figure()
  #
  #   ax = f.add_subplot(121)
  #   plt.title('Mean precip for max scale')
  #   plt.scatter(center, pmaxx1)
  #
  #   plt.title('Rain Mean nozero')
  #
  #   ax = f.add_subplot(122)
  #   plt.scatter(center, pmaxx2)
  #
  #   plt.title('Probability Extreme nozero')
  #
  #
  #   ############## PLOT SYSTEM IN AREA BINS, MAX SCALE PER BIN <70
  #
  #   # In[110]:
  #   sizes = np.array(
  #       [15000, 25000, 35000, 50000, 75000, 100000, 150000, 250000, 300000, 400000,
  #        500000])
  #   center = sizes[0:-1] + (sizes[1::] - sizes[0:-1]) / 2
  #
  #   max_ll1 = []
  #   max_ll2 = []
  #   max_ll3 = []
  #   max_ll4 = []
  #   max_ll5 = []
  #   max_ll6 = []
  #   pmaxx1 = []
  #   pmaxx2 = []
  #   pmaxx3 = []
  #   pmaxx4 = []
  #   max_sc = []
  #   for ind, siz in enumerate(sizes):
  #       ll1 = []
  #       ll2 = []
  #       ll3 = []
  #       ll4 = []
  #       ll5 = []
  #       ll6 = []
  #
  #       msc = []
  #       if siz == 15000:
  #           continue
  #       for iind, s in enumerate(scales):
  #           if iind == 0:
  #               continue
  #           sums = np.nansum(df['sum'][
  #                                (df['scale'] >= scales[iind - 1]) & (df['scale'] < s) & (df['area'] * 25 >= sizes[ind - 1]) & (df['area'] * 25 < siz) & (
  #                                df['tmin'] <= -68)])
  #           sums30 = np.nansum(
  #               df['sum30'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['area'] * 25 >= sizes[ind - 1]) & (df['area'] * 25 < siz) & (
  #                   df['tmin'] <= -68)])
  #           nz = np.nansum(df['sumnz'][
  #                              (df['scale'] >= scales[iind - 1]) & (df['scale'] < s) & (df['area'] * 25 >= sizes[ind - 1]) & (df['area'] * 25 < siz) & (
  #                                  df['tmin'] <= -68)])
  #           valid = np.nansum(
  #               df['sumvalid'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['area'] * 25 >= sizes[ind - 1]) & (df['area'] * 25 < siz) & (
  #                   df['tmin'] <= -68)])
  #
  #           #   tt =  df['tmin'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['area'] * 25 >= sizes[ind - 1]) & (df['area'] * 25 < siz) & (df['tmin'] <= 70)])
  #
  #           if (nz == 0) or (valid == 0):
  #               nz = np.nan
  #               valid = np.nan
  #
  #           mean = sums / nz
  #           exmean = sums30 / nz
  #           meanval = sums / valid
  #           exmeanval = sums30 / valid
  #           ll1.append(mean)
  #           ll2.append(exmean)
  #           ll3.append(meanval)
  #           ll4.append(exmeanval)
  #           ll5.append(sums30)
  #           ll6.append(nz)
  #
  #       count, nb = np.histogram(msc)
  #       max_scale = nb[np.argmax(count)]
  #
  #       max_ll1.append(scenter[np.argmax(ll1)])
  #       max_ll2.append(scenter[np.argmax(ll2)])
  #       max_ll3.append(scenter[np.argmax(ll3)])
  #       max_ll4.append(scenter[np.argmax(ll4)])
  #       max_ll5.append(ll5[np.argmax(ll2)])
  #       max_ll6.append(ll6[np.argmax(ll2)])
  #       pmaxx1.append(np.max(ll1))
  #       pmaxx2.append(np.max(ll2))
  #       pmaxx3.append(np.max(ll3))
  #       pmaxx4.append(np.max(ll4))
  #
  #
  #   f = plt.figure()
  #
  #   ax = f.add_subplot(221)
  #   plt.scatter(center, max_ll1)
  #   plt.title('Mean nozero | le -65')
  #
  #   ax = f.add_subplot(222)
  #   plt.scatter(center, max_ll2)
  #   plt.title('Extreme nozero | le -65')
  #
  #
  #   ax = f.add_subplot(223)
  #   plt.scatter(center, max_ll5)
  #   plt.title('Prob at max scale(nozero)')
  #
  #   ax = f.add_subplot(224)
  #   plt.scatter(center, max_ll6)
  #   plt.title('Nz pix at max scale')
  #
  #   ############## PLOT SYSTEM IN AREA BINS, MAX SCALE PER BIN, >=50
  #
  #   # In[110]:
  #   sizes = np.array(
  #       [15000, 25000, 35000, 50000, 75000, 100000, 150000, 250000, 300000, 400000,
  #        500000])
  #   center = sizes[0:-1] + (sizes[1::] - sizes[0:-1]) / 2
  #
  #   max_ll1 = []
  #   max_ll2 = []
  #   max_ll3 = []
  #   max_ll4 = []
  #   max_ll5 = []
  #   max_ll6 = []
  #   pmaxx1 = []
  #   pmaxx2 = []
  #   pmaxx3 = []
  #   pmaxx4 = []
  #   max_sc = []
  #   for ind, siz in enumerate(sizes):
  #       ll1 = []
  #       ll2 = []
  #       ll3 = []
  #       ll4 = []
  #       ll5 = []
  #       ll6 = []
  #
  #       msc = []
  #       if siz == 15000:
  #           continue
  #       for iind, s in enumerate(scales):
  #           if iind == 0:
  #               continue
  #           sums = np.nansum(df['sum'][
  #                                (df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['area'] * 25 >= sizes[ind - 1]) & (df['area'] * 25 < siz) & (
  #                                    df['tmin'] >= -60)])
  #           sums30 = np.nansum(
  #               df['sum30'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['area'] * 25 >= sizes[ind - 1]) & (df['area'] * 25 < siz) & (
  #                   df['tmin'] >= -60)])
  #           nz = np.nansum(df['sumnz'][
  #                              (df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['area'] * 25 >= sizes[ind - 1]) & (df['area'] * 25 < siz) & (
  #                                  df['tmin'] >= -60)])
  #           valid = np.nansum(
  #               df['sumvalid'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['area'] * 25 >= sizes[ind - 1]) & (df['area'] * 25 < siz) & (
  #                   df['tmin'] >= -60)])
  #
  #           #   tt =  df['tmin'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['area'] * 25 >= sizes[ind - 1]) & (df['area'] * 25 < siz) & (df['tmin'] <= 70)])
  #
  #           if (nz == 0) or (valid == 0):
  #               nz = np.nan
  #               valid = np.nan
  #
  #           mean = sums / nz
  #           exmean = sums30 / nz
  #           meanval = sums / valid
  #           exmeanval = sums30 / valid
  #           ll1.append(mean)
  #           ll2.append(exmean)
  #           ll3.append(meanval)
  #           ll4.append(exmeanval)
  #           ll5.append(sums30)
  #           ll6.append(nz)
  #
  #       count, nb = np.histogram(msc)
  #       max_scale = nb[np.argmax(count)]
  #
  #       max_ll1.append(scenter[np.argmax(ll1)])
  #       max_ll2.append(scenter[np.argmax(ll2)])
  #       max_ll3.append(scenter[np.argmax(ll3)])
  #       max_ll4.append(scenter[np.argmax(ll4)])
  #       max_ll5.append(ll5[np.argmax(ll2)])
  #       max_ll6.append(ll6[np.argmax(ll2)])
  #       pmaxx1.append(np.max(ll1))
  #       pmaxx2.append(np.max(ll2))
  #       pmaxx3.append(np.max(ll3))
  #       pmaxx4.append(np.max(ll4))
  #
  #   f = plt.figure()
  #
  #   ax = f.add_subplot(221)
  #   plt.scatter(center, max_ll1)
  #   plt.title('Mean nozero | ge -60')
  #
  #   ax = f.add_subplot(222)
  #   plt.scatter(center, max_ll2)
  #   plt.title('Extreme nozero | ge -60')
  #
  #   ax = f.add_subplot(223)
  #   plt.scatter(center, max_ll5)
  #   plt.title('Prob at max scale(nz) | ge -60')
  #
  #   ax = f.add_subplot(224)
  #   plt.scatter(center, max_ll6)
  #   plt.title('Nz pix at max scale | ge -60')
  #
  #   # In[88]:
  #   ############## PLOT SYSTEM IN LARGER AREA BINS Fig 4
  #   ll1 = []
  #   for iind, s in enumerate(scales):
  #       if iind == 0:
  #           continue
  #       sums = np.sum(df['sum30'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['area'] * 25 >= 50000) & (df['area'] * 25 <= 250000) & (df['tmin'] <= -65)])
  #       nz = np.sum(df['sumnz'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['area'] * 25 >= 50000) & (df['area'] * 25 <= 250000) & (df['tmin'] <= -65)])
  #       if (nz == 0):
  #           nz = np.nan
  #       mean = sums / nz
  #       ll1.append(mean)
  #
  #   ll2 = []
  #   for iind, s in enumerate(scales):
  #       if iind == 0:
  #           continue
  #       sums = np.sum(df['sum30'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['area'] * 25 < 50000) & (df['area'] * 25 >= 25000)& (df['tmin'] <= -65)])
  #       nz = np.sum(df['sumnz'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['area'] * 25 < 50000) & (df['area'] * 25 >= 25000)& (df['tmin'] <= -65)])
  #       if (nz == 0):
  #           nz = np.nan
  #       mean = sums / nz
  #       ll2.append(mean)
  #
  #   ll3 = []
  #   for iind, s in enumerate(scales):
  #       if iind == 0:
  #           continue
  #       sums = np.sum(df['sum30'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['area'] * 25 >= 250000)& (df['tmin'] <= -65)])
  #       nz = np.sum(df['sumnz'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['area'] * 25 >= 250000)& (df['tmin'] <= -65)])
  #       if (nz == 0):
  #           nz = np.nan
  #       mean = sums / nz
  #       ll3.append(mean)
  #
  #   ll4 = []
  #   for iind, s in enumerate(scales):
  #       if iind == 0:
  #           continue
  #       sums = np.sum(df['sum30'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['area'] * 25 < 25000) & (df['area'] * 25 >= 15000)& (df['tmin'] <= -65)])
  #       nz = np.sum(df['sumnz'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['area'] * 25 < 25000) & (df['area'] * 25 >= 15000)& (df['tmin'] <= -65)])
  #       if (nz == 0):
  #           nz = np.nan
  #       mean = sums / nz
  #       ll4.append(mean)
  #
  #
  #   f = plt.figure()
  #
  #   ax = f.add_subplot(221)
  #   plt.scatter(scenter, ll4, label='ge 15000 lt 25000 | lt-65')
  #   plt.text(0.8, 0.9, str(scales[np.argmax(ll4)]), transform=ax.transAxes, fontsize=10)
  #   plt.legend()
  #   plt.title('Precip sum / nozero points')
  #   ax = f.add_subplot(222)
  #   plt.scatter(scenter, ll2, label='ge 25000 lt 50000 | lt-65')
  #   plt.text(0.8, 0.9, str(scales[np.argmax(ll2)]), transform=ax.transAxes, fontsize=10)
  #   plt.legend()
  #   plt.title('Precip sum / nozero points')
  #   ax = f.add_subplot(223)
  #   plt.scatter(scenter, ll1, label='ge 50000 lt 250000 | lt-65')
  #   plt.text(0.8, 0.9, str(scales[np.argmax(ll1)]), transform=ax.transAxes, fontsize=10)
  #   plt.legend()
  #   ax = f.add_subplot(224)
  #   plt.scatter(scenter, ll3, label='gt 250000')
  #   plt.text(0.8, 0.9, str(scales[np.argmax(ll3)]), transform=ax.transAxes, fontsize=10)
  #   plt.legend()
  #
  #   # In[88]:
    ############## PLOT SYSTEM IN LARGER AREA BINS Fig 4
    ll1 = []
    for iind, s in enumerate(scales):
        if iind == 0:
            continue
        sums = np.sum(df['sum30'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['area'] * 25 >= 50000) & (df['area'] * 25 <= 75000) ])
        nz = np.sum(df['sumnz'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['area'] * 25 >= 50000) & (df['area'] * 25 <= 75000) ])
        if (nz == 0):
            nz = np.nan
        mean = sums / nz
        ll1.append(mean)

    ll2 = []
    for iind, s in enumerate(scales):
        if iind == 0:
            continue
        sums = np.sum(df['sum30'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['area'] * 25 < 50000) & (df['area'] * 25 >= 25000)])
        nz = np.sum(df['sumnz'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['area'] * 25 < 50000) & (df['area'] * 25 >= 25000)])
        if (nz == 0):
            nz = np.nan
        mean = sums / nz
        ll2.append(mean)

    ll3 = []
    for iind, s in enumerate(scales):
        if iind == 0:
            continue
        sums = np.sum(df['sum30'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['area'] * 25 >= 250000)])
        nz = np.sum(df['sumnz'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['area'] * 25 >= 250000)])
        if (nz == 0):
            nz = np.nan
        mean = sums / nz
        ll3.append(mean)

    ll4 = []
    for iind, s in enumerate(scales):
        if iind == 0:
            continue
        sums = np.sum(df['sum30'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['area'] * 25 < 25000) & (df['area'] * 25 >= 15000)])
        nz = np.sum(df['sumnz'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['area'] * 25 < 25000) & (df['area'] * 25 >= 15000)])
        if (nz == 0):
            nz = np.nan
        mean = sums / nz
        ll4.append(mean)

    ll5 = []
    for iind, s in enumerate(scales):
        if iind == 0:
            continue
        sums = np.sum(df['sum'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['area'] * 25 < 100000) & (df['area'] * 25 >= 75000)])
        nz = np.sum(df['sumnz'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['area'] * 25 < 100000) & (df['area'] * 25 >= 75000)])
        if (nz == 0):
            nz = np.nan
        mean = sums / nz
        ll5.append(mean)

    ll6 = []
    for iind, s in enumerate(scales):
        if iind == 0:
            continue
        sums = np.sum(df['sum'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['area'] * 25 < 250000) & (df['area'] * 25 >= 100000)])
        nz = np.sum(df['sumnz'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['area'] * 25 < 250000) & (df['area'] * 25 >= 100000)])
        if (nz == 0):
            nz = np.nan
        mean = sums / nz
        ll6.append(mean)


    f = plt.figure()

    ax = f.add_subplot(231)
    plt.scatter(scenter, ll4, label='ge 15000 lt 25000 | ge -60')
    plt.text(0.8, 0.9, str(scenter[np.argmax(ll4)]), transform=ax.transAxes, fontsize=10)
    plt.legend()
    plt.title('Precip sum / nozero points')
    ax = f.add_subplot(232)
    plt.scatter(scenter, ll2, label='ge 25000 lt 50000')
    plt.text(0.8, 0.9, str(scenter[np.argmax(ll2)]), transform=ax.transAxes, fontsize=10)
    plt.legend()
    plt.title('Precip sum / nozero points')
    ax = f.add_subplot(233)
    plt.scatter(scenter, ll1, label='ge 50000 lt 75000')
    plt.text(0.8, 0.9, str(scenter[np.argmax(ll1)]), transform=ax.transAxes, fontsize=10)
    plt.legend()
    ax = f.add_subplot(236)
    plt.scatter(scenter, ll3, label='gt 250000')
    plt.text(0.8, 0.9, str(scenter[np.argmax(ll3)]), transform=ax.transAxes, fontsize=10)
    plt.legend()
    ax = f.add_subplot(234)
    plt.scatter(scenter, ll5, label='gt 75000 lt 100000')
    plt.text(0.8, 0.9, str(scenter[np.argmax(ll5)]), transform=ax.transAxes, fontsize=10)
    plt.legend()
    ax = f.add_subplot(235)
    plt.scatter(scenter, ll6, label='gt 100000 lt 250000')
    plt.text(0.8, 0.9, str(scenter[np.argmax(ll6)]), transform=ax.transAxes, fontsize=10)
    plt.legend()
  #
  #   # In[21]:
  #
    ############## PLOT SYSTEMS FOR DIURNAL SLICES, NO AREA DIFFERENTIATION Fig 5
    ll1 = []
    ll11=[]
    for iind, s in enumerate(scales):
        if iind == 0:
            continue
        sums = np.sum(df['sum'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & ( (df['hour'] < 6)) & (df['clat'] >=12) ])
        sums30 = np.sum(df['sum30'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & ((df['hour'] < 6)) & (df['clat'] >= 12)])
        nz = np.sum(df['sumnz'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & ( (df['hour'] < 6)) & (df['clat'] >=12) ])
        if (nz == 0):
            nz = np.nan
        mean = sums / nz
        mean30 = sums30 / nz
        ll1.append(mean)
        ll11.append(mean30)

    ll = []
    lll=[]
    for iind, s in enumerate(scales):
        if iind == 0:
            continue
        sums = np.sum(df['sum'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['hour'] >= 15) & (df['hour'] < 21) & (df['clat'] >=12) ])
        sums30 = np.sum(df['sum30'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['hour'] >= 15) & (df['hour'] < 21) & (df['clat'] >= 12)])
        nz = np.sum(df['sumnz'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['hour'] >= 15) & (df['hour'] < 21) & (df['clat'] >=12) ])
        if (nz == 0):
            nz = np.nan
        mean = sums / nz
        mean30 = sums30 / nz
        ll.append(mean)
        lll.append(mean30)

    f = plt.figure()
    ax = f.add_subplot(221)
    plt.scatter(scenter, lll)
    plt.title('Prob sum30 / nz 15-21h, Sahel')
    ax = f.add_subplot(222)
    plt.scatter(scenter, ll11)
    plt.title('Prob sum30 / nz 0-6h, Sahel')
    ax = f.add_subplot(223)
    plt.scatter(scenter, ll)
    plt.title('Prob sum / nz 15-21h, Sahel')
    ax = f.add_subplot(224)
    plt.scatter(scenter, ll1)
    plt.title('Prob sum / nz 0-6h, Sahel')
  #
  #   ############## PLOT HISTOGRAMS FOR SYSTEM NUMBER PER DIURNAL SLICES, OVERALL SYSTEM NUMBER Fig 6
  #
  #   bla = df['area'][((df['hour'] >= 21) ^ (df['hour'] < 3)) & (df['clat'] > 12) ] * 25
  #   bla1 = df['area'][(df['hour'] >= 15) & (df['hour'] < 21) & (df['clat'] > 12)] * 25
  #   weight = np.ones_like(bla.unique()) / float(len(bla.unique()))
  #   weight1 = np.ones_like(bla1.unique()) / float(len(bla1.unique()))
  #
  #   bla2 = df['area'][((df['hour'] >= 21) ^ (df['hour'] < 3)) & (df['clat'] > 8) & (df['clat'] < 12) ] * 25
  #   bla3 = df['area'][(df['hour'] >= 15) & (df['hour'] < 21) & (df['clat'] > 8) & (df['clat'] < 12)] * 25
  #   weight2 = np.ones_like(bla2.unique()) / float(len(bla2.unique()))
  #   weight3 = np.ones_like(bla3.unique()) / float(len(bla3.unique()))
  #
  #   bla4 = df['area'][((df['hour'] >= 21) ^ (df['hour'] < 3)) & (df['clat'] <=8 ) & (df['clat'] >=4 ) ] * 25
  #   bla5 = df['area'][(df['hour'] >= 15) & (df['hour'] < 21) & (df['clat'] <=8) & (df['clat'] >=4 ) ] * 25
  #   weight4 = np.ones_like(bla4.unique()) / float(len(bla4.unique()))
  #   weight5 = np.ones_like(bla5.unique()) / float(len(bla5.unique()))
  #
  #
  #   f = plt.figure()
  #   ax = f.add_subplot(221)
  #   plt.hist(bla.unique(), bins=30, range=(15000, 500000), weights=weight, label='21-3h')
  #   plt.hist(bla1.unique(), bins=30, color='r', alpha=0.5, range=(15000, 500000), weights=weight1, label='15-21h')
  #   plt.title('Sahel: gt12N')
  #   plt.legend()
  #
  #   ax = f.add_subplot(222)
  #   plt.hist(bla2.unique(), bins=30, range=(15000, 500000), weights=weight2, label='21-3h')
  #   plt.hist(bla3.unique(), bins=30, color='r', alpha=0.5, range=(15000, 500000), weights=weight3, label='15-21h')
  #   plt.title('Sudano-Sahel 8-12N')
  #
  #   ax = f.add_subplot(223)
  #   plt.hist(bla4.unique(), bins=30, range=(15000, 500000), weights=weight4, label='21-3h')
  #   plt.hist(bla5.unique(), bins=30, color='r', alpha=0.5, range=(15000, 500000), weights=weight5, label='15-21h')
  #   plt.title('Coast 4-8N')
  #
  #   ax = f.add_subplot(224)
  #   plt.hist(df['area'].unique() * 25, bins=30, range=(15000, 500000), cumulative=True)
  #
  #   f = plt.figure()
  #   ax = f.add_subplot(221)
  #   plt.hist(bla.unique(), bins=30, range=(15000, 500000),  label='21-3h')
  #   plt.hist(bla1.unique(), bins=30, color='r', alpha=0.5, range=(15000, 500000),  label='15-21h')
  #   plt.title('Sahel: gt12N')
  #   plt.legend()
  #
  #   ax = f.add_subplot(222)
  #   plt.hist(bla2.unique(), bins=30, range=(15000, 500000),  label='21-3h')
  #   plt.hist(bla3.unique(), bins=30, color='r', alpha=0.5, range=(15000, 500000),  label='15-21h')
  #   plt.title('Sudano-Sahel 8-12N')
  #
  #   ax = f.add_subplot(223)
  #   plt.hist(bla4.unique(), bins=30, range=(15000, 500000),  label='21-3h')
  #   plt.hist(bla5.unique(), bins=30, color='r', alpha=0.5, range=(15000, 500000),  label='15-21h')
  #   plt.title('Coast 8-12N')
  #
  #   ax = f.add_subplot(224)
  #   plt.hist(df['area'].unique() * 25, bins=30, range=(15000, 500000), cumulative=True)
  #
  #
  #
  #
    ############## PLOT AREA BINNED DIURNAL SLICES, Fig 7
    sizes = np.array(
        [15000, 25000, 35000, 50000, 75000, 100000, 150000, 250000, 300000, 400000,
         500000])
    center = sizes[0:-1] + (sizes[1::] - sizes[0:-1]) / 2

    max_ll1 = []
    max_ll2 = []
    max_ll3 = []
    max_ll4 = []
    pmaxx1 = []
    pmaxx2 = []
    pmaxx3 = []
    pmaxx4 = []
    for ind, siz in enumerate(sizes):
        ll1 = []
        ll2 = []
        ll3 = []
        ll4 = []
        if siz == 15000:
            continue
        for iind, s in enumerate(scales):
            if iind == 0:
                continue
            sums = np.sum(df['sum'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['area'] * 25 >= sizes[ind - 1]) & (df['area'] * 25 < siz) & (df['hour'] >= 15) & (df['hour'] < 21) & (df['clat'] > 10) ])
            sums30 = np.sum(
                df['sum30'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['area'] * 25 >= sizes[ind - 1]) & (df['area'] * 25 < siz) & (df['hour'] >= 15) & (df['hour'] < 21) & (df['clat'] > 10) ])
            nz = np.sum(df['sumnz'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['area'] * 25 >= sizes[ind - 1]) & (df['area'] * 25 < siz) & (df['hour'] >= 15) & (df['hour'] < 21) & (df['clat'] > 10) ])
            valid = np.sum(
                df['sumvalid'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['area'] * 25 >= sizes[ind - 1]) & (df['area'] * 25 < siz) & (df['hour'] >= 15) & (df['hour'] < 21) & (df['clat'] > 10) ])

            if (nz == 0) or (valid == 0):
                nz = np.nan
                valid = np.nan

            mean = sums / nz
            exmean = sums30 / nz
            meanval = sums / valid
            exmeanval = sums30 / valid
            ll1.append(mean)
            ll2.append(exmean)
            ll3.append(meanval)
            ll4.append(exmeanval)

        max_ll1.append(scenter[np.argmax(ll1)])
        max_ll2.append(scenter[np.argmax(ll2)])
        max_ll3.append(scenter[np.argmax(ll3)])
        max_ll4.append(scenter[np.argmax(ll4)])
        pmaxx1.append(np.max(ll1))
        pmaxx2.append(np.max(ll2))
        pmaxx3.append(np.max(ll3))
        pmaxx4.append(np.max(ll4))


    sizes = np.array(
        [15000, 25000, 35000, 50000, 75000, 100000, 150000, 250000, 300000, 400000,
         500000])
    center = sizes[0:-1] + (sizes[1::] - sizes[0:-1]) / 2

    max_ll11 = []
    max_ll22 = []
    max_ll33 = []
    max_ll44 = []
    max_ll55 = []
    max_ll66 = []
    pmaxx11 = []
    pmaxx22 = []
    pmaxx33 = []
    pmaxx44 = []
    for ind, siz in enumerate(sizes):
        ll1 = []
        ll2 = []
        ll3 = []
        ll4 = []
        ll5 = []
        ll6 = []
        if siz == 15000:
            continue
        for iind, s in enumerate(scales):
            if iind == 0:
                continue
            sums = np.sum(df['sum'][
                              (df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['area'] * 25 >= sizes[ind - 1]) & (df['area'] * 25 < siz) &  ((df['hour'] < 6) ) & (df['clat'] > 10) ])
            sums30 = np.sum(
                df['sum30'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['area'] * 25 >= sizes[ind - 1]) & (df['area'] * 25 < siz) & ( (df['hour'] < 6) ) & (df['clat'] > 10) ])
            nz = np.sum(df['sumnz'][
                            (df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['area'] * 25 >= sizes[ind - 1]) & (df['area'] * 25 < siz) & ( (df['hour'] < 6)) & (df['clat'] > 10)])
            valid = np.sum(
                df['sumvalid'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['area'] * 25 >= sizes[ind - 1]) & (df['area'] * 25 < siz) & ( (df['hour'] < 6))& (df['clat'] > 10) ])

            if (nz == 0) or (valid == 0):
                nz = np.nan
                valid = np.nan

            mean = sums / nz
            exmean = sums30 / nz
            meanval = sums / valid
            exmeanval = sums30 / valid
            ll1.append(mean)
            ll2.append(exmean)
            ll3.append(meanval)
            ll4.append(exmeanval)


        max_ll11.append(scenter[np.argmax(ll1)])
        max_ll22.append(scenter[np.argmax(ll2)])
        max_ll33.append(scenter[np.argmax(ll3)])
        max_ll44.append(scenter[np.argmax(ll4)])

        pmaxx11.append(np.max(ll1))
        pmaxx22.append(np.max(ll2))
        pmaxx33.append(np.max(ll3))
        pmaxx44.append(np.max(ll4))

    f = plt.figure()

    ax = f.add_subplot(121)
    plt.title('Precip sum / nb: all scales')
    plt.scatter(center, max_ll1, label='15-21h Sahel')
    plt.scatter(center, max_ll11, color='r', label='21-3h Sahel')
    plt.legend()

    plt.title('Mean nozero')
    plt.legend()

    ax = f.add_subplot(122)
    plt.scatter(center, max_ll2)
    plt.scatter(center, max_ll22, color='r')

    plt.title('Extreme nozero')
    plt.legend()

  #
  #
  #
  #
  #   ############# diurnal cycle Fig 9
  #
  #   # In[110]:
  #   sizes = np.array([0, 3, 6, 9, 12, 15, 18, 21, 24])
  #   center = np.arange(0, 24, 3) + 1.5
  #
  #   max_ll1 = []
  #   max_ll2 = []
  #   max_ll3 = []
  #   max_ll4 = []
  #
  #
  #   for ind, siz in enumerate(sizes):
  #       print(siz)
  #       ll1 = []
  #       ll2 = []
  #       ll3 = []
  #       ll4 = []
  #
  #
  #       if ind == 0:
  #           continue
  #       for iind, s in enumerate(scales):
  #           if iind == 0:
  #               continue
  #           sums = np.nansum(
  #               df['sum'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['hour'] >= sizes[ind - 1]) & (df['hour'] < siz) & (df['clat'] > 8)  ])
  #           sums30 = np.nansum(
  #               df['sum30'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['hour'] >= sizes[ind - 1]) & (df['hour'] < siz) & (df['clat'] > 8)])
  #           nz = np.nansum(
  #               df['sumnz'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['hour'] >= sizes[ind - 1]) & (df['hour'] < siz) & (df['clat'] > 8)])
  #           valid = np.nansum(
  #               df['sumvalid'][(df['scale']>=scales[iind-1]) & (df['scale']<s) & (df['hour'] >= sizes[ind - 1]) & (df['hour'] < siz) & (df['clat'] > 8) ])
  #
  #           if (nz == 0) or (valid == 0):
  #               nz = np.nan
  #               valid = np.nan
  #
  #           mean = sums / nz
  #           exmean = sums30 / nz
  #           meanval = sums / valid
  #           exmeanval = sums30 / valid
  #           ll1.append(mean)
  #           ll2.append(exmean)
  #           ll3.append(meanval)
  #           ll4.append(exmeanval)
  #
  #
  #       max_ll1.append(scenter[np.argmax(ll1)])
  #       max_ll2.append(scenter[np.argmax(ll2)])
  #       max_ll3.append(scenter[np.argmax(ll3)])
  #       max_ll4.append(scenter[np.argmax(ll4)])
  #
  #   f = plt.figure()
  #
  #   ax = f.add_subplot(121)
  #   plt.scatter(center, max_ll1)
  #   plt.title('Mean nozero | le -80 Sahel')
  #
  #   ax = f.add_subplot(122)
  #   plt.scatter(center, max_ll2)
  #   plt.title('Extreme nozero | le -80 Sahel')


  #
  #   snb = []
  #   smean=[]
  #   for ind, siz in enumerate(sizes):
  #
  #       if ind == 0:
  #           continue
  #
  #       area = np.unique(df['area'][(df['hour'] >= sizes[ind - 1]) & (df['hour'] < siz)] * 25)
  #       area = np.unique(df['area'][(df['hour'] >= sizes[ind - 1]) & (df['hour'] < siz)] * 25)
  #       mean = df['scale'][(df['hour'] >= sizes[ind - 1]) & (df['hour'] < siz) & (df['clat'] > 13) & (df['tmin'] < -60)].size
  #       snb.append(area.size)
  #       smean.append(mean)
  #   snb = np.array(snb)
  #   smean = np.array(smean)
  #   smean_perc = smean/np.sum(smean)
  #
  #   f = plt.figure()
  #   plt.scatter(center, snb)
  #
  #   f = plt.figure()
  #
  #   ax = f.add_subplot(131)
  #   colors=['r', 'b', 'g', 'y', 'black']
  #   lab = ['[15-20[', '[20-40[','[40-60[', '[60-90[', '[90-200[']
  #   for l, c in zip(ll1,colors):
  #       plt.plot(center,l, color=c)
  #   plt.plot(center, smean)
  #
  #   ax = f.add_subplot(132)
  #   for l, c, la in zip(ll1, colors, lab):
  #       plt.plot(center, l / np.sum(l), color=c, label=la)
  #   plt.plot(center,smean_perc)
  #
  #
  #   ax = f.add_subplot(133)
  #   for l, c, la in zip(ll1,colors, lab):
  #       plt.plot(center,l/np.sum(l)-smean_perc, color=c, label=la)
  #   plt.legend()
  #   #plt.plot(center, smean_perc)
  #   return
  #
    ############# diurnal cycle size Fig 9

    # In[110]:
    #sizes = np.array([0, 1,2,3,4,5,6,7,8,9,10,11,12,13,14, 15,16,17, 18,19,20, 21,22,23])
    # sizes = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23])
    # center = (np.arange(23) + 1)-0.5
    #
    # ll1 = []
    # std1 = []
    # ll_sahel = []
    # ll_soud = []
    # ll_coast = []
    # for ind, siz in enumerate(sizes):
    #
    #     if ind == 0:
    #         continue
    #
    #     area = np.mean(np.unique(df['area'][(df['hour'] >= sizes[ind - 1]) & (df['hour'] < siz)]) * 25)
    #
    #     area_sa = np.mean(np.unique(df['area'][(df['hour'] >= sizes[ind - 1])  & (df['hour'] < siz) & (
    #     df['clat'] > 12) & (df['clon'] < 10)])*25)
    #     area_su = np.mean(np.unique(df['area'][(df['hour'] >= sizes[ind - 1])  & (df['hour'] < siz) & (
    #     df['clat'] >= 8) & (df['clat'] < 12) & (df['clon'] < 10)])*25)
    #     area_co = np.mean(np.unique(df['area'][(df['hour'] >= sizes[ind - 1])  & (df['hour'] < siz) & (
    #     df['clat'] < 8) & (df['clon'] < 10)])*25)
    #
    #     ll1.append(area)
    #
    #     ll_sahel.append(area_sa)
    #     ll_soud.append(area_su)
    #     ll_coast.append(area_co)
    #
    # f = plt.figure(figsize=(15, 7))
    # ax = f.add_subplot(121)
    # plt.title('Area diurn cycle')
    # plt.scatter(center, ll1, s=30)
    #
    # ax = f.add_subplot(122)
    # plt.title('Regions diurn cycle')
    # plt.scatter(center, ll_sahel, label='Sahel', s=30)
    # plt.scatter(center, ll_soud, label='Sudano-Sahel', color='r', s=30)
    # plt.scatter(center, ll_coast, label='Coast', color='g', s=30)
    # plt.legend()

  #   ############# diurnal cycle size Fig 9
  #
  #   # In[110]:
  #   sizes = np.array([0, 1,2,3,4,5,6,7,8,9,10,11,12,13,14, 15,16,17, 18,19,20, 21,22,23])
  #   center = (np.arange(23) + 1)-0.5
  #
  #   ll1 = []
  #
  #   for ind, siz in enumerate(sizes):
  #
  #       if ind == 0:
  #           continue
  #
  #       area = np.unique(df['area'][ (df['hour'] >= sizes[ind - 1]) & (df['hour'] < siz) ]*25)
  #
  #       ll1.append(area.size)
  #
  #
  #   f = plt.figure()
  #   ax = f.add_subplot(111)
  #   plt.title('System number diurn cycle')
  #   plt.scatter(center, ll1)
  #   print(np.sum(ll1))
  #
  #   ############# diurnal cycle size Fig 9
  #
  #   # In[110]:
  #   sizes = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23])
  #   center = (np.arange(23) + 1)-0.5
  #
  #   ll1 = []
  #   std1=[]
  #   ll_sahel = []
  #   ll_soud = []
  #   ll_coast = []
  #   for ind, siz in enumerate(sizes):
  #
  #       if ind == 0:
  #           continue
  #
  #       area = np.mean(df['tmin'][(df['hour'] >= sizes[ind - 1]) & (df['hour'] < siz) ])
  #       std = np.std(df['tmin'][(df['hour'] >= sizes[ind - 1]) & (df['hour'] < siz) ])
  #
  #       area_sa = np.mean(df['tmin'][(df['hour'] >= sizes[ind - 1])  & (df['hour'] < siz) & (df['clat'] > 12) &  (df['clon'] < 10)  ])
  #       area_su = np.mean(df['tmin'][(df['hour'] >= sizes[ind - 1])  & (df['hour'] < siz) & (df['clat'] >= 8) & (df['clat'] < 12) &  (df['clon'] < 10)  ])
  #       area_co = np.mean(df['tmin'][(df['hour'] >= sizes[ind - 1])  & (df['hour'] < siz) & (df['clat'] < 8) &  (df['clon'] < 10)  ])
  #
  #
  #       ll1.append(area)
  #       std1.append(std)
  #       ll_sahel.append(area_sa)
  #       ll_soud.append(area_su)
  #       ll_coast.append(area_co)
  #
  #   f = plt.figure()
  #
  #   ax = f.add_subplot(121)
  #   plt.title('Tmin diurn cycle')
  #   plt.scatter(center, ll1 , s=30)
  #
  #   ax = f.add_subplot(122)
  #   plt.title('Regions diurn cycle')
  #   plt.scatter(center, ll_sahel, label='Sahel', s=30)
  #   plt.scatter(center, ll_soud, label='Sudano-Sahel', color='r', s=30)
  #   plt.scatter(center, ll_coast, label='Coast', color='g', s=30)
  #   plt.legend()
  #  # plt.errorbar(center, ll1, yerr=np.array(std1)*3)
  #
  #   ############# diurnal cycle size Fig 9
  #
  #   # In[110]:
  #   sizes = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23])
  #   center = (np.arange(23) + 1) - 0.5
  #
  #   ll1 = []
  #   std1 = []
  #   ll_sahel = []
  #   ll_soud = []
  #   ll_coast = []
  #   for ind, siz in enumerate(sizes):
  #
  #       if ind == 0:
  #           continue
  #
  #       area = np.mean(df['tmean'][(df['hour'] >= sizes[ind - 1]) & (df['hour'] < siz)])
  #       std = np.std(df['tmean'][(df['hour'] >= sizes[ind - 1]) & (df['hour'] < siz)])
  #
  #       area_sa = np.mean(df['tmean'][(df['hour'] >= sizes[ind - 1])  & (df['hour'] < siz) & (
  #       df['clat'] > 12) & (df['clon'] < 10)])
  #       area_su = np.mean(df['tmean'][(df['hour'] >= sizes[ind - 1]) & (df['hour'] < siz) & (
  #       df['clat'] >= 8) & (df['clat'] < 12) & (df['clon'] < 10)])
  #       area_co = np.mean(df['tmean'][(df['hour'] >= sizes[ind - 1])  & (df['hour'] < siz) & (
  #       df['clat'] < 8) & (df['clon'] < 10)])
  #
  #       ll1.append(area)
  #       std1.append(std)
  #       ll_sahel.append(area_sa)
  #       ll_soud.append(area_su)
  #       ll_coast.append(area_co)
  #
  #   f = plt.figure()
  #
  #   ax = f.add_subplot(121)
  #   plt.title('Tmean diurn cycle')
  #   plt.scatter(center, ll1, s=30)
  #
  #   ax = f.add_subplot(122)
  #   plt.title('Regions diurn cycle')
  #   plt.scatter(center, ll_sahel, label='Sahel', s=30)
  #   plt.scatter(center, ll_soud, label='Sudano-Sahel', color='r', s=30)
  #   plt.scatter(center, ll_coast, label='Coast', color='g', s=30)
  #   plt.legend()
  #   #  plt.errorbar(center, ll1, yerr=np.array(std1)*3)
  #
  #   ############# diurnal cycle size Fig 9
  #
  #   # In[110]:
  #   sizes = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23])
  #   center = (np.arange(23) + 1)-0.5
  #
  #   ll1 = []
  #
  #   for ind, siz in enumerate(sizes):
  #
  #       if ind == 0:
  #           continue
  #
  #       area = np.mean(df['shape'][(df['hour'] >= sizes[ind - 1]) & (df['hour'] < siz) ])
  #
  #       ll1.append(area)
  #
  #   f = plt.figure()
  #
  #   ax = f.add_subplot(111)
  #   plt.title('Nb composite diurn cycle')
  #   plt.scatter(center, ll1)
  #
  #   ############# diurnal cycle size Fig 9
  #
    # In[110]:
    sizes = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23])
    center = (np.arange(23) + 1) - 0.5

    ll1 = []
    ll2 = []
    ll3 = []

    for ind, siz in enumerate(sizes):

        if ind == 0:
            continue

        sum = np.nansum(df['sum'][(df['hour'] == siz) & (df['scale']<200) & (df['tmin']<-40)])
        sumnz = np.nansum(df['shape'][(df['hour'] == siz)& (df['scale']<200) & (df['tmin']<-40)])

        p = sumnz
        p1 = sum

        ll1.append(p1)
        ll2.append(p)
        ll3.append(p1/p)
        #
    f = plt.figure()
    #
    ax = f.add_subplot(131)
    #   plt.title('Nz diurn cycle')
    plt.plot(center, ll1)
    ax = f.add_subplot(132)
    #   plt.title('Nz diurn cycle')
    plt.plot(center, ll2)
    ax = f.add_subplot(133)
    #   plt.title('Nz diurn cycle')
    plt.plot(center, ll3)
    plt.show()

    ll = []

    sizes = np.unique(df['scale'])
    for ind, siz in enumerate(sizes):


        sum = df['scale'][df['scale']==siz].size

        ll.append(sum)

    f = plt.figure()
    #
    ax = f.add_subplot(111)
    plt.plot(sizes, ll)

    #   plt.title('Probability diurn cycle')
    #   plt.scatter(center, ll2)
    # #
    #   sum = np.array(df['sum30'][df['sumnz']>0])
    #   nz = np.array(df['sumnz'][df['sumnz'] > 0])
    #   t = np.array( df['tmin'][df['sumnz'] > 0])
    #   scale = np.array(df['scale'][df['sumnz'] > 0])
    # #
    #   f = plt.figure()
    #   # ax = f.add_subplot(111, projection='3d')
    #   # plt.scatter(scale, t, zs=(sum/nz), c='r', marker='o')
    #
    #   ax = f.add_subplot(111)
    #   plt.scatter(t, sum/nz)
    #
    #
  # #

