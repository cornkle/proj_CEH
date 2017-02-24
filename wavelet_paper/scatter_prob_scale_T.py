import seaborn as sns
pal = sns.color_palette('Blues')
sns.set_context("paper", font_scale=1.5)
sns.set_style("ticks")
import ipdb
import matplotlib.cm as cm

from utils import u_statistics as ug

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl


def probability():


    dic = pkl.load(open('/users/global/cornkle/C_paper/wavelet/saves/pandas/3dmax_gt15000_percircle_fake.p', 'rb'))


    ids = np.array(dic['id'])
    scales_all = np.array(dic['scale'])

    udscale = np.unique(scales_all)
    udscale = np.sort(udscale)

    pall = dic['p']
    tmin = np.array(dic['tmin'])


    l_id = []
    l_scale = []
    l_p = []
    l_t = []

    for id, scale, p, t in zip(ids, scales_all, pall, tmin):


        l_id.append(id)
        l_scale.append(scale)
        l_p.append(np.sum(p>38)/np.sum(p>1))
        l_t.append(t)

    sort = np.argsort(np.array(l_p))

    l_p = np.array(l_p)[sort]
    l_scale = np.array(l_scale)[sort]
    l_t = np.array(l_t)[sort]
    l_id = np.array(l_id)[sort]

    f = plt.figure(figsize=(15, 8), dpi=400)
    ax = f.add_subplot(111)

    mappable = ax.scatter(l_scale, l_t, c=l_p, marker="o", s=70, zorder=2, edgecolor='black',
                               linewidth=1, cmap='viridis_r', vmax = 0.5)

    plt.xlabel('Scale (km)')
    plt.ylabel('T_center (degC)')


    cbar = f.colorbar(mappable, label = 'Probability | Rain (mm h-1)>30')
    plt.tight_layout()
    plt.savefig('/users/global/cornkle/C_paper/wavelet/figs/paper/prob_scale_t_fake.png')


def max():
    dic = pkl.load(open('/users/global/cornkle/C_paper/wavelet/saves/pandas/3dmax_gt15000_percircle_fake.p', 'rb'))

    ids = np.array(dic['id'])
    scales_all = np.array(dic['scale'])

    udscale = np.unique(scales_all)
    udscale = np.sort(udscale)

    pall = dic['p']
    tmin = np.array(dic['tmin'])

    l_id = []
    l_scale = []
    l_p = []
    l_t = []

    for id, scale, p, t in zip(ids, scales_all, pall, tmin):
        l_id.append(id)
        l_scale.append(scale)
        l_p.append(np.max(p))
        l_t.append(t)

    sort = np.argsort(np.array(l_p))

    l_p = np.array(l_p)[sort]
    l_scale = np.array(l_scale)[sort]
    l_t = np.array(l_t)[sort]
    l_id = np.array(l_id)[sort]

    f = plt.figure(figsize=(15, 8), dpi=400)
    ax = f.add_subplot(111)

    mappable = ax.scatter(l_scale, l_t, c=l_p, marker="o", s=70, zorder=2, edgecolor='black',
                          linewidth=1, cmap='viridis_r', vmax=90)

    plt.xlabel('Scale (km)')
    plt.ylabel('T_center (degC)')

    cbar = f.colorbar(mappable, label='Max Rain (mm h-1)')
    plt.tight_layout()
    plt.savefig('/users/global/cornkle/C_paper/wavelet/figs/paper/prob_scale_t_max_fake.png')


def frequency():
    dic = pkl.load(open('/users/global/cornkle/C_paper/wavelet/saves/pandas/3dmax_gt15000_percircle.p', 'rb'))

    ids = np.array(dic['id'])
    scales_all = np.array(dic['scale'])

    udscale = np.unique(scales_all)
    udscale = np.sort(udscale)

    pall = dic['p']
    tmin = np.array(dic['tmin'])

    l_id = []
    l_scale = []
    l_p = []
    l_t = []

    for id, scale, p, t in zip(ids, scales_all, pall, tmin):
        l_id.append(id)
        l_scale.append(scale)
        l_p.append(np.max(p))
        l_t.append(t)

    sort = np.argsort(np.array(l_p))

    l_p = np.array(l_p)[sort]
    l_scale = np.array(l_scale)[sort]
    l_t = np.array(l_t)[sort]
    l_id = np.array(l_id)[sort]

    f = plt.figure(figsize=(15, 8), dpi=400)
    ax = f.add_subplot(111)

    mappable = ax.scatter(l_scale, l_t, c=l_p, marker="o", s=70, zorder=2, edgecolor='black',
                          linewidth=1, cmap='viridis_r', vmax=90)

    plt.xlabel('Scale (km)')
    plt.ylabel('T_center (degC)')

    cbar = f.colorbar(mappable, label='Max Rain (mm h-1)')
    plt.tight_layout()
    plt.savefig('/users/global/cornkle/C_paper/wavelet/figs/paper/prob_scale_t_max_fake.png')


if __name__ == "__main__":
    max()