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
import statsmodels.stats.proportion as stats


def probability():


    dic = pkl.load(open('/users/global/cornkle/C_paper/wavelet/saves/pandas/3dmax_gt15000_shuffle.p', 'rb'))


    ids = np.array(dic['id'])
    scales_all = np.array(dic['scale'])

    udscale = np.unique(scales_all)
    udscale = np.sort(udscale)

    psum = dic['circle_g30']
    pnz = dic['circle_nz']
    p = dic['circle_p']
    tmin = np.array(dic['circle_Tcentre'])

    ipdb.set_trace()
    l_id = []
    l_scale = []
    l_p = []
    l_t = []

    for id, scale, p, pnz, t in zip(ids, scales_all, p, pnz, tmin):


        l_id.append(id)
        l_scale.append(scale)
        l_p.append(np.sum(p>30)/np.sum(p>0.1))
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
    plt.savefig('/users/global/cornkle/C_paper/wavelet/figs/paper/prob_scale_t_shuffle.png')
    plt.close('all')


def max():
    dic = pkl.load(open('/users/global/cornkle/C_paper/wavelet/saves/pandas/3dmax_gt15000.p', 'rb'))


    ids = np.array(dic['id'])
    scales_all = np.array(dic['scale'])

    udscale = np.unique(scales_all)
    udscale = np.sort(udscale)

    pmax = dic['circle_max']
    tmin = np.array(dic['circle_Tcentre'])


    l_id = []
    l_scale = []
    l_p = []
    l_t = []

    for id, scale, p, t in zip(ids, scales_all, pmax, tmin):
        l_id.append(id)
        l_scale.append(scale)
        l_p.append(p)
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
    plt.savefig('/users/global/cornkle/C_paper/wavelet/figs/paper/prob_scale_t_max.png')
    plt.close('all')

def tmin():


    dic = pkl.load(open('/users/global/cornkle/C_paper/wavelet/saves/pandas/3dmax_gt15000_no.p', 'rb'))


    ids = np.array(dic['id'])
    scales = np.array(dic['scale'])

    uids, uinds = np.unique(dic['id'], return_index=True)

    udscale = np.unique(scales)

    print(np.percentile(scales, np.arange(0,101,20)))

    psum = np.array(dic['circle_p'])
    tmin = np.array(dic['circle_Tcentre'])
    tbulk_min = np.array(dic['bulk_tmin_p'])
    tbulk_mean = np.array(dic['bulk_tmean_p'])
    pbulk_max = np.array(dic['bulk_pmax'])
    pbulk_mean = np.array(dic['bulk_pmean'])
    pmax = np.array(dic['circle_max'])
    pbulk_g30 = np.array(dic['bulk_g30'])

    pp = []
    for p in psum:
        pp.extend(p)

    pp = np.array(pp)

    pp15 = []
    for p, s in zip(psum, scales):
        if s <= 40:
            pp15.extend(p)

    pp = np.array(pp)
    pp15= np.array(pp15)

    print('Nb 30mm pixel in scales', np.sum(pp>=30))
    print('Nb 30mm pixel all', np.sum(pbulk_g30[uinds]))
    print('Nb 30mm pixel 15km', np.sum(pp15>=30))
    print('Nb 30mm pixel 15km fraction', np.sum(pp15>=30) / np.sum(pp>=30))
    print('Nb 30mm pixel 15km fraction all', np.sum(pp15 >= 30) / np.sum(pbulk_g30[uinds]))
    print('Nb 30mm pixel fraction all', np.sum(pp >= 30) / np.sum(pbulk_g30[uinds]))

    bins = np.array(list(range(-95, -39, 5)))  # compute probability per temperature range (1degC)
    print(bins)
    ranges = [10, 29, 50, 205]
    outrange = [ 29, 50,  205]
    # #
    # ranges = [15, 30, 60, 202]
    # outrange = [    30, 60, 202]

    path = '/users/global/cornkle/C_paper/wavelet/figs/paper/'
    fig = plt.figure(figsize=(15, 5), dpi=400)
    cc = 0.8
    width = 0.7 * (bins[1] - bins[0])

    center = (bins[:-1] + bins[1:]) / 2

    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)
    colors = cm.viridis_r(np.linspace(0, 1, len(outrange)))

    for id, r in enumerate(ranges):
        if id == 0:
            continue

        c = colors[id-1]
        start = ranges[id-1]

        t = tmin[(scales <= r) & (scales > ranges[id - 1])]
        p = pmax[(scales <= r) & (scales > ranges[id - 1])]

        to30 = t[p>=30]

        # bins = np.percentile(t, np.arange(0,101,5))
        # center = (bins[:-1] + bins[1:]) / 2

        print('Tbins', bins)
        H1, bins1 = np.histogram(to30, bins=bins, range=(-95, -40))
        H, bins = np.histogram(t, bins=bins, range=(-95, -40))
        H = H.astype(float)
        H1 = H1.astype(float)

        #ipdb.set_trace()

        H[H < 10] = np.nan

        histo = H1 / H * 100.


        lower, upper = stats.proportion_confint(H1, H)

        ax1.plot(center, histo, color=c, linewidth=1.5, label=str(start)+'-'+str(r) + ' km',marker='o' )
        ax1.title('Probability Precip>30mm')
        ax1.fill_between(center, lower*100, upper*100, color=c, alpha=0.3)
        ax2.plot(center, H, color=c, linewidth=1.5, label=str(start) + '-' + str(r) + ' km',marker='o')
        ax3.title('Number of rainfall pixel >30mm (nP)')
        #ax2.set_ylim(0,160)

        ax3.plot(center, H1, color=c, linewidth=1.5, label=str(start) + '-' + str(r) + ' km',marker='o')
        ax3.title('Number of rainfall pixel >30mm (nP)')
    tmean=[]
    tmin = []
    tcmean = []
    for iid in uids:

        pos = np.where(ids == iid)




    ax1.set_xlabel('Min. Temperature (5 $^{\degree}C$ bins)')
    ax1.set_ylabel('Probability (% | Max. precip $>$ 30 $mm\ h^{-1}$)')
    plt.text(0.03, 0.9, 'b', transform=ax1.transAxes, fontsize=20)

    plt.legend()
    plt.tight_layout()
    plt.savefig(path + 'wavelet_scale_p_T_no.png')
   # plt.savefig(path + 'wavelet_scale_p_T.pdf')
    plt.close('all')

if __name__ == "__main__":
    tmin()