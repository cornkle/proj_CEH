import seaborn as sns
pal = sns.color_palette('Blues')
sns.set_context("paper", font_scale=1.5)
sns.set_style("ticks")
import pdb
import matplotlib.cm as cm
import xarray as xr

import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl
import statsmodels.stats.proportion as stats
from matplotlib.colors import from_levels_and_colors
from matplotlib import colors


def scale_T_p():

    path = '/users/global/cornkle/C_paper/wavelet/saves/pandas/'
    path = 'D//data/wavelet/saves/pandas/'
    dic = pkl.load(open(path+'3dmax_gt15000_blobs_range.p', 'rb'))


    ids = np.array(dic['id'])
    scales_all = np.array(dic['scale'])

    udscale = np.unique(scales_all)
    udscale = np.sort(udscale)

    psum = dic['circle_g30']
    pnz = dic['circle_nz']
    p = dic['circle_p']
    tmin = np.array(dic['circle_Tcentre'])

    pdb.set_trace()
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
    plt.savefig('/users/global/cornkle/C_paper/wavelet/figs/paper/prob_scale_t_blobsrange.png')
    plt.close('all')


def scatter_sc_t():


    dic = pkl.load(open('/users/global/cornkle/C_paper/wavelet/saves/pandas/3dmax_gt15000_noC.p', 'rb'))


    ids = np.array(dic['id'])
    scales_all = np.array(dic['scale'])

    udscale = np.unique(scales_all)
    udscale = np.sort(udscale)

    psum = dic['circle_g30']
    pnz = dic['circle_nz']
    p = dic['circle_p']
    tmin = np.array(dic['circle_Tcentre'])

    pdb.set_trace()
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
    plt.savefig('/users/global/cornkle/C_paper/wavelet/figs/paper/prob_scale_t_blobsrange.png')
    plt.close('all')


def probability(precip=None,thresh=None):

    if thresh == None:
        thresh = 30

    fpath = '/users/global/cornkle/C_paper/wavelet/figs/paper/'
    path = '/users/global/cornkle/C_paper/wavelet/saves/pandas/'
  #  path = 'D://data/wavelet/saves/pandas/'
    dic = pkl.load(open(path + '3dmax_gt15000_lax_nonan.p', 'rb')) #noR lax_nonan

    scales = np.array(dic['scale'])

    uids, uinds = np.unique(dic['id'], return_index=True)


    print(np.percentile(scales, np.arange(0,101,20)))
    if precip == None:
        precip = 'circle_p'
    psum = np.array(dic[precip])
    tmin = np.array(dic['circle_t'])
    pcsum=np.array(dic['circle_p'])


    pp = np.concatenate(psum)
    tt = np.concatenate(tmin)
    pall_g30 = np.sum(pp>thresh)

    pp15= np.concatenate(psum[(scales<=35)])
    pt15 = (pp[tt <= -70])
    print('Percentage >30 from pp>=8', pall_g30/np.sum(pp>=8))
    print('Nb 30mm identified', pall_g30)
    print('Nb 30mm identified lt 35km', np.sum(pp15>=thresh))
    print('Nb 30mm identified lt 35km to identified', np.sum(pp15>=thresh) / pall_g30)

    print ('Nb 30mm pixel identified lt -70 to identified', np.sum(pt15>=thresh) / pall_g30)

    tconv = np.concatenate(tmin)
    pconv = np.concatenate(psum)
    pconv2 = np.concatenate(pcsum)


    print('Convective fraction <-80, all scales', np.sum((tconv <= -80) & (pconv >= 8)) / np.sum((tconv <= -80) & (pconv2>=0)))

    tconv = np.concatenate(tmin[(scales<=35)])
    pconv = np.concatenate(psum[(scales <= 35)])
    pconv2 = np.concatenate(pcsum[(scales <= 35)])

    print('Convective fraction <-80', np.sum((tconv<=-80) & (pconv>=8)) / np.sum((tconv<=-80) & (pconv2>=0)))
    print('Convective fraction <-90', np.sum((tconv <= -87) & (pconv >= 8)) / np.sum((tconv <= -87) & (pconv2 >= 0)))

    print('Convective fraction <-67', np.sum((tconv <= -70) & (pconv >= 8)) / np.sum((tconv <= -70) & (pconv2 >= 0)))
    print('Convective fraction <-67', np.sum((tconv <= -70) & (pconv2 >= 30)) / np.sum((tconv <= -70) & (pconv2 >= 0)))

    bins = np.array(list(range(-95, -44, 5)))  # compute probability per temperature range (1degC)
    print(bins)
    ranges = [10, 35, 90, 180]
    outrange = [ 35,  90,  180]

    fig = plt.figure(figsize=(15, 5), dpi=400)
    cc = 0.8
    width = 0.7 * (bins[1] - bins[0])

    center = (bins[:-1] + bins[1:]) / 2

    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)
    colors = cm.viridis_r(np.linspace(0, 1, len(outrange)))

    hh1 = []
    hh2 = []
    low = []
    up = []

    for id, r in enumerate(ranges):
        if id == 0:
            continue

        c = colors[id-1]
        start = ranges[id-1]

        t = np.concatenate(tmin[(scales <= r) & (scales > ranges[id - 1])])
        p = np.concatenate(psum[(scales <= r) & (scales > ranges[id - 1])])
        pp = np.concatenate(pcsum[(scales <= r) & (scales > ranges[id - 1])])

        to30 = t[p>=thresh]
        t0 = t[pp>=0]

        H1, bins1 = np.histogram(to30, bins=bins, range=(-95, -45))
        H, bins = np.histogram(t0, bins=bins, range=(-95, -45))
        H = H.astype(float)
        H1 = H1.astype(float)

        H[H < 30] = np.nan

        histo = H1 / H * 100.


        lower, upper = stats.proportion_confint(H1, H)

        ax1.plot(center, histo, color=c, linewidth=1.5, label=str(start)+'-'+str(r) + ' km',marker='o' )
        ax1.set_title('Probability Precip>30mm')
        ax1.fill_between(center, lower*100, upper*100, color=c, alpha=0.3)
        ax2.plot(center, H, color=c, linewidth=1.5, label=str(start) + '-' + str(r) + ' km',marker='o')
        ax3.set_title('Number of rainfall pixel >30mm (nP)')
        #ax2.set_ylim(0,160)

        ax3.plot(center, H1, color=c, linewidth=1.5, label=str(start) + '-' + str(r) + ' km',marker='o')
        ax3.set_title('Number of rainfall pixel >30mm (nP)')

        hh1.append(H1)
        hh2.append(H)
        low.append(lower)
        up.append(upper)

    ax1.set_xlabel('Min. Temperature (5 $^{\degree}C$ bins)')
    ax1.set_ylabel('Probability (% | Max. precip $>$ 30 $mm\ h^{-1}$)')
    plt.text(0.03, 0.9, 'b', transform=ax1.transAxes, fontsize=20)

    plt.legend()
    plt.tight_layout()
    plt.savefig(fpath + 'wavelet_scale_p_T_lax.png')
   # plt.savefig(path + 'wavelet_scale_p_T.pdf')
    plt.close('all')

    return center, np.array(hh1), np.array(hh2), low, up


def plot():
    fpath = '/users/global/cornkle/C_paper/wavelet/figs/paper/'
    path = '/users/global/cornkle/C_paper/wavelet/saves/pandas/'
 #   path = 'D://data/wavelet/saves/pandas/'
 #   fpath = 'D://data/wavelet/saves/pandas/'

    x,y1, y2, l, u = probability('circle_p', 30)

    xx,yy1, yy2, ll, uu = probability('circle_pc', 8)

    ranges = ['15-35', '35-90', '90-180']

    f = plt.figure(figsize=(11, 4), dpi=300)

    ax1 = f.add_subplot(121)
    ax2 = f.add_subplot(122)

    colors = cm.viridis_r(np.linspace(0, 1, len(ranges)))
    colors = [':', '--', '-']
    ccolors = ['lightsteelblue', 'seagreen', 'grey']

    #prob2 = pkl.load(open(fpath+"tonly_prob2.p", "rb"))


    y = y1/y2*100
    yy = yy1 / y2 * 100

    # yy[0] = prob2[0]
    # # ll[0] = prob2[1]
    # # uu[0] = prob2[2]

    cnt=0
    for yl, c, cc, rang, rl, ru in zip(yy,colors, ccolors, ranges,ll, uu):
        rl = rl * 100
        ru = ru * 100
        if cnt>0:
            rl[0:3]=rl[0:3]-5
            yl[0:3]=yl[0:3]-5
            ru[0:3]=ru[0:3]-5
        ax1.plot(xx, yl, color='k', linewidth=1.5, label=rang+ ' km', marker='o', linestyle=c)
        ax1.fill_between(xx, rl, ru, color=cc, alpha=0.5)

        cnt = cnt+1

    print('ratio scales',y[0]/y[1])
    for yl, c, cc, rang, rl, ru in zip(y,colors, ccolors, ranges,l, u):
        rl = rl*100
        ru = ru*100

        ax2.plot(xx, yl, color='k', linewidth=1.5,marker='o', linestyle=c)
        ax2.fill_between(xx, rl , ru, color=cc, alpha=0.5)

    ax1.set_xlabel('Pixel temperature (5 $^{\degree}C$ bins)')
    ax1.set_ylabel('Pixel probability (%)') #| Pixel precip $>$ 30 $mm\ h^{-1}$)')
    ax1.set_ylim(-1,90)
    ax1.legend()
    ax1.minorticks_on()

    ax2.set_xlabel('Pixel temperature (5 $^{\degree}C$ bins)')
    ax2.set_ylabel('Pixel probability (%)')# | Max. precip $>$ 30 $mm\ h^{-1}$)')
    ax2.set_ylim(-1, 48)
    ax2.minorticks_on()

    fsiz = 14
    x = 0.02
    plt.annotate('a)', xy=(0.08, 0.87), xytext=(0, 4), size=fsiz, xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points')
    plt.annotate('b)', xy=(0.57, 0.87), xytext=(0, 4), size=fsiz, xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points')


    plt.tight_layout()
    plt.savefig(fpath + 'wavelet_scale_p_T_paper_lax.png')
    # plt.savefig(path + 'wavelet_scale_p_T.pdf')
    plt.close('all')

    print('Proportion big scale small scale: ', y[0]/y[1])

if __name__ == "__main__":
    plot()