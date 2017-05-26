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
from wavelet_paper import latitude_var as lv

def comp_t():
    fpath = '/users/global/cornkle/C_paper/wavelet/figs/paper/'
   # path = 'D://data/wavelet/saves/pandas/'
    path = '/users/global/cornkle/C_paper/wavelet/saves/pandas/'
    dic = pkl.load(open(path+'3dmax_gt15000_TR.p', 'rb'))
    dic2 = pkl.load(open(path+'3dmax_gt15000_noR.p', 'rb'))

    ids = np.array(dic['id'])
    scales = np.array(dic['scale'])
    clat = np.array(dic['clat'])
    hour = np.array(dic['hour'])
    hour2 = np.array(dic2['hour'])

    uids, uinds = np.unique(dic['id'], return_index=True)
    uids2, uinds2 = np.unique(dic2['id'], return_index=True)

    udscale = np.unique(scales)
    pbulk_g30 = np.nansum(np.array(dic['bulk_g30'])[uinds])
    pbulk_g302 = np.nansum(np.array(dic2['bulk_g30'])[uinds2])

    print(np.percentile(scales, np.arange(0, 101, 20)))

    ids2 = np.array(dic2['id'])
    scales2 = np.array(dic2['scale'])
    clat2 = np.array(dic2['clat'])


    psum = np.concatenate(np.array(dic['circle_pc'])) #[(hour>15) & (hour<23)]
    psumm = np.concatenate(np.array(dic['circle_p']))  # [(hour>15) & (hour<23)]
    tmin = np.concatenate(np.array(dic['circle_t']))

    psum2 = np.concatenate(np.array(dic2['circle_pc'])[scales2<=35])
    psumm2 = np.concatenate(np.array(dic2['circle_p'])[scales2 <= 35])
    tmin2 = np.concatenate(np.array(dic2['circle_t'])[scales2<=35])

    print('T', np.sum((psum>=30) ))
    print('S', np.sum((psum2>=30) ))

    print('Tid', np.unique(ids).shape)
    print('Sid', np.unique(ids2).shape)

    pall_g30 = np.sum(psum >= 30)
    pall_g302 = np.sum(np.concatenate(np.array(dic2['circle_p'])[(clat2>=0)]) >= 30)


    pt15 = np.sum((tmin <= -65) & (psum>=30) )
    pp15 = np.sum(psum2>=30)


    print('Nb 30mm bulk T', pbulk_g30)
    print('Nb 30mm bulk S', pbulk_g302)
    print('Nb 30mm identified T', pall_g30)
    print('Nb 30mm identified S', pall_g302)

    print('Nb 30mm identified to bulk T', pall_g30 / pbulk_g30)
    print('Nb 30mm identified to bulk S', pall_g302 / pbulk_g302)
    print('Nb 30mm identified T65', pt15 / pall_g30)
    print('Nb 30mm identified S40', pp15 / pall_g302)

    print('-80 conv rain T', np.sum((tmin <= -80) & (psum>=8) )/np.sum((tmin <= -80) & (psum>=0)))
    print('-80 conv rain S', np.sum((tmin2 <= -80) & (psum2>=8) )/np.sum((tmin2 <= -80)& (psum2>=0)))

    pdb.set_trace()

    print('Number of points', np.sum(np.isfinite(tmin2)))

    bins = np.array(list(range(-95, -44, 5)))  # compute probability per temperature range (1degC)
    print(bins)



    fig = plt.figure(figsize=(15, 10), dpi=400)
    cc = 0.8
    width = 0.7 * (bins[1] - bins[0])

    center = (bins[:-1] + bins[1:]) / 2

    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)


    print('Tbins', bins)
    H1, bins1 = np.histogram(tmin[(psum>=8) ], bins=bins, range=(-95, -45))
    H, bins = np.histogram(tmin[psumm>=0], bins=bins, range=(-95, -45))

    H = H.astype(float)
    H1 = H1.astype(float)

    H[H < 10] = np.nan

    histo = H1 / H * 100.

    H12, bins12 = np.histogram(tmin2[psum2>=30], bins=bins, range=(-95, -45))
    H2, bins2 = np.histogram(tmin2[psumm2>=0], bins=bins, range=(-95, -45))
    H2 = H2.astype(float)
    H12 = H12.astype(float)

    H2[H2 < 10] = np.nan

    histo2 = H12 / H2 * 100.

    lower, upper = stats.proportion_confint(H1, H)
    lower2, upper2 = stats.proportion_confint(H12, H2)


    ax1.plot(center, histo,  linewidth=1.5 , marker='o', label='Temperature only')
    ax1.plot(center, histo2,  linewidth=1.5 , marker='o', color='r', label='Scales < 35km')
    ax1.legend()
    ax1.set_title('Probability Precip>30mm')
    ax1.fill_between(center, lower * 100, upper * 100, alpha=0.3)
    ax1.fill_between(center, lower2 * 100, upper2 * 100, alpha=0.3, color='r')

    print ((histo2-histo)/histo)
    print(histo2-histo)


    ax2.plot(center, H,  linewidth=1.5, marker='o')
    ax2.plot(center, H2,  linewidth=1.5, marker='o', color='r')
    ax2.set_title('Number of valid pixels')
    # ax2.set_ylim(0,160)

    b1 = []
    b2 = []
    # for id, b in enumerate(bins):
    #
    #     if id == 0:
    #         continue
    #
    #     p1 = psum[(tmin>=bins[id-1]) & (tmin<b) & (psum>=30)]
    #     p2 = psum2[(tmin2>=bins[id-1]) & (tmin2<b) & (psum2>=30)]

    ax4.scatter(tmin, psum, color='b')
    ax4.scatter(tmin2, psum2, color='r')

        # b1.append(np.nanmean(p1))
        # b2.append(np.nanmean(p2))


       # ax4.fill_between(center, np.nanmin(p1) , np.nanmax(p1), alpha=0.3)
      #  ax4.fill_between(center, np.nanmin(p2) , np.nanmax(p2), alpha=0.3, color='r')



    ax3.plot(center, H1,  linewidth=1.5,  marker='o')
    ax3.plot(center, H12,  linewidth=1.5,  marker='o', color='r')
    ax3.set_title('Number of rainfall pixel >30mm')
    tmean = []
    tmin = []
    tcmean = []
    for iid in uids:
        pos = np.where(ids == iid)

    ax1.set_xlabel('Min. Temperature (5 $^{\degree}C$ bins)')
    ax1.set_ylabel('Probability (% | Max. precip $>$ 30 $mm\ h^{-1}$)')
    plt.text(0.03, 0.9, 'b', transform=ax1.transAxes, fontsize=20)


    plt.tight_layout()
    plt.savefig(fpath + 'wavelet_scale_p_no.png')

    plt.close('all')

    return center, histo, histo2, upper, lower, upper2, lower2


def plot():
    fpath = '/users/global/cornkle/C_paper/wavelet/figs/paper/'
   # fpath = 'D://data/wavelet/saves/pandas/'
    #path = '/users/global/cornkle/C_paper/wavelet/saves/pandas/'

    center, hist, hist2, up, low, up2, low2 = comp_t()

    cent, prob1, prob2, lower, upper, lower2, upper2, area = lv.comp_lat()


    f = plt.figure(figsize=(12, 5), dpi=400)

    ax1 = f.add_subplot(121)
    ax2 = f.add_subplot(122)

    ax1.plot(center, hist,  linewidth=1.5 , marker='o', label='Temperature only')
    ax1.plot(center, hist2,  linewidth=1.5 , marker='o', color='r', label='Scales$\leq$35km')
    ax1.legend()
    ax1.minorticks_on()
    ax1.set_ylabel('Pixel probability (%)') # | Pixel precip $>$ 30 $mm\ h^{-1}$)')
    ax1.set_xlabel('Pixel temperature (5 $^{\degree}C$ bins)')
    ax1.set_ylim(-1, 85)
    #ax11 = ax1.twinx()
   # ax11.plot(center, (hist2-hist)/hist*100, linestyle='', marker='o', color='grey', label='Percentage change', mec='black', mew=0.5)
   # ax11.set_ylim(-1, 250)
   # ax11.minorticks_on()
    #ax11.set_ylabel('Percentage change (%)')

    ax1.fill_between(center, low * 100, up * 100, alpha=0.3)
    ax1.fill_between(center, low2 * 100, up2 * 100, alpha=0.3, color='r')

    prob1 = np.array(prob1)
    prob2 = np.array(prob2)

    ax2.plot(cent, prob1* 100,  linewidth=1.5 , marker='o', label='Temperature only | -80$^{\degree}C$')
    ax2.plot(cent, prob2* 100,  linewidth=1.5 , marker='o', color='r', label='Scale$\leq$35km | -80$^{\degree}C$')
    ax2.legend()
    ax2.minorticks_on()
    ax2.set_ylabel('  ')# | Pixel precip $>$ 30 $mm\ h^{-1}$)')
    ax2.set_xlabel('Latitude ($^{\degree}N$)')
    ax2.set_ylabel('Pixel probability (%)')
    ax22 = ax2.twinx()
    ax22.plot(cent, np.array(area)/1000, linestyle='', marker='o', color='grey') #(prob2 - prob1) / prob1 * 100
    ax22.set_ylabel('Average MCS area ($10^{3}km^{2}$)')
    ax22.minorticks_on()
    ax2.fill_between(cent, np.array(lower) * 100, np.array(upper) * 100, alpha=0.3)
    ax2.fill_between(cent, np.array(lower2) * 100, np.array(upper2) * 100, alpha=0.3, color='r')


    plt.text(0.03, 0.93, 'a)', transform=ax1.transAxes, fontsize=16)
    plt.text(0.03, 0.93, 'b)', transform=ax2.transAxes, fontsize=16)

    plt.legend()
    plt.tight_layout()
    plt.savefig(fpath + 'tonly_paperR.png')
    # plt.savefig(path + 'wavelet_scale_p_T.pdf')
    plt.close('all')

if __name__ == "__main__":
    plot()