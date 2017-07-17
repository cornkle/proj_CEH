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


    df = pkl.load(open('/users/global/cornkle/C_paper/wavelet/saves/pandas/3dmax_gt15000_0.5.p', 'rb'))

    df2 = pkl.load(open('/users/global/cornkle/C_paper/wavelet/saves/pandas/3dmax_gt15000_0.5.p', 'rb'))

    ids = np.array(df['id'])
    scales = np.array(df['scale'])
    scales2 = np.array(df2['scale'])
    uscales = np.unique(scales)

    tmin = np.array(df['circle_Tcentre'])
    pmax = np.array(df['circle_max'])
    p = np.array(df['circle_p'])

    tmin2 = np.array(df2['circle_Tcentre'])
    pmax2 = np.array(df2['circle_max'])
    p2 = np.array(df2['circle_p'])

    ranges = np.arange(-90,-49,10)
    #ranges=[-90,-40]
    scaler = [15, 20, 30, 40, 50, 60, 70, 80, 100, 120, 150, 205]
    #scaler = np.unique(df['scale'])
    dic = {}
    dic2 = {}
    dic3 = {}
    dic4 = {}

    f = plt.figure(figsize=(15, 8), dpi=400)

   # ax1 = f.add_subplot(221)
   # ax2 = f.add_subplot(222)
    ax3 = f.add_subplot(121)
    ax4 = f.add_subplot(122)
    colors = cm.viridis_r(np.linspace(0, 1, len(ranges)))
    for id, r in enumerate(ranges):
        if id == 0:
            continue
        filter = (tmin <= r) & (tmin > ranges[id - 1]) & (pmax > 0.1)
        filter2 = (tmin2 <= r) & (tmin2 > ranges[id - 1]) & (pmax2 > 0.1)
        sc = (scales[filter])
        sc2 = (scales2[filter2])

        pp = p[filter]
        pp2 = p2[filter2]
        #psum = psum[filter]
        #pnz = pnz[filter]

        dic[r]=[]
        dic2[r] = []
        dic3[r] = []
        dic4[r] = []
        for ids, usc in enumerate(scaler):
            if ids == 0:
                continue
            ffilter = (sc <= usc) & (sc > scaler[ids-1])
            ffilter2 = (sc2 <= usc) & (sc2 > scaler[ids - 1])

            ppf = np.concatenate(pp[ffilter])
            ppf2 = np.concatenate(pp2[ffilter2])

            #dic[r].append(np.nansum(ppf[ppf>30])/np.nansum(ppf>=0))
            #dic[r].append(np.nansum(psum[ffilter])/np.nansum(pnz[ffilter]))


            cnt = 0
            for maxi in pp[ffilter]:
                if np.nanmax(maxi) >= 30:
                    cnt += 1
            #dic[r].append(cnt/len(ffilter)  )

            dic[r].append(np.nansum(ffilter)) #

            dic2[r].append(cnt / np.nansum(ffilter) )

            dic3[r].append(np.nansum(ppf2>30)/np.nansum(ppf2>=0.))
             # dic[r].append(np.nanmax(ppf))

            dic4[r].append(np.percentile(ppf2[ppf2 >= 0.1], 99))


        # ax1.plot(scaler[0:-1], (dic[r]), color=colors[id], label=str(ranges[id-1])+' to '+str(r) + ' C')
        # ax2.plot(scaler[0:-1], (dic2[r]), color=colors[id], label=str(ranges[id - 1]) + ' to ' + str(r) + ' C')
        ax3.plot(scaler[0:-1], (dic3[r]), color=colors[id], label=str(ranges[id - 1]) + ' to ' + str(r) + ' C')
        ax4.plot(scaler[0:-1], (dic4[r]), color=colors[id], label=str(ranges[id - 1]) + ' to ' + str(r) + ' C')

    # ax1.set_xlabel('Scale (km)')
    # ax1.set_ylabel('Probability')
    # #ax1.title("Nb of circles with Rain_max > 30 / nb of circles (per scale/Trange)")
    #
    # ax2.set_xlabel('Scale (km)')
    # ax2.set_ylabel('95th Percentile')

    plt.legend()
    plt.tight_layout()
  #  plt.savefig('/users/global/cornkle/C_paper/wavelet/figs/paper/prob_scatter_0.5.png')
   # plt.close('all')


if __name__ == "__main__":
    probability()