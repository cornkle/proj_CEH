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
from wavelet_paper import fig9_latitude_var as lv


fpath = '/users/global/cornkle/C_paper/wavelet/figs/paper/'
#path = 'D://data/wavelet/saves/pandas/'
path = '/users/global/cornkle/C_paper/wavelet/saves/pandas/'
dic = pkl.load(open(path+'3dmax_gt15000_lax_nonan_dominant.p', 'rb'))

# bulk = pkl.load( open ('/users/global/cornkle/C_paper/wavelet/saves/bulk_40big_size_zR.p', 'rb'))
#
# print('Sys check', len(bulk['pmax'])), len(np.unique(dic['id']))


psum = np.array(dic['circle_max'])
tmin = np.array(dic['circle_Tcentre'])
# tmin = np.array(dic['circle_t'])
# for id, tt in enumerate(tmin):
#     tmin[id] = np.nanmean(tt)

scales = np.array(dic['scale'])

bins = np.array(list(range(-95, -49, 5)))  # compute probability per temperature range (1degC)
print(bins)
ranges = [10, 60, 90, 180]
outrange = [60, 90, 180]
# #
# ranges = [15, 30, 60, 202]
# outrange = [    30, 60, 202]

fig = plt.figure(figsize=(8, 5), dpi=400)
cc = 0.8
width = 0.7 * (bins[1] - bins[0])

center = (bins[:-1] + bins[1:]) / 2

ax1 = fig.add_subplot(111)

colors = cm.viridis_r(np.linspace(0, 1, len(outrange)))

for id, r in enumerate(ranges):
    if id == 0:
        continue

    c = colors[id - 1]
    start = ranges[id - 1]

    t = tmin[(scales <= r) & (scales > ranges[id - 1])]
    p = psum[(scales <= r) & (scales > ranges[id - 1])]

    to30 = t[p > 30]
    t = t[p>1]

    # bins = np.percentile(t, np.arange(0,101,5))
    # center = (bins[:-1] + bins[1:]) / 2

    print('Tbins', bins)
    H1, bins1 = np.histogram(to30, bins=bins, range=(-95, -50))
    H, bins = np.histogram(t, bins=bins, range=(-95, -50))
    H = H.astype(float)
    H1 = H1.astype(float)

    #H[H < 10] = np.nan

    histo = H1 / H * 100.

    print('Total percentage', r, np.nansum(H1)/ np.nansum(H) * 100.)

    lower, upper = stats.proportion_confint(H1, H)

    ax1.plot(center, histo, color=c, linewidth=1.5, label=str(start) + '-' + str(r) + ' km', marker='o')
    ax1.fill_between(center, lower * 100, upper * 100, color=c, alpha=0.3)

ax1.set_xlabel('SCF T$_{Mean}$ (5 $^{\degree}C$ bins)')
ax1.set_ylabel('SCF proportion with rainfall (%)')
plt.text(0.03, 0.9, 'b', transform=ax1.transAxes, fontsize=20)

plt.legend()
plt.tight_layout()

#plt.savefig(fpath + 'rainy_SCF.png')
# plt.savefig(path + 'wavelet_scale_p_T.pdf')
#plt.close('all')
