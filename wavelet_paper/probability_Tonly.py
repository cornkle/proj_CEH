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


dic = pkl.load(open('/users/global/cornkle/C_paper/wavelet/saves/pandas/3dmax_gt15000_T.p', 'rb'))
dic2 = pkl.load(open('/users/global/cornkle/C_paper/wavelet/saves/pandas/3dmax_gt15000.p', 'rb'))

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


psum = np.concatenate(np.array(dic['circle_p'])) #[(hour>15) & (hour<23)]
tmin = np.concatenate(np.array(dic['circle_t']))

psum2 = np.concatenate(np.array(dic2['circle_p'])[(scales2<=210)   ])
tmin2 = np.concatenate(np.array(dic2['circle_t'])[(scales2<=210)   ])

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

bins = np.array(list(range(-95, -39, 5)))  # compute probability per temperature range (1degC)
print(bins)


path = '/users/global/cornkle/C_paper/wavelet/figs/paper/'
fig = plt.figure(figsize=(15, 10), dpi=400)
cc = 0.8
width = 0.7 * (bins[1] - bins[0])

center = (bins[:-1] + bins[1:]) / 2

ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)


print('Tbins', bins)
H1, bins1 = np.histogram(tmin[(psum>=30) ], bins=bins, range=(-95, -40))
H, bins = np.histogram(tmin[(psum>=0) ], bins=bins, range=(-95, -40))

H = H.astype(float)
H1 = H1.astype(float)

H[H < 10] = np.nan

histo = H1 / H * 100.

H12, bins12 = np.histogram(tmin2[psum2>=30], bins=bins, range=(-95, -40))
H2, bins2 = np.histogram(tmin2[psum2>=0], bins=bins, range=(-95, -40))
H2 = H2.astype(float)
H12 = H12.astype(float)

H2[H2 < 10] = np.nan

histo2 = H12 / H2 * 100.

lower, upper = stats.proportion_confint(H1, H)
lower2, upper2 = stats.proportion_confint(H12, H2)


ax1.plot(center, histo,  linewidth=1.5 , marker='o', label='temperature only')
ax1.plot(center, histo2,  linewidth=1.5 , marker='o', color='r', label='scale < 40km')
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
plt.savefig(path + 'wavelet_scale_p.png')
# plt.savefig(path + 'wavelet_scale_p_T.pdf')
plt.close('all')
