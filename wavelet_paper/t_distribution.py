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


df = pkl.load(open('/users/global/cornkle/C_paper/wavelet/saves/pandas/3dmax_gt15000_blobs_range.p', 'rb'))

ids = np.array(df['id'])
scales = np.array(df['scale'])

udscale = np.sort(np.unique(scales))
tmin = np.array(df['circle_Tcentre'])



# ranges = [15, 20, 30, 40, 50, 60, 70, 80, 100, 120, 150, 205 ]
# outrange = [20, 30, 40, 50, 60, 70, 80, 100, 120, 150, 205 ]
ranges = [15, 29, 50, 205]
outrange = [ 29, 50,  205]

dic = {}
for id, r in enumerate(ranges):
    if id == 0:
        continue

    t = tmin[(scales <= r) & (scales>ranges[id-1])]
    dic[r] = t



f = plt.figure(figsize=(15, 7), dpi=400)
ax = f.add_subplot(131)

colors = cm.viridis_r(np.linspace(0,1,len(outrange)))

for id,k in enumerate(outrange):  #[::-1]
    c = colors[id]
    hist, h = np.histogram(dic[k], bins=np.arange(-90,-39,5), range=(-90,-40)) # weights=weights,

    ax.plot(h[1::]-0.5, hist, color=c, lw=2, label=str(ranges[id])+'-'+str(k) + ' km')
    plt.legend(fontsize=7)
    plt.ylabel('Frequency of T(power maximum)')
    plt.xlabel('Tmin per circle')
    plt.title('Sub-system temperature minima, >15000km2')


ax = f.add_subplot(132)

for id,k in enumerate(outrange):  #
    c = colors[id]
    weights = np.ones_like(dic[k]) / float(len(dic[k]))
    hist, h = np.histogram(dic[k], bins=np.arange(-90,-39,5), range=(-90,-40), weights=weights) # weights=weights,

    ax.plot(h[1::]-0.5, hist, color=c, lw=2, label=str(ranges[id])+'-'+str(k) + ' km')
    plt.legend(fontsize=7)
    plt.ylabel('Normalised frequency of T(power maximum)')
    plt.xlabel('Tmin per circle')
    plt.title('Sub-system temperature minima, >15000km2')



arr_list = []
for id,k in enumerate(outrange):  #

    weights = np.ones_like(dic[k]) / float(len(dic[k]))
    hist, h = np.histogram(dic[k], bins=np.arange(-90,-39,5), range=(-90,-40)) # weights=weights,
    arr_list.append(np.array(hist,  dtype=float))

arr = np.sum(np.vstack(arr_list), 0)

arr[arr==0] = np.nan

ax = f.add_subplot(133)
for id,k in enumerate(outrange):
    a = arr_list[id]
    hhist = a / arr
    c = colors[id]
    ax.plot(h[1::]-0.5, hhist, color=c, lw=2, label=str(ranges[id])+'-'+str(k) + ' km')
    plt.legend(fontsize=7)
    plt.ylabel('Normalised frequency of T(power maximum)')
    plt.xlabel('Tmin per circle')
    plt.title('Sub-system temperature minima, >15000km2')



plt.tight_layout()
plt.savefig('/users/global/cornkle/C_paper/wavelet/figs/paper/t_histogram_blobs.png')
plt.close('all')