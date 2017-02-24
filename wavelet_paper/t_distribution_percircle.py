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


dic = pkl.load(open('/users/global/cornkle/C_paper/wavelet/saves/pandas/3dmax_gt15000_percircle.p', 'rb'))

ids = np.array(dic['id'])
scales_all = np.array(dic['scale'])

udscale = np.unique(scales_all)
udscale = np.sort(udscale)
tmin = np.array(dic['tmin'])



dic = {}

for k in udscale:
 #   dic[k] = [item for sublist in dic[k] for item in sublist]
    p = tmin[(scales_all==k)]
    dic[k] = p

f = plt.figure(figsize=(15, 5), dpi=400)
ax = f.add_subplot(131)

colors = cm.viridis(np.linspace(0,1,len(udscale)))

for k,c in zip(udscale, colors):  #[::-1]

    hist, h = np.histogram(dic[k], bins=np.arange(-100,-39,1), range=(-100,-40)) # weights=weights,
    print(h)

    ax.plot(h[1::]-0.5, hist, color=c, lw=2, label=str(k))
    plt.legend(fontsize=7)
    plt.ylabel('Frequency of T(power maximum)')
    plt.xlabel('Tmin per circle')
    plt.title('Sub-system temperature minima, >15000km2')



ax = f.add_subplot(132)

colors = cm.viridis_r(np.linspace(0,1,len(udscale)))

for k,c in zip(udscale[::-1], colors):  #
    weights = np.ones_like(dic[k]) / float(len(dic[k]))
    hist, h = np.histogram(dic[k], bins=np.arange(-100,-39,1), range=(-100,-40), weights=weights) # weights=weights,
    print(h)

    ax.plot(h[1::]-0.5, hist, color=c, lw=2, label=str(k))
    plt.legend(fontsize=7)
    plt.ylabel('Normalised frequency of T(power maximum)')
    plt.xlabel('Tmin per circle')
    plt.title('Sub-system temperature minima, >15000km2')


ax = f.add_subplot(133)


for k, c in zip(udscale[::-1], colors):

    p = dic[k]

    plt.scatter(np.zeros(len(p))+k, p, color='darkseagreen')
    plt.xlabel('Scale')
    plt.ylabel('T(power maximum)')


plt.tight_layout()
plt.savefig('/users/global/cornkle/C_paper/wavelet/figs/paper/t_histogram_scales_count.png')