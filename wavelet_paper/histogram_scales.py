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


df = pkl.load(open('/users/global/cornkle/C_paper/wavelet/saves/pandas/3dmax_gt15000_blobs.p', 'rb'))

p = np.array(df['circle_p'])
ids = np.array(df['id'])
scales = np.array(df['scale'])

udscale = np.sort(np.unique(scales))


ranges = [15, 20, 30, 40, 50, 60, 70, 80, 100, 120, 150, 205 ]
outrange = [20, 30, 40, 50, 60, 70, 80, 100, 120, 150, 205 ]
dic = {}
for id, r in enumerate(ranges):
    if id == 0:
        continue

    t = p[(scales <= r) & (scales>ranges[id-1])]
    t = np.concatenate(t)
    dic[r] = t

for k in outrange:
 #   dic[k] = [item for sublist in dic[k] for item in sublist]
    arr = np.array(dic[k])
    arr = arr[(np.isfinite(arr)) & (arr>=0.1)]

    dic[k] = arr

f = plt.figure(figsize=(15, 8), dpi=400)
ax = f.add_subplot(121)

colors = cm.rainbow(np.linspace(0,1,len(outrange)))

for k,c in zip(outrange[::-1], colors):  #[::-1]
    weights = np.ones_like(dic[k]) / float(len(dic[k]))
    hist, h = np.histogram(dic[k], bins=np.arange(0.1,100+1,1), weights=weights, range=(0.1,100))

    line, = ax.semilogy(hist, color=c, lw=2, label=str(k))
    plt.ylabel('ln(normalised frequency of non-zero rain)')
    plt.xlim((0,120))
    plt.xlabel('rainfall (mm h-1)')
    plt.title('Sub system features of storms, r=scale, >15000km2')

plt.legend(fontsize=8)

ax = f.add_subplot(122)



# colors = cm.rainbow(np.linspace(0, 1, len(keys)))
#
for k, c in zip(outrange[::-1], colors):

    p = dic[k][dic[k] > 0.1]

    plt.scatter(np.zeros(len(p))+k, p)
    plt.xlabel('Scale')
    plt.ylabel('Rain intensity (mm h-1)')
#     weights = np.ones_like(dic[k]) / float(len(dic[k]))
#     hist, h = np.histogram(dic[k], bins=np.arange(0.1, 100 + 1, 1), weights=weights, range=(0.1, 100))
#     print(h)
#
#     plt.plot(hist, color=c, lw=2, label=str(k))
#     plt.ylabel('normalised frequency of non-zero rain')
#     plt.xlabel('rainfall (mm h-1)')

plt.legend()
plt.tight_layout()
plt.savefig('/users/global/cornkle/C_paper/wavelet/figs/paper/precip_histogram_scale_blobs.png')