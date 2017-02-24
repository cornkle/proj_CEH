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


dic = pkl.load(open('/users/global/cornkle/C_paper/wavelet/saves/pandas/3dmax_gt15000_30km_percircle.p', 'rb'))

dic2 = pkl.load(open('/users/global/cornkle/C_paper/wavelet/saves/bulk_40big.p', 'rb'))

ids = np.array(dic['id'])
scales_all = np.array(dic['scale'])

udscale = np.unique(scales_all)
udscale = np.sort(udscale)

pall = np.array(dic['p'])
tmin = np.array(dic['tmin'])

d2 = np.concatenate(np.array(dic2['p']))


dic = {}

for k in udscale:
 #   dic[k] = [item for sublist in dic[k] for item in sublist]
    p = pall[(scales_all==k) & (tmin <= -70)].flatten()
    try:
        p = np.concatenate(p)
    except ValueError:
        pass

    dic[k] = p

f = plt.figure(figsize=(15, 8), dpi=400)
ax = f.add_subplot(121)

colors = cm.rainbow(np.linspace(0,1,len(udscale)))

for k,c in zip(udscale[::-1], colors):  #[::-1]
    weights = np.ones_like(dic[k]) / float(len(dic[k]))
    hist, h = np.histogram(dic[k], bins=np.arange(0.1,100+1,1), weights=weights, range=(0.1,100))
    print(h)

    line, = ax.semilogy(hist, color=c, lw=2, label=str(k))
    plt.ylabel('ln(normalised frequency of non-zero rain)')
    plt.xlim((0,120))
    plt.xlabel('rainfall (mm h-1)')
    plt.title('Sub system features of storms <-70degC, r=30km, >15000km2')

weights = np.ones_like(d2) / float(len(d2))
hist, h = np.histogram(d2, bins=np.arange(0.1,100+1,1), weights=weights, range=(0.1,100))
print(h)

line, = ax.semilogy(hist, color='black', lw=2)
plt.legend(fontsize=8)

ax = f.add_subplot(122)



# colors = cm.rainbow(np.linspace(0, 1, len(keys)))
#
for k, c in zip(udscale[::-1], colors):

    p = dic[k][dic[k] > 30]

    plt.scatter(np.zeros(len(p))+k, p)
    plt.xlabel('Scale')
    plt.ylabel('Rain intensity (mm h-1)')

plt.legend()
plt.tight_layout()
plt.savefig('/users/global/cornkle/C_paper/wavelet/figs/paper/precip_histogram_scale_-70_30km.png')

#     weights = np.ones_like(dic[k]) / float(len(dic[k]))
#     hist, h = np.histogram(dic[k], bins=np.arange(0.1, 100 + 1, 1), weights=weights, range=(0.1, 100))
#     print(h)
#
#     plt.plot(hist, color=c, lw=2, label=str(k))
#     plt.ylabel('normalised frequency of non-zero rain')
#     plt.xlabel('rainfall (mm h-1)')

# f = plt.figure(figsize=(15, 8), dpi=400)
# ax = f.add_subplot(121)
#
# weights = np.ones_like(d2) / float(len(d2))
# hist, h = np.histogram(d2, bins=np.arange(0.1,100+1,1), weights=weights, range=(0.1,100))
# print(h)
#
# line, = ax.semilogy(hist, color=c, lw=2)
# plt.ylabel('ln(normalised frequency of non-zero rain)')
# plt.xlim((0,120))
# plt.xlabel('rainfall (mm h-1)')
#
#
# ax = f.add_subplot(122)
# p = d2[d2[k] > 30]
#
# plt.plot(p)
# plt.xlabel('Scale')
# plt.ylabel('Rain intensity (mm h-1)')
#
# plt.show()