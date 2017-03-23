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
import pdb


df = pkl.load(open('/users/global/cornkle/VERA/blobs/trmm_blobs_1000km2.p', 'rb'))
df2 = pkl.load(open('/users/global/cornkle/VERA/blobs/trmm_blobs_all.p', 'rb'))

print(df.keys())

p = np.array(df['p'])
lon = np.array(df['lon_c'])
lat = np.array(df['lat_c'])

p2 = np.array(df2['p'])
lon2 = np.array(df2['lon_c'])
lat2 = np.array(df2['lat_c'])



sahel = (-12,12,15,18)
south = (-12,5,15,11.9)

name = ['Sahel', 'South']

regions = [sahel, south]

dic = {}
dic2 = {}
for n, r in zip(name, regions):

    t = p[(lon > r[0]) & (lon < r[2])& (lat > r[1]) & (lat < r[3])]
    t = np.concatenate(t)
    dic[n] = t

    t2 = p2[(lon2 > r[0]) & (lon2 < r[2]) & (lat2 > r[1]) & (lat2 < r[3])]
    t2 = np.concatenate(t2)
    dic2[n] = t2


f = plt.figure(figsize=(15, 10), dpi=400)
ax = f.add_subplot(111)

colors = cm.rainbow(np.linspace(0,1,len(name)))

for k,c in zip(name, colors):  #[::-1]
    weights = np.ones_like(dic[k]) / float(len(dic[k]))
    hist, h = np.histogram(dic[k], bins=np.arange(0,151,2), weights=weights, range=(0,150))

    line, = ax.semilogy(hist, color=c, lw=2, label=str(k),  marker='o')

    weights = np.ones_like(dic2[k]) / float(len(dic2[k]))
    hist, h = np.histogram(dic2[k], bins=np.arange(0, 151, 2), weights=weights, range=(0, 150))

    line, = ax.semilogy(hist, color=c, lw=2, label=str(k), linestyle='--',   marker='o')

plt.ylabel('ln(normalised frequency of >1mm rain)')
plt.xlim((0,100))
plt.xlabel('Rainfall (mm h-1)')

plt.title('solid: Areas > 1000km2, dashed: all contiguous areas')

plt.legend(fontsize=12)


plt.tight_layout()
plt.savefig('/users/global/cornkle/VERA/plots/precip_histogram_gt1000.png')

f = plt.figure(figsize=(15, 10), dpi=400)
ax = f.add_subplot(111)

colors = cm.rainbow(np.linspace(0,1,len(name)))

for k,c in zip(name, colors):  #[::-1]

    hist, h = np.histogram(dic[k], bins=np.arange(0,151,2),  range=(0,150))

    line, = ax.semilogy(hist, color=c, lw=2, label=str(k), marker='o')

    weights = np.ones_like(dic2[k]) / float(len(dic2[k]))
    hist, h = np.histogram(dic2[k], bins=np.arange(0, 151, 2),  range=(0, 150))


    line, = ax.semilogy(hist, color=c, lw=2, label=str(k), linestyle='--',  marker='o')
plt.ylabel('ln(Count of >1mm rain)')
plt.xlim((0,100))
plt.xlabel('Rainfall (mm h-1)')

plt.title('solid: Areas > 1000km2, dashed: all contiguous areas')

plt.legend(fontsize=12)


plt.tight_layout()
plt.savefig('/users/global/cornkle/VERA/plots/precip_histogram_gt1000_count.png')

plt.close('all')