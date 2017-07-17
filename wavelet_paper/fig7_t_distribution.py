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


#df = pkl.load(open('/users/global/cornkle/C_paper/wavelet/saves/pandas/3dmax_gt15000_no.p', 'rb'))

fpath = '/users/global/cornkle/C_paper/wavelet/figs/paper/'
path = '/users/global/cornkle/C_paper/wavelet/saves/pandas/'
# path = 'D://data/wavelet/saves/pandas/'

df = pkl.load(open(path+'3dmax_gt15000_noR.p', 'rb'))

ids = np.array(df['id'])
scales = np.array(df['scale'])

udscale = np.sort(np.unique(scales))
tmin = np.array(df['circle_Tcentre'])
tmin = np.array(df['circle_t'])
lat = np.array(df['clat'])
tmin2 = np.array(df['circle_t'])

# for id, tt in enumerate(tmin):
#     tmin[id] = np.nanmean(tt)

# ranges = [15, 20, 30, 40, 50, 60, 70, 80, 100, 120, 150, 205 ]
# outrange = [20, 30, 40, 50, 60, 70, 80, 100, 120, 150, 205 ]
ranges = [10, 35, 90, 180]
outrange = [ 35, 90,  180]

dic = {}
dic_l = {}
for id, r in enumerate(ranges):
    if id == 0:
        continue

    llat = lat[(scales <= r) & (scales > ranges[id - 1])]
    t = np.concatenate(tmin[(scales <= r) & (scales>ranges[id-1])])

    print('Number valid pixel', np.sum(np.isfinite(np.concatenate(tmin2[(scales <= r) & (scales > ranges[id - 1])]))))

    dic[r] = t
    dic_l[r] = llat



f = plt.figure(figsize=(10, 4), dpi=300)
# ax = f.add_subplot(131)
#
colors = cm.viridis_r(np.linspace(0,1,len(outrange)))
#
# for id,k in enumerate(outrange):  #[::-1]
#     c = colors[id]
#     hist, h = np.histogram(dic[k], bins=np.arange(-90,-44,3), range=(-90,-45)) # weights=weights,
#
#     ax.plot(h[1::]-0.5, hist, color=c, lw=2, label=str(ranges[id])+'-'+str(k) + ' km', marker='o')
#     plt.legend(fontsize=7)
#     plt.ylabel('Frequency of T(power maximum)')
#     plt.xlabel('Tmin per circle')
#     plt.title('Sub-system temperature minima, >15000km2')

linestyles = [':', '--', '-']
ax = f.add_subplot(121)

for id,k in enumerate(outrange):  #
    c = colors[id]
    ll = linestyles[id]
    weights = np.ones_like(dic[k]) / float(len(dic[k]))
    start = ranges[id]
    if start == 10:
        start = 15

    hist, h = np.histogram(dic[k], bins=np.arange(-95, -44, 5), weights=weights, range=(-95,-45)) # weights=weights,
    print(np.sum(hist))
    print(k)

    ax.plot(h[1::]-0.5, hist, lw=2, label=str(start)+'-'+str(k) + ' km', marker='o', linestyle=ll, color='k')

    plt.legend(fontsize=9, handlelength=2.5, loc='upper left')
    plt.ylabel('Normalised frequency')
    plt.xlabel('Pixel temperature (5 $^{\degree}C$ bins)')
    plt.annotate('a)', xy=(0.04,0.94),xytext=(0,4), size=15, xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points') #transform=ax.transAxes,
    #plt.title('Sub-system temperature minima, >15000km2')



arr_list = []
for id,k in enumerate(outrange):  #

    weights = np.ones_like(dic[k]) / float(len(dic[k]))
    hist, h = np.histogram(dic[k], bins=np.arange(-95, -44, 5), range=(-95,-45)) # weights=weights,
    arr_list.append(np.array(hist,  dtype=float))


arr = np.sum(np.vstack(arr_list), 0)

arr[arr==0] = np.nan
print('-80C share', np.sum(arr_list[0][0:3])/np.sum(arr[0:3]))

ax = f.add_subplot(122)
for id,k in enumerate(outrange):
    a = arr_list[id]
    hhist = a / arr
    c = colors[id]
    ll = linestyles[id]
    start = ranges[id]
    if start == 10:
        start = 15

    ax.plot(h[1::]-0.5, hhist*100, lw=2, label=str(start)+'-'+str(k) + ' km', marker='o', linestyle=ll, color='k')

   # plt.legend(fontsize=8)
    plt.ylabel('Fraction (%)')
    plt.xlabel('Pixel temperature (5 $^{\degree}C$ bins)')
   # plt.title('Sub-system temperature minima, >15000km2')
    plt.annotate('b)', xy=(0.52, 0.94), xytext=(0, 4), size=15, xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points')

# ax = f.add_subplot(133)
# weights = np.ones_like(lat) / float(len(lat))
# histm, h = np.histogram(lat, bins=np.arange(4, 20, 3), weights=weights, range=(4, 20))  # weights=weights,
#
# for id,k in enumerate(outrange):
#     c = colors[id]
#     weights = np.ones_like(dic_l[k]) / float(len(dic_l[k]))
#     hist, h = np.histogram(dic_l[k], bins=np.arange(4, 20, 3), weights=weights, range=(4, 20))  # weights=weights,
#
#     print(np.sum(hist))
#     print(k)
#     start = ranges[id]
#     if start == 10:
#         start = 15
#     ax.plot(h[1::] - 0.5, hist-histm, color=c, lw=2, label=str(start) + '-' + str(k) + ' km', marker='o')
#     plt.legend(fontsize=8)
#     plt.ylabel('Normalised frequency')
#     plt.xlabel('Latitude')
#     plt.annotate('a)', xy=(0.08, 0.94), xytext=(0, 4), size=15, xycoords=('figure fraction', 'figure fraction'),
#                  textcoords='offset points')  # transform=ax.transAxes,
#     # plt.title('Sub-system temperature minima, >15000km2')


plt.tight_layout()
plt.savefig(fpath+'t_histogram.png')
plt.close('all')