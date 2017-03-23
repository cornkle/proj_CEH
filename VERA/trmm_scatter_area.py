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
from scipy.stats import gaussian_kde


df = pkl.load(open('/users/global/cornkle/VERA/blobs/trmm_blobs_1000km2.p', 'rb'))


print(df.keys())

p = np.array(df['p'])
pmean = np.array(df['pmean'])
pmax = np.array(df['pmax'])
lon = np.array(df['lon_c'])
lat = np.array(df['lat_c'])
hour = np.array(df['area'])

sahel = (-12,12,15,18)
south = (-12,5,15,11.9)

name = ['Sahel', 'South']

regions = [sahel, south]

dic = {}
dic2 = {}
for n, r in zip(name, regions):

    t = pmean[(lon > r[0]) & (lon < r[2])& (lat > r[1]) & (lat < r[3])]
    tmax = pmax[(lon > r[0]) & (lon < r[2]) & (lat > r[1]) & (lat < r[3])]
    tt = p[(lon > r[0]) & (lon < r[2]) & (lat > r[1]) & (lat < r[3])]

    h = hour[(lon > r[0]) & (lon < r[2])& (lat > r[1]) & (lat < r[3])]
    hcop = tt.copy()
    for id, hh in enumerate(h):
        hcop[id]=hcop[id]*0+hh


    tt = np.concatenate(tt)
    hcop = np.concatenate(hcop)

    #t = np.concatenate(t)
    dic[n] = (t, tmax, h)

    dic2[n] = (tt, hcop)

area = np.arange(0,400001,1000)
centre = area[1::] + (area[0:-1]-area[1::])/2

Sahelcontr = []
Southcontr = []
tmeanSa = []
tmaxSa = []
tmeanSo = []
tmaxSo = []


for id, a in enumerate(area):

    if id == 0:
        continue

    filterSa = (dic['Sahel'][2] <= a) & (dic['Sahel'][2] > area[id-1])
    filterSo = (dic['South'][2] <= a) & (dic['South'][2] > area[id - 1])

    filterSa2 = (dic2['Sahel'][1] <= a) & (dic2['Sahel'][1] > area[id - 1])
    filterSo2 = (dic2['South'][1] <= a) & (dic2['South'][1] > area[id - 1])


    dic['Sahel'][2][filterSa] = area[id-1] + (a-area[id-1])/2
    dic['South'][2][filterSo] = area[id - 1] + (a - area[id - 1]) / 2

    dic2['Sahel'][1][filterSa2] = area[id - 1] + (a - area[id - 1]) / 2
    dic2['South'][1][filterSo2] = area[id - 1] + (a - area[id - 1]) / 2


    sa = np.sum(dic2['Sahel'][0][filterSa2]) / np.sum(dic2['Sahel'][0]) *100
    so = np.sum(dic2['South'][0][filterSo2]) / np.sum(dic2['South'][0]) * 100

    Sahelcontr.append(sa)
    Southcontr.append(so)


print(np.sum(Sahelcontr))

f = plt.figure(figsize=(12, 7), dpi=400)
ax = f.add_subplot(231)


xy = np.vstack([np.sqrt(dic['Sahel'][2]/np.pi), dic['Sahel'][0]])
z = gaussian_kde(xy)(xy)
test = z / (z.max() - z.min())

plt.scatter(np.sqrt(dic['Sahel'][2]/np.pi),dic['Sahel'][0], edgecolor='black', c=test, cmap='viridis' )
plt.title('Sahel: Average rainfall per blob')
plt.xlabel('Radius')
plt.ylabel('Rainfall intensity (mm h-1)')
ax = f.add_subplot(232)


xy = np.vstack([np.sqrt(dic['Sahel'][2]/np.pi), dic['Sahel'][1]])
z = gaussian_kde(xy)(xy)
test = z / (z.max() - z.min())
plt.scatter(np.sqrt(dic['Sahel'][2]/np.pi),dic['Sahel'][1], edgecolor='black', c=test, cmap='viridis' )
plt.title('Sahel: Maximum rainfall per blob')
plt.xlabel('Radius')
plt.ylabel('Rainfall intensity (mm h-1)')
ax1 = f.add_subplot(233)


map = ax1.scatter(np.sqrt(dic2['Sahel'][1]/np.pi),dic2['Sahel'][0], edgecolor='black' )
plt.title('Sahel: All rainfall pixels')
plt.xlabel('Radius')
plt.ylabel('Rainfall intensity (mm h-1)')
ax11 = ax1.twinx()
ax11.scatter(np.sqrt(np.array(centre)/np.pi), Sahelcontr, color='r')
ax11.set_ylabel('Contribution to total rainfall (%)')
#ax11.xlim(0,100000)


ax = f.add_subplot(234)

xy = np.vstack([np.sqrt(dic['South'][2]/np.pi), dic['South'][0]])
z = gaussian_kde(xy)(xy)
test = z / (z.max() - z.min())
plt.scatter(np.sqrt(dic['South'][2]/np.pi),dic['South'][0], edgecolor='black', c=test, cmap='viridis' )
plt.title('South: Average rainfall per blob')
plt.xlabel('Radius')
plt.ylabel('Rainfall intensity (mm h-1)')
ax = f.add_subplot(235)
xy = np.vstack([np.sqrt(dic['South'][2]/np.pi), dic['South'][1]])
z = gaussian_kde(xy)(xy)
test = z / (z.max() - z.min())
plt.scatter(np.sqrt(dic['South'][2]/np.pi),dic['South'][1], edgecolor='black', c=test, cmap='viridis' )
plt.title('South: Maximum rainfall per blob')
plt.xlabel('Radius')
plt.ylabel('Rainfall intensity (mm h-1)')
ax2 = f.add_subplot(236)
map = ax2.scatter(np.sqrt(dic2['South'][1]/np.pi),dic2['South'][0], edgecolor='black' )
plt.title('South: All rainfall pixels ')
plt.xlabel('Radius')
plt.ylabel('Rainfall intensity (mm h-1)')
ax22 = ax2.twinx()
ax22.scatter(np.sqrt(np.array(centre)/np.pi), Southcontr, color='r')
ax22.set_ylabel('Contribution to total rainfall (%)')

plt.legend(fontsize=12)

plt.tight_layout()
plt.savefig('/users/global/cornkle/VERA/plots/precip_scatter_area.png')
plt.close('all')