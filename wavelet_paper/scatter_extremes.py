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
from scipy import stats


df = pd.read_pickle('/users/global/cornkle/C_paper/wavelet/saves/pandas/3dmax_gt15000_30km.pkl')

ids = np.array(df['id'])
scales_all = np.array(df['scale'])

p30 = np.array(df['pmax'])

max = []
id = []
uids = np.unique(ids)
uscales = np.unique(scales_all)
for i in uids:
    try:
        maxi = np.max(p30[(ids==i) & (p30 >=5)])

    except ValueError:
        continue
    max.append(maxi)
    id.append(i)


id_sort = np.array(id)[np.argsort(max)]

weak = id_sort[0:200]
strong = id_sort[-200:]

weak_scales = []
strong_scales = []
w_std =[]
s_std = []
for s in uscales:
    sumup_weak = []
    sumup_strong = []
    for i in weak:
        sumup_weak.append(np.sum((df['id'] == i) & (df['scale'] == s)))
        print(df['scale'][(df['id'] == i)])

    for i in strong:
        print(df['scale'][(df['id'] == i)])
        sumup_strong.append( np.sum((df['id'] == i) & (df['scale'] == s)))

    ipdb.set_trace()
    weak_scales.append(np.sum(sumup_weak))
    strong_scales.append(np.sum(sumup_strong))


    w_std.append(np.std(sumup_weak))
    s_std.append(np.std(sumup_strong))


f = plt.figure()
plt.plot(uscales, weak_scales)
plt.errorbar(uscales, weak_scales, xerr = w_std*2)
plt.plot(uscales, strong_scales, color='r')
plt.errorbar(uscales, strong_scales, xerr = s_std*2)
print(np.sum(weak_scales/np.sum(weak_scales)))

# f = plt.figure()
# ax1 = f.add_subplot(121)
# ax2 = f.add_subplot(122)
# #ax3 = f.add_subplot(223)
# #ax4 = f.add_subplot(224)
#
# ax1.scatter(scales_all, p30)


