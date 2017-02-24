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

df = pkl.load(open('/users/global/cornkle/C_paper/wavelet/saves/pandas/3dmax_gt15000_percircle.p', 'rb'))

scales = np.unique(df['scale'])
scales_all = np.array(df['scale'])
id = np.array(df['id'])
p = np.array(df['p'])
tmin = np.array(df['tmin'])
# nbval = np.array(df['sumnz'])
# p30 = np.array(df['sum'])
# area = np.array(df['sumvalid'])

f = plt.figure(figsize=(15, 8), dpi=400)
ax1 = f.add_subplot(121)
ax2 = f.add_subplot(122)
#ax3 = f.add_subplot(223)
#ax4 = f.add_subplot(224)

colors = cm.rainbow(np.linspace(0,1,len(scales)))
x = []
y = []
for k,c in zip(scales[::-1], colors): #

    pos = np.where((scales_all == k) & (tmin <=-60))# )

    dummy = p[(scales_all == k)& (tmin <=-60)]
    ids = id[(scales_all == k)& (tmin <=-60)]
    nb = np.sum((scales_all==k) &  (tmin <=-60))#len(np.unique(ids))

    try:
        dummy2 = np.concatenate(dummy)

    except ValueError:
        pass

    pval = np.sum(dummy2>0.1)/np.float(nb)
    pp30 = np.sum(dummy2>30)/np.float(nb)

    ax1.scatter(pval, pp30, color=c, label=str(k))
   # ax1.set_xlim((0,10000))
   # ax1.set_ylim((0,1000))

    ax1.set_ylabel('average nb pix > 30 mm h-1')
    ax1.set_xlabel('average nb pix > 0.1 mm h-1')

    ax2.scatter(k, pp30/pval, color=c, label=str(k))
    ax1.set_ylabel('nb pix > 30 mm h-1')
    ax1.set_xlabel('nb pix > 0.1 mm h-1')


    x.append(pval)
    y.append(pp30)


x= np.array(x)
y = np.array(y)

gradient, intercept, r_value, p_value, std_err = stats.linregress(x,y)
print('rval', r_value)

ax1.plot(x, gradient*x+intercept)
ax1.text(50, 12, str(np.round(gradient,decimals=2))+'x + '+str(np.round(intercept,decimals=2)))



ax1.legend(fontsize = 9)
ax2.legend(fontsize = 9)

plt.tight_layout()
plt.savefig('/users/global/cornkle/C_paper/wavelet/figs/paper/scatter_scales_percircle_test2_-60.png')
