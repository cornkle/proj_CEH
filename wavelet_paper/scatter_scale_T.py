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
import scipy.stats as ss
from scipy.stats import gaussian_kde

fpath = '/users/global/cornkle/C_paper/wavelet/figs/paper/'
#path = 'D://data/wavelet/saves/pandas/'
path = '/users/global/cornkle/C_paper/wavelet/saves/pandas/'

dic2 = pkl.load(open(path+'3dmax_gt15000_no.p', 'rb'))

# psum = np.array(dic['circle_p'])[(hour>17) & (hour<=23)] #[(hour>15) & (hour<23)]
# tmin = np.array(dic['circle_t'])[(hour>17) & (hour<=23)]
# lat = np.array(dic['clat'])[(hour>17) & (hour<=23)]

scales = np.array(dic2['scale'])

tmin = np.array(dic2['circle_t'])


t = []
s = []

for ss, tt in zip(scales,tmin):

    dummy = np.zeros_like(tt.flat)+ss

    dummyt = tt.flat

    t.extend(dummyt)
    s.extend(dummy)

t = np.array(t)
s = np.array(s)

print(t.shape, s.shape)

t = t[np.isfinite(t)]
s = s[np.isfinite(t)]

xy = np.vstack([s, t])
z = gaussian_kde(xy)(xy)
test = z / (z.max() - z.min())

f = plt.figure()
plt.plot(t)

f = plt.figure()
plt.scatter(s, t, c=test, cmap = 'viridis')

plt.show()



