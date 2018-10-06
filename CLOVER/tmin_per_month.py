import seaborn as sns
import pickle as pkl
pal = sns.color_palette('Blues')
sns.set_context("paper", font_scale=1.5)
sns.set_style("ticks")
from scipy.stats import gaussian_kde
import numpy as np
import matplotlib.pyplot as plt
from utils import u_statistics as ustat
import pdb
import itertools
from utils import u_plot as uplot
import scipy.stats as stats
import numpy.ma as ma

# In[2]:

dic = pkl.load( open ('/users/global/cornkle/data/CLOVER/saves/bulk_-50_zeroRain_gt1k_shear_CP4.p', 'rb')) #MSG_TRMM_temp_pcp_300px2004-2013_new.p', 'rb'))

cc=0.8

pp = np.array(dic['pmax'])
sh = np.array(dic['shearmin']) * (-1)
qq = np.array(dic['qmax']) * 1000
tt = np.array(dic['tmin'])
month = np.array(dic['month'])
area = np.array(dic['area'])

for r in range(3,11):

    pos = np.where((pp >= 1) & (sh >= 7) &  (sh <= 30) & (month==r) )   # 5 + 10 look nicest
    try:
        print(r, np.percentile(tt[pos], 10))
    except IndexError:
        continue



####### Cloud 50th percentile temperature:
#3 -69.8496306552
#4 -71.0854993535
#5 -70.5316552914
#6 -66.3084593063
#8 -63.6200749564
#9 -64.9375465304
#10 -67.3169104604
