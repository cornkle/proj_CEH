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

out = '/users/global/cornkle/C_paper/wavelet/saves/pandas/'
df = pkl.load(open(out + '3dmax_gt15000_percircle_fake.p','rb'))

ids = np.array(df['id'])
scales_all = np.array(df['scale'])

pall = np.array(df['p'])
tmin = np.array(df['tmin'])

max = []
id = []
uids = np.unique(ids)
uscales = np.unique(scales_all)

cnt = 0
for i in uids:
    try:
        pp = pall[(ids==i) & (tmin<=-50)].flatten()
        try:
            pp = np.concatenate(pp)
        except ValueError:
            pass

        #maxi = np.nanmax
        maxi = np.sum(pp>30)/np.sum(pp>1)
        cnt = cnt+1
        print(cnt)
    except ValueError:#ValueError:

        continue

    max.append(maxi)
    id.append(i)

f= plt.figure()
plt.plot(np.sort(max))
plt.title('Sorted probability for extreme rain')
plt.ylabel('Probability Rain >30mm')
plt.xlabel('System number')
plt.show()

print('Nb storms: ', len(uids))
id_sort = np.array(id)[np.argsort(max)]

weak = id_sort[0:300]
strong = id_sort[-300:]

weak_scales = []
strong_scales = []
mean_scales = []
w_std =[]
s_std = []


for s in uscales:
    sumup_weak = []
    sumup_strong = []
    sumup_mean = []
    for i in weak:
        sumup_weak.append(np.sum((ids == i) & (scales_all == s)))
      #  print(scales_all[(ids == i)])

    for i in strong:
       # print(scales_all[(ids == i)])
        sumup_strong.append( np.sum((ids == i) & (scales_all == s)))

    for i in uids:
        sumup_mean.append(np.sum((ids == i) & (scales_all == s)) )



    weak_scales.append(np.sum(sumup_weak))
    strong_scales.append(np.sum(sumup_strong))
    mean_scales.append(np.sum(sumup_mean))



    w_std.append(np.std(sumup_weak))
    s_std.append(np.std(sumup_strong))

ipdb.set_trace()
f = plt.figure()
plt.plot(uscales, weak_scales/np.sum(weak_scales) ) #- mean_scales/np.sum(mean_scales)
#plt.errorbar(uscales, np.array(weak_scales)/np.sum(weak_scales), yerr = w_std, label='low probability')
plt.title('Composites of 200 systems with high/low probability for extreme rain')
plt.xlabel('Scale')
plt.ylabel('Frequency')
plt.plot(uscales, strong_scales/np.sum(strong_scales), color='r') #- mean_scales/np.sum(mean_scales)
#plt.errorbar(uscales, np.array(strong_scales)/np.sum(strong_scales), yerr = s_std, color='r', label='high probability')
print(np.sum(weak_scales/np.sum(weak_scales)))
plt.legend()

# f = plt.figure()
# ax1 = f.add_subplot(121)
# ax2 = f.add_subplot(122)
# #ax3 = f.add_subplot(223)
# #ax4 = f.add_subplot(224)
#
# ax1.scatter(scales_all, p30)


