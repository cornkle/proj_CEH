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
import statsmodels.stats.weightstats as stats
import pickle as pkl

df = pd.read_pickle('/users/global/cornkle/C_paper/wavelet/saves/pandas/3dmax_gt15000_no.p')

ids = np.array(df['id'])
scales_all = np.array(df['scale'])

p30 = np.array(df['circle_g30'])
p30 = np.array(df['circle_sum'])
nz = np.array(df['circle_nz'])
t = np.array(df['bulk_tmean_p'])
tmin = np.array(df['circle_Tcentre'])

max = []
id = []
uids = np.unique(ids)
uscales = np.unique(scales_all)
for i in uids:
    if np.min(t[(ids==i)]) > -50:
        continue

    try:
        maxi = np.nansum(p30[(ids==i)])/np.nansum(nz[(ids==i)])
    # maxi = np.nansum(p30[(ids == i)]) / np.nansum(nz[(ids == i)])
    #     if np.nansum(nz[(ids==i)]) < 5:
    #         continue

    except ValueError:
        continue
    max.append(maxi)
    id.append(i)

id_sort = np.array(id)[np.argsort(max)]

weak = id_sort[0:500]
strong = id_sort[-500:]


weak_scales = []
strong_scales = []
mean_scales = []
wupper =[]
wlower = []
supper = []
slower = []

for s in uscales:
    sumup_weak = []
    sumup_strong = []
    sumup_mean = []
    for i in weak:

        sumup_weak.append(tmin[(ids == i) & (scales_all == s) & (tmin <=-50)])


    for i in strong:
       # print(scales_all[(ids == i)])
        sumup_strong.append( tmin[(ids == i) & (scales_all == s) & (tmin <=-50)])

    for i in uids:
       # print(scales_all[(ids == i)])
        sumup_mean.append( tmin[(ids == i) & (scales_all == s) & (tmin <=-50)])

    weak_scales.append(np.nanmean(np.concatenate(sumup_weak)))
    strong_scales.append(np.nanmean(np.concatenate(sumup_strong)))
    mean_scales.append(np.nanmean(np.concatenate(sumup_mean)))

    wupper.append(stats.zconfint(np.concatenate(sumup_weak))[1])
    wlower.append(stats.zconfint(np.concatenate(sumup_weak))[0])
    supper.append(stats.zconfint(np.concatenate(sumup_strong))[1])
    slower.append(stats.zconfint(np.concatenate(sumup_strong))[0])


ipdb.set_trace()
f = plt.figure()
plt.plot(uscales, weak_scales, label = 'Lowest Probability')
plt.fill_between(uscales, wlower, wupper, alpha=0.3)
#plt.errorbar(uscales, weak_scales, xerr = w_std*2)
plt.plot(uscales, strong_scales, color='r', label = 'Highest probability')
plt.fill_between(uscales, slower, supper, color='r', alpha=0.3)
plt.plot(uscales, mean_scales, color='g', label = 'Average distribution')
plt.xlabel('Scales (km)')
plt.ylabel('Tmean(power max)')
plt.title('MCS > 15000km2, sub-cloud features')
#plt.errorbar(uscales, strong_scales, xerr = s_std*2)

plt.legend()
#f = plt.figure()
# ax1 = f.add_subplot(121)
# ax2 = f.add_subplot(122)
# #ax3 = f.add_subplot(223)
# #ax4 = f.add_subplot(224)
#
# ax1.scatter(scales_all, p30)


