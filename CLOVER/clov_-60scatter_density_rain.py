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
import scipy.stats as stats
from utils import constants as cnst

# In[2]:

dic = pkl.load( open (cnst.CLOVER_SAVES + 'bulk_-50_zeroRain_gt1k_shear_CP4.p', 'rb')) #MSG_TRMM_temp_pcp_300px2004-2013_new.p', 'rb')

fig = plt.figure(figsize=(12, 6))
cc=0.8

t = np.array(dic['pmax'])
p = (np.array(dic['umax_srfc'])/np.array(dic['shearmin']))*-100
p[np.isinf(p)]=0

month = np.array(dic['month'])
pos = np.where((t >= 1) & (p >= 4))  #

t = t[pos]
p = p[pos]

xy = np.vstack([t, p])
z = gaussian_kde(xy)(xy)
test = z / (z.max() - z.min())

r = stats.pearsonr(t, p)

bins = np.percentile(p[p>=8], np.arange(0,101,10))#np.arange(4, 29, 2)  # compute probability per temperature range (1degC)
shearlist = []
pprob = []
plist = []
#ax5 = fig.add_subplot(236)
plt.figure()
colours = [ 'b', 'b',  'c', 'c', 'g', 'g', 'r', 'r', 'k', 'k']

for id, c in enumerate(bins[0:-1]):
    pos = np.where((p >= c) & (p < bins[id+1]))
    print('bigger than',c * 2 )
    print('smaller than', bins[id+1])

    try:
        cmean = np.percentile(t[pos], 90)
    except IndexError:
        cmean = np.nan

    H, binz = np.histogram(t[pos], bins=np.arange(1,111,10))
    #H, bins = ustat.histo_frequency(t[pos])

    plt.plot(binz[0:-1]+(binz[1::]-binz[0:-1]), H, 'o-',  label=str(c)+'-'+str(bins[id+1]) , color=colours[id])

    prob = np.sum(t[pos]>=60) / np.sum(t[pos]>=1)
    pprob.append(prob)
    plist.append(cmean)
    shearlist.append(((bins[id+1])-c)/2)

plt.xlabel('Maximum rainfall intensity')
plt.ylabel('Number of storms')
plt.plot()
plt.legend()
xtick = bins[0:-1]

xtickwidth= (bins[1::]-bins[0:-1])

# bins=np.arange(4, 29, 2)  # compute probability per temperature range (1degC)
# ppo30=np.where(t > 50)
# to30=p[ppo30]
#
# H1, bins1 = np.histogram(to30, bins=bins, range=(t.min(), t.max()))
# H, bins = np.histogram(t, bins=bins, range=(t.min(), t.max()))
# H=H.astype(float)
# H1=H1.astype(float)
# histo=H1/H*100.
# width = 0.7 * (bins[1] - bins[0])
# center = (bins[:-1] + bins[1:]) / 2
#
# print('Bins included, [bin[', bins[0:11])
# print('Fraction of <-80s to be >30', np.sum(H1[0:11])/np.sum(H[0:11]))
#
# print('Bins included, [bin[', bins[-31::])
# print('Fraction of >-40s to be >30', np.sum(H1[-31::])/np.sum(H[-31::]))
#
# # In[7]:
#
# H, xedges, yedges = np.histogram2d(t, p, bins=(25, 25))
# hh = (H / np.sum(H))*100.
# np.sum(hh)

ax1 = fig.add_subplot(231)

#rarea/1000 # z / (z.max() - z.min()) #sarea #
mappable = ax1.scatter(p, t, c=test, edgecolor='', cmap='viridis_r', s=20) # viridis_r
#ax1.set_title('~13400 contiguous cold clouds (< -10$^{\degree}C$, > 325 km$^2$ )')
ax1.set_ylabel('Max. rainfall (mm h$^{-1}$)')
ax1.set_xlabel('Max. u-shear')
ax1.set_xlim(-10,50)
ax1.set_title('Pearsonr: '+str(np.round(r[0], decimals=2)))
cbar = fig.colorbar(mappable)
cbar.set_label('Density')
#plt.text(0.03, 0.9, 'a', transform=ax1.transAxes, fontsize=20)


ax2 = fig.add_subplot(232)
t = np.array(dic['pmax'])
p = np.array(dic['qmax'])
xy = np.vstack([t,p])
z = gaussian_kde(xy)(xy)
test = z / (z.max() - z.min())

r = stats.pearsonr(t,p)
mappable = ax2.scatter(p, t, c=test, edgecolor='', cmap='viridis_r', s=15) # viridis_r

ax2.set_ylabel('Min. temperature ($^{\degree}C$)')
ax2.set_xlabel('Max. u-shear')
ax2.set_title('Pearsonr: '+str(np.round(r[0], decimals=2)))
cbar = fig.colorbar(mappable)
cbar.set_label('Density')


ax3 = fig.add_subplot(233)

t = np.array(dic['qmax'])
p = np.array(dic['shearmin'])
p[np.isinf(p)]=0

xy = np.vstack([t,p])
z = gaussian_kde(xy)(xy)
test = z / (z.max() - z.min())

r = stats.pearsonr(t,p)
mappable = ax3.scatter(p, t, c=test, edgecolor='', cmap='viridis_r', s=15) # viridis_r

ax3.set_ylabel('max. q')
ax3.set_xlabel('Max. u-shear')
ax3.set_ylim(0.012, 0.022)
ax3.set_title('Pearsonr: '+str(np.round(r[0], decimals=2)))
cbar = fig.colorbar(mappable)
cbar.set_label('Density')

ax4 = fig.add_subplot(234)

ax4.bar(xtick, plist, xtickwidth, align='edge', ec='black')
ax4.set_xlabel('u-shear bins (deciles)')
ax4.set_ylabel('99th centile rain')

ax4 = fig.add_subplot(235)

ax4.bar(xtick, pprob, xtickwidth, align='edge', ec='black')
ax4.set_xlabel('u-shear bins (deciles)')
ax4.set_ylabel('Probability Rainfall > 60mm h-1')


# ax5 = fig.add_subplot(236)
# t = np.array(dic['tmin'])
# p = np.array(dic['pmax'])
#
# pos = np.where((p >= 1))
# t = t[pos]
# p = p[pos]
#
# xy = np.vstack([t,p])
# z = gaussian_kde(xy)(xy)
# test = z / (z.max() - z.min())
#
# r = stats.pearsonr(t,p)
# mappable = ax5.scatter(t, p, c=test, edgecolor='', cmap='viridis_r', s=15) # viridis_r
#
# ax5.set_ylabel('Min. temperature ($^{\degree}C$)')
# ax5.set_xlabel('Max. rain')
# ax5.set_title('Pearsonr: '+str(np.round(r[0], decimals=2)))
# cbar = fig.colorbar(mappable)
# cbar.set_label('Density')
#

plt.tight_layout()
plt.show()
