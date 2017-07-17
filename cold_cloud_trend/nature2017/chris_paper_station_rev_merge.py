
# coding: utf-8

import pickle as pkl

import seaborn as sns

pal = sns.color_palette('Blues')
sns.set_context("paper", font_scale=1.5)
sns.set_style("ticks")
#sns.set(color_codes=True)
from scipy.stats import gaussian_kde
from cold_cloud_trend.nature2017 import stations_cc
import numpy as np
from scipy.stats.stats import pearsonr
import matplotlib.pyplot as plt



# In[2]:

dic40 = stations_cc.readStation(40)
dic60 = stations_cc.readStation(60)
dic70 = stations_cc.readStation(70)


# In[3]:

alog = dic40['area']
log=np.log10(alog)
drain = dic40['rain']
drain25 = drain[dic40['area']>25000]
Tmin25 = dic40['minT'][dic40['area']>25000]

thresh = np.squeeze(dic40['thresh'])
rthresh = np.array([thresh]*500).transpose()
nthresh = rthresh[dic40['area']>25000]

ok = np.where(~np.isnan(drain))
drain=drain[ok]
alog = alog[ok]
log = log[ok]
dthresh = rthresh.copy()[ok]

nok = np.where(np.isfinite(log))
drain=drain[nok]
alog = alog[nok]
#log[nok]=0.
log = log[nok]
dthresh = dthresh[nok]


# In[4]:

myDicts = pkl.load( open ('/users/global/cornkle/data/OBS/MSG_TRMM_temp_pcp_300px2004-2013_new.p', 'rb'))


# In[5]:

dic=myDicts[0]
p=np.array(dic['pp']) # 98th perc per MCS
t=np.array(dic['t'])  #mean T
print('Number MCSs:', p.size)
print('98th perc:', np.nanmean(np.array(dic['perc'])))
pc = p[p<81]
tm = t[p<81]


# In[6]:

r=pearsonr(pc,tm)
print(r[0]**2)




# In[7]:

i=np.percentile(Tmin25, np.arange(0,101,10) )

tprob = []
i=np.array(list(reversed(i)))
for cnt , ii, in enumerate(i):

    if cnt == 0:
        continue

    print(i[cnt-1], ii)
    pos = np.where((Tmin25 < i[cnt-1] ) & (Tmin25 >= ii ))
    nex= np.sum(drain25>nthresh)
    ptresh = drain25[pos]
    ttresh = nthresh[pos]
    prob=np.sum(ptresh>ttresh)/nex#/ptresh.size

    tprob.append(prob*100.)
tprob = np.round(np.array(tprob), decimals=1)
print(tprob)
xtick = i[0:-1]

xtickwidth= (i[1::]-i[0:-1])

###############
i=np.percentile(log, np.arange(0,101,10) )
atprob = []

for cnt , ii, in enumerate(i):

    if cnt == 0:
        continue

    pos = np.where((log > i[cnt-1] ) & (log <= ii ))
    nex= np.sum(drain>dthresh)
    ptresh = drain[pos]
    ttresh = dthresh[pos]
    prob=np.sum(ptresh>ttresh)/nex#/ptresh.size

    atprob.append(prob*100.)

atprob = np.round(np.array(atprob), decimals=1)

axtick = i[0:-1]
atickwidth= (i[1::]-i[0:-1])


pos = np.where(drain>dthresh)
aarea = alog[pos]
aprob = aarea[aarea>=25000].size/aarea.size
print('gt 25000', aprob*100.)

#####################
i=np.percentile(tm, np.arange(0,101,10) )
print(i)
ttprob = []
i=np.array(list(reversed(i)))
for cnt , ii, in enumerate(i):

    if cnt == 0:
        continue

    pos = np.where((tm < i[cnt-1] ) & (tm >= ii ))
    nex= np.sum(pc>=38)
    ptresh = pc[pos]

    prob=np.sum(ptresh>=38)/nex#/ptresh.size

    ttprob.append(prob*100.)


ttprob = np.round(np.array(ttprob), decimals=1)

txtick = i[0:-1]
tickwidth= (i[1::]-i[0:-1])



# In[27]:

path = '/users/global/cornkle/C_paper/chris2016/figs/'
fig = plt.figure(figsize=(15, 6), dpi=400)
# definitions for the axes
left, width = 0.061, 0.20
bottom, height = 0.11, 0.7

rect_scatter1 = [left, bottom, width, height]
rect_scatter2 = [left+width+left+left, bottom, width, height]
rect_scatter3 = [5*left+2*width, bottom, width, height]

ax1 = plt.axes(rect_scatter1)
ax11 = ax1.twinx()

ax11.set_ylim((0, 30))
#ax33.set_xlim((-85, -40))

ax11.bar(axtick, atprob, atickwidth, fill=False, color='black', linewidth=1, alpha = 0.8)
plt.text(1, 70, '> 25,000 $km^{2}$: 96%', fontsize=10)

xy = np.vstack([log, drain])
z = gaussian_kde(xy)(xy)
test = z / (z.max() - z.min())
mappable = ax1.scatter(log, drain, c=test, edgecolor='', cmap='viridis_r')
c = 0.8
ax1.set_ylabel('Daily rain ($mm$)')
ax1.set_xlabel('log(Cold cloud area)($km^{2}$) -40$^{\degree}C$')
ax1.vlines(np.log10(25000), np.min(drain), np.max(drain), linestyles='dashed', label='25,000 $km^{2}$', linewidth=c)
ax1.text(4.5, 123, '25,000 $km^{2}$',  fontsize=10, rotation=90)

ax1.hlines(np.mean(thresh), np.min(log), np.max(log), linestyles='dashed', label='threshold', linewidth=c)
# #ax1.set_title('$98^{th}$ perc Pcp  vs meanT per MCS')
# #cbar = fig.colorbar(mappable)
# #cbar.set_label('Normalised probability density')
plt.text(0.03, 0.9, 'a', transform=ax1.transAxes, fontsize=20)
# #text(0.5, 0.35, 'Probability of rain > 38 $mm$:', transform=ax1.transAxes, fontsize=10)
# #text(0.5, 0.3, 'Area > 25,000 $km^{2}$: 17.5%', transform=ax1.transAxes, fontsize=10)

ax1.text(2, 41, '38 $mm$',  fontsize=10)






######################
ax2 = plt.axes(rect_scatter2)
ax22 = ax2.twinx()

ax22.set_ylim((0, 30))
#ax33.set_xlim((-85, -40))

ax22.bar(xtick, tprob, xtickwidth, fill=False, color='black', linewidth=1, alpha = 0.8)

xy = np.vstack([Tmin25, drain25])
z = gaussian_kde(xy)(xy)
test = z / (z.max() - z.min())
mappable = ax2.scatter(Tmin25, drain25,  c=test, edgecolor='', cmap='viridis_r')
ax2.set_ylabel('Daily rain ($mm$)')
ax2.set_xlabel('Local daily min. temperature ($^{\degree}C$)')
ax2.hlines(np.mean(thresh), np.min(Tmin25), np.max(Tmin25), linestyles='dashed', label='threshold', linewidth=c)
#ax1.set_title('$98^{th}$ perc Pcp  vs meanT per MCS')
#cbar = fig.colorbar(mappable)
#cbar.set_label('Normalised probability density')

ax2.text(-42, 41, '38 $mm$',  fontsize=10)
plt.text(0.03, 0.9, 'b', transform=ax2.transAxes, fontsize=20)
plt.text(0.6, 0.92, '< -70$^{\degree}C$: 87%', fontsize=10, transform=ax22.transAxes,)


#     #####################################

ax3 = plt.axes(rect_scatter3)
ax3.set_ylim((0,80))
ax33 = ax3.twinx()

ax33.set_ylim((0, 30))
#ax33.set_xlim((-85, -40))

ax33.bar(txtick, ttprob, tickwidth,  fill=False, color='black', linewidth=1, alpha = 0.8)
ax33.text(-60, 28, '< -55$^{\degree}C$: 81%', fontsize=10)

xy = np.vstack([tm,pc])
z = gaussian_kde(xy)(xy)
test = z / (z.max() - z.min())
mappable = ax3.scatter(tm, pc, c=test, edgecolor='', cmap='viridis_r')

ax3.set_ylabel('Rain ($mm\ h^{-1}$)')
ax3.set_xlabel('Temperature ($^{\degree}C$)')
#ax1.set_title('$98^{th}$ perc Pcp  vs meanT per MCS')
plt.text(0.03, 0.9, 'c', transform=ax3.transAxes, fontsize=20)
ax3.hlines(30, np.min(tm), np.max(tm), linestyles='dashed', label='threshold', linewidth=c)
ax3.text(-41, 31.5, '30 $mm$',  fontsize=10)

#
# from mpl_toolkits.axes_grid1 import make_axes_locatable
# divider = make_axes_locatable(plt.gca())
# cax = divider.append_axes("right", "5%", pad="50%")
# plt.colorbar(mappable, cax=cax, label='Normalised kernel density')
#
#
#



#plt.tight_layout()
plt.savefig(path+'station_TRMM_scatter_merge.png')
#plt.savefig(path+'station_TRMM_scatter.eps')

plt.close('all')

