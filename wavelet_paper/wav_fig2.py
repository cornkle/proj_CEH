
import os
import seaborn as sns
import pickle as pkl
pal = sns.color_palette('Blues')
sns.set_context("paper", font_scale=1.5)
sns.set_style("ticks")
import numpy as np
import matplotlib.pyplot as plt
import pdb
from scipy import stats
import matplotlib.cm as cm



# In[173]:

dic = pkl.load( open ('/users/global/cornkle/C_paper/wavelet/saves/bulk_40big_zR.p', 'rb')) #MSG_TRMM_temp_pcp_300px2004-2013_new.p', 'rb'))
dic2 = pkl.load( open ('/users/global/cornkle/C_paper/wavelet/saves/bulk_40big_size_zR.p', 'rb'))

# In[174]:

_p=np.array(dic['pmax']) # 98th perc per MCS
_t=np.array(dic['tmin'])  #mean T
_clat = np.array(dic['clat'])
_area = np.array(dic['area'])*25
_isfin = np.array(dic['isfin'])
_po30 = np.array(dic['po30'])
_perc = np.array(dic['pperc'])
_pp = np.array(dic['p'])
print('Number MCSs:', _p.size)

_p2=np.array(dic2['pmax']) # 98th perc per MCS
_area2 = np.array(dic2['area'])*25
_pp2 = np.array(dic2['p'])


# In[175]:

pthresh = 500
athresh = 5000000
p = _p[(_p<=pthresh)&(_area<=athresh)]
pp = _pp[(_p<=pthresh)&(_area<=athresh)]
t = _t[(_p<=pthresh)&(_area<=athresh)]
clat = _clat[(_p<=pthresh)&(_area<=athresh)]
area = _area[(_p<=pthresh)&(_area<=athresh)]
isfin = _isfin[(_p<=pthresh)&(_area<=athresh)]
po30 = _po30[(_p<=pthresh)&(_area<=athresh)]

p2 = _p2[(_p2<=pthresh)&(_area2<=athresh)]


# In[176]:

print(np.sum(po30))
print('Percentile: ', np.percentile(_perc[_perc>=1], 99))
print(area.min(), area.max())
print(p.min(), p.max())

arsum = np.sum(area)/10
sarea = np.sort(area)


# In[177]:

bins=list(range(0, 900000,5000))   # compute probability per temperature range (1degC) 
bins=[300,1000,2000,5000,10000,15000,25000,50000,100000,250000,500000]
#bins=np.percentile(area,np.arange(0,101,10))
bins=[10000,50000,100000,150000,230000,300000,500000]

apo30=np.where(p > 30)  
area30=area[apo30]   

aH1, abins1 = np.histogram(area30, bins=bins)
aH, abins = np.histogram(area, bins=bins)
aH=aH.astype(float)
aH1=aH1.astype(float)
ahisto=aH1/aH*100.
awidth = 0.7 * (abins[1] - abins[0])
acenter = (abins[:-1] + abins[1:]) / 2
awidth = (abins[1:] - abins[:-1])


# In[178]:

probs=[]
nb = []
cnt = []
for binn in abins:
        prob = np.sum(po30[(area>=binn)])
        nbb = np.sum(area>=binn)
        ccnt = np.sum(isfin[(area>=binn)])
        probs.append(prob)
        nb.append(nbb)
        cnt.append(ccnt)

maxlist =[]
stddev=[]
maxlistt = []
sttddev=[]
arrain = []
for idd, binn in enumerate(bins):
    if idd==0:
        continue
    #pval = np.concatenate(pp[(area>bins[idd-1]) & (area<=binn)])
    pvals=np.concatenate(pp[(area > bins[idd - 1]) & (area <= binn)])
    pval = np.mean((p[(area > bins[idd - 1]) & (area <= binn)]))

    print(np.sum(area[(area > bins[idd - 1]) & (area <= binn)]))

    weights = np.ones_like(pvals[pvals>=0.1]) / float(len(pvals[pvals>=0.1]))
    hist, hplot = np.histogram(pvals[pvals>=0.1], bins=np.arange(30, 150 + 1, 5), range=(30, 150), weights=weights)

    tval = t[(area > bins[idd - 1]) & (area <= binn)]
    imax = np.mean(pval[pval>1])
    std = np.std(pval[pval>1])
    iimax = np.mean(tval[tval<-50])
    tstd = np.std(tval[tval<-50])
    maxlist.append(imax)
    stddev.append(std)
    sttddev.append(tstd)
    maxlistt.append(iimax)
    arrain.append(hist)


thre=0.5
weights = np.ones_like(p[p>=thre]) / float(len(p[p>=thre]))
hist, h = np.histogram(p[p>=thre], bins=np.arange(1,100+1,1), range=(1,100))

weights = np.ones_like(p2[p2>=thre]) / float(len(p2[p2>=thre]))
hist2, h = np.histogram(p2[p2>=thre], bins=np.arange(1,100+1,1),  range=(1,100))



# In[179]:

pprobs = np.array(probs) / np.sum(po30) * 100 # portion of included >30 pixel
nnb = np.array(nb) / len(area) * 100 # portion of systems that size vs all system count
cntt = np.array(probs)/np.array(cnt) # probability of pixel nb > 30 given all pixels in systems

np.sum(po30[area>25000])


# In[180]:

path = '/users/global/cornkle/C_paper/wavelet/figs/paper/'
fig = plt.figure(figsize=(5, 4), dpi=300)
cc=0.8

# ax1 = fig.add_subplot(121)
# mappable = ax1.scatter(acenter, ahisto, c=aH, marker="o",color='#5ea1d4', s=60, zorder=2,edgecolor = 'black', linewidth=1, cmap='viridis_r')
# cbar = fig.colorbar(mappable)
# cbar.set_label('Nb clouds')
# #ax1.bar(abins[0:-1], ahisto, width=awidth, color='#5ea1d4')
#
# ax1.set_xlabel('Area bins (km$^2$) ')
# #ax1.vlines(25000, 0, 100, linestyles='dashed', label='99$^{th}$ percentile', linewidth=1.5, color='red')
# ax1.set_ylabel('Probability (%) | Max. Precip $>$ 30 $mm\ h^{-1}$')
# plt.text(0.03, 0.9, 'a', transform=ax1.transAxes, fontsize=24)
######################

ax1 = fig.add_subplot(111)

ax1.scatter(np.log10(abins), pprobs, marker="o",color='#5ea1d4', s=60, zorder=2, edgecolor = 'black', linewidth=1, label='Rain pixels $>$30 mm h$^{-1}$')
ax1.set_xlabel('log10(Area threshold [km$^2$])')
ax1.set_ylabel('Fraction (%)')

#ax1.set_xticklabels(abins)
#ax2 = ax1.twinx()
ax1.scatter(np.log10(abins), nnb, marker="o",color='red', s=60, zorder=2, edgecolor = 'black', linewidth=1, label='Cold clouds')
#ax2.set_ylabel('Included nb of cold clouds $>$ 30 $mm\ h^{-1}$ (%)') 
#plt.text(0.03, 0.9, 'b', transform=ax1.transAxes, fontsize=24)
ax1.vlines(np.log10(15000), 20, 99, linestyles='dashed',  linewidth=1.5, color='black')
ax1.hlines(32.77, 2.5, 5, linestyles='dashed',  linewidth=1.5, color='black')
ax1.hlines(91.72,2.5, 5, linestyles='dashed',  linewidth=1.5, color='black')
plt.text(4.6, 35, '33 %', fontsize=10)
plt.text(4.6, 94, '91 %', fontsize=10)
plt.text(4.25, 55, '15000 km$^2$', fontsize=10, rotation=90)
ax1.legend(loc='lower left')


plt.tight_layout()
plt.savefig(path+'area-40_test.png')
plt.close('all')

fig = plt.figure(figsize=(8, 5), dpi=300)

ax2 = fig.add_subplot(221)

ax2.scatter(acenter, maxlist)
ax2.set_xlabel('log10(Area threshold [km$^2$])')
ax2.set_ylabel('Average maximum rainfall')
gradient, intercept, r_value, p_value, std_err = stats.linregress(np.log10(acenter),maxlist)
print('rval', r_value)

#ax2.plot(np.log10(acenter), gradient*np.log10(acenter)+intercept)
#ax2.text(3, 50, str(np.round(gradient,decimals=2))+'x + '+str(np.round(intercept,decimals=2)))
#ax2.errorbar(np.log10(acenter), maxlist, yerr=stddev)


# line, = ax2.semilogy(hist, color='r', lw=2, label='Small')
# line, = ax2.semilogy(hist2, color='b', lw=2, label='Big')

ax1 = fig.add_subplot(222)

ax1.scatter(acenter, maxlistt, color='r')
#ax1.errorbar(acenter, maxlistt, yerr=sttddev, color='r')
gradient, intercept, r_value, p_value, std_err = stats.linregress(np.log10(acenter),maxlistt)
print('rval', r_value)
#ax1.plot(np.log10(acenter), gradient*np.log10(acenter)+intercept)
#ax1.text(3, -75, str(np.round(gradient,decimals=2))+'x + '+str(np.round(intercept,decimals=2)))
ax1.set_xlabel('log10(Area threshold [km$^2$])')
ax1.set_ylabel('Average minimum temperature')
ax2.set_title('~ 13000 Cold clouds < -40degC, > 320km2')

ax3 = fig.add_subplot(223)
colors = cm.viridis(np.linspace(0,1,len(arrain)))
for hh, c , b in zip(arrain, colors, bins):
    ax3.semilogy(hplot[0:-1], hh, color=c , label=str(b))

ax3.legend()


plt.tight_layout()
plt.savefig(path+'area_p_t.png')
plt.close('all')


# In[190]:

print(np.sum(isfin[(area>=25000)&(area>=50000)])/25)
print(np.sum(isfin[(area<5000)])/25)


# In[182]:

print('Percentages >30', list(zip(abins, pprobs)))


# In[75]:

print('Nb systems', list(zip(abins, nb)))


# In[76]:

print('Percentages nb systems', list(zip(abins, nnb)))


# In[ ]:

"Only 670 > 25000 left?! -> 940 systems (42% of ) > 10000 deliver 96% of heavy rain"

