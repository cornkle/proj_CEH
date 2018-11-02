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
from utils import u_statistics as u_stat
import pandas as pd
from utils import constants as cnst

# In[2]:

dic = pkl.load( open (cnst.CLOVER_SAVES + 'bulk_-50_zeroRain_gt1k_shear_CP4.p', 'rb')) #MSG_TRMM_temp_pcp_300px2004-2013_new.p', 'rb'))

cc=11

pp = np.array(dic['pmax'])
sh = np.array(dic['shearmin']) * (-1)
umin = np.array(dic['umin_mid']) * (-1)
umax = np.array(dic['umax_srfc'])
qq = np.array(dic['qmax']) * 1000
tt = np.array(dic['tmin'])
month = np.array(dic['month'])
area = np.array(dic['area'])*(4.4**2)
lat = np.array(dic['clat'])

dummy = (umax)*qq/(umin)* * 100

# rainy_area= []
# new_area = []
# for isp in dic['p']:
#     all = np.isfinite(dic['t'])
#     num = np.sum(isp>=0.01)
#     rainy_area.append(num)
#     new_area.append(all)
# new_area = np.array(new_area)
# new_area = np.array(new_area)
#
# pos = np.where((rainy_area>=1) &  (pp >= 0.01) & (month==10) )
#
# # plt.figure()
# # plt.hist(umin[pos], bins=10, range=(3,30))
# # plt.title('smaller')
# # plt.show()
# #
# fig = plt.figure()
# ax1 = fig.add_subplot(121)
# plt.scatter(umin[pos], rainy_area[pos]*(4.4**2))
# plt.title('')
# plt.ylabel('rainy area')
# plt.xlabel('u650hPa')
# plt.show()
#
# ax2 = fig.add_subplot(122)
# plt.scatter(umin[pos], area[pos])
# plt.title('')
# plt.ylabel('storm area')
# plt.xlabel('u650hPa')
# plt.show()


#& (dummy<75) & (dummy>-150)
pos = np.where( (pp >= 8)  & (dummy<12) )# np.where((pp >= 3) & (sh >= 8) &  (sh <= 30) &  (area<=700000) & ((month<=5) | (month>=10)) & (lat>=5) )   # 5 + 10 look nicest

tt = tt[pos]
pp = pp[pos]
qq = qq[pos]
shi = sh[pos]
sh = dummy[pos]
umin = umin[pos]
umax = umax[pos]
lats = lat[pos]
area = area[pos]

#new_array = np.vstack([tt,pp,qq,shi,sh,umin,umax]).T
struct = { 'qq' : qq,  'umax': umax, 'umin': umin, 'pp' : pp}
pd_frame = pd.DataFrame(data=struct) #,

parray = u_stat.multi_partial_correlation(pd_frame)


struct = {'pp' : pp, 'qq' : qq, 'sh': shi}
pd_frame = pd.DataFrame(data=struct) #,

psarray = u_stat.multi_partial_correlation(pd_frame)

#new_array = np.vstack([tt,pp,qq,shi,sh,umin,umax]).T
struct = { 'qq' : qq,  'umax': umax, 'umin': umin, 'tt' : tt}
pd_frame = pd.DataFrame(data=struct) #,

tarray = u_stat.multi_partial_correlation(pd_frame)


struct = {'tt' : tt, 'qq' : qq, 'sh': shi}
pd_frame = pd.DataFrame(data=struct) #,

tsarray = u_stat.multi_partial_correlation(pd_frame)

print(parray)
print(psarray)
print(tarray)
print(tsarray)

tr = tt[tt<=-60]
pr = pp[tt<=-60]

print('TTPP', stats.pearsonr(tr,pr))

pdb.set_trace()

# plt.figure()
# plt.hist(umin, bins=10, range=(3,30))
# plt.title('bigger')
# plt.show()
#
# plt.figure()
# plt.scatter(area,pp)
# plt.title('bigger')
# plt.show()

pall = (np.array(dic['p']))[pos]

rainy_area = np.log10(area)

nbs= 7
nbq= 7
nba = 7


shearbins = np.percentile(sh[(sh>=1) & (sh<=9.5)], np.linspace(0,100,nbs))

realshearbins = np.percentile(sh[(sh>=8) & (sh<=30)], np.linspace(0,100,nbs))

#np.percentile(p[(p>=7) & (p<=29)], np.arange(0,101,10)) #np.percentile(p[p>=8], np.arange(0,101,10)) #np.linspace(p[p>=8].min(), p[p>=8].max(),nbshear)
qbins = np.linspace(15.5,18.5, nbq) #np.linspace(16,19, nbq)#np.linspace(16,19, nbq) #np.percentile(qq[(qq>=16) & (qq<=19)], np.linspace(0,100,nbq))
areabins = np.percentile(rainy_area, np.linspace(0,100,nba) ) #np.linspace(50, 5000,nba) #np.linspace(14.5,19, nbq) #np.percentile(qq[(qq>=14.5) & (qq<=19)], np.linspace(0,100,nbq))#np.linspace(14.5,19, nbq)

shearlist = []
pprob = []
plist = []
#ax5 = fig.add_subplot(236)

colours = [ 'b', 'b',  'c', 'c', 'g', 'g', 'r', 'r', 'k', 'k']
for id, c in enumerate(shearbins[0:-1]):
    posi = np.where((sh >= c) & (sh < shearbins[id+1]))
    print('bigger than',c )
    print('smaller than', shearbins[id+1])

    try:
        cmean = np.percentile(pp[posi], 90)
    except IndexError:
        cmean = np.nan

    H, binz = np.histogram(pp[posi], bins=np.arange(1,111,10))

    prob = np.sum(pp[posi]>=60) / np.sum(pp[posi]>=1)
    pprob.append(prob)
    plist.append(cmean)
    shearlist.append(((shearbins[id+1])-c)/2)


xtick = shearbins[0:-1]
xtickwidth= (shearbins[1::]-shearbins[0:-1])

fig = plt.figure(figsize=(25, 45), dpi=70)

ax1 = fig.add_subplot(331)

xy = np.vstack([sh, pp])
z = gaussian_kde(xy)(xy)
test = z / (z.max() - z.min())

r = u_stat.pcor(sh,pp, qq)

mappable = ax1.scatter(sh, pp, c=test, edgecolor='', cmap='viridis_r', s=20) # viridis_r
ax1.set_ylabel('Max. rainfall (mm h$^{-1}$)')
ax1.set_xlabel('u925hPa')
ax1.set_title('P-corr. u925hPa/rain | q removed: '+str(np.round(r[0], decimals=2)), fontsize=cc)
ax1.tick_params(direction='in')
cbar = fig.colorbar(mappable)
cbar.set_label('Density')

print('Partial correlation umin_mid:', u_stat.pcor(umin,pp, qq) )
####################################################################################
ax2 = fig.add_subplot(332)

xy = np.vstack([pp,qq])
z = gaussian_kde(xy)(xy)
test = z / (z.max() - z.min())

r = u_stat.pcor(qq,pp, sh)
print(r)
mappable = ax2.scatter(qq, pp, c=test, edgecolor='', cmap='viridis_r', s=15) # viridis_r
ax2.set_xlim(14, 20)
ax2.set_ylabel('Max. rainfall')
ax2.set_xlabel('max. q925hPa')
ax2.set_title('P-corr. q/rain | u925hPa removed: '+str(np.round(r[0], decimals=2)), fontsize=cc)
ax2.tick_params(direction='in')
cbar = fig.colorbar(mappable)
cbar.set_label('Density')


###################################################################################

ax3 = fig.add_subplot(333)

xy = np.vstack([qq,sh])
z = gaussian_kde(xy)(xy)
test = z / (z.max() - z.min())

r = u_stat.pcor(qq,sh, umin)
mappable = ax3.scatter(sh, qq, c=test, edgecolor='', cmap='viridis_r', s=15) # viridis_r

ax3.set_ylabel('max. q925hPa')
ax3.set_xlabel('Max. u925hPa')
ax3.set_ylim(14, 20)
ax3.tick_params(direction='in')
ax3.set_title('P-corr. q/u925hPa | u600hPa removed: '+str(np.round(r[0], decimals=2)), fontsize=cc)
cbar = fig.colorbar(mappable)
cbar.set_label('Density')

ax3 = fig.add_subplot(334)

xy = np.vstack([qq,umin])
z = gaussian_kde(xy)(xy)
test = z / (z.max() - z.min())
r = u_stat.pcor(qq,umin, sh)

mappable = ax3.scatter(umin, qq, c=test, edgecolor='', cmap='viridis_r', s=15) # viridis_r

ax3.set_ylabel('max. q925hPa')
ax3.set_xlabel('Min. u600hPa')
ax3.set_ylim(14, 20)
ax3.tick_params(direction='in')
ax3.set_title('P-corr. q/u600hPa | u925hPa removed: '+str(np.round(r[0], decimals=2)), fontsize=cc)
cbar = fig.colorbar(mappable)
cbar.set_label('Density')

ax3 = fig.add_subplot(335)

xy = np.vstack([pp,umin])
z = gaussian_kde(xy)(xy)
test = z / (z.max() - z.min())
r = u_stat.pcor(pp,umin, qq)

rr = stats.pearsonr(pp,umin)

mappable = ax3.scatter(umin, pp, c=test, edgecolor='', cmap='viridis_r', s=15) # viridis_r

ax3.set_ylabel('Max. rain')
ax3.set_xlabel('u650hPa')  #'Max wind shear (600hPa-925hPa)'
ax3.tick_params(direction='in')
ax3.set_title('P-corr. u650hPa/rain|u925hPa removed: '+str(np.round(r[0], decimals=2)), fontsize=cc)
cbar = fig.colorbar(mappable)
cbar.set_label('Density')

ax3 = fig.add_subplot(336)

xy = np.vstack([tt,umin])
z = gaussian_kde(xy)(xy)
test = z / (z.max() - z.min())
r = stats.pearsonr(tt,umin)

mappable = ax3.scatter(umin, tt, c=test, edgecolor='', cmap='viridis_r', s=15) # viridis_r

ax3.set_ylabel('Min T')
ax3.set_xlabel('Max. u925hPa')
ax3.set_title('P-corr. u925hPa/minT|q removed: '+str(np.round(r[0], decimals=2)))
cbar = fig.colorbar(mappable)
cbar.set_label('Density')

####################################################################################
ax4 = fig.add_subplot(337)

ax4.bar(xtick, pprob, xtickwidth, align='edge', ec='black')
ax4.set_xlabel('Max. u925hPa bins (equally populated)')
ax4.set_ylabel('Probability Rainfall > 60mm h-1')
ax4.set_title('')


#####################################################################################

pp = np.array(dic['pmax'])
sh = np.array(dic['umax_srfc'])# * (-1)
dummy = np.array(dic['shearmin'])
qq = np.array(dic['qmax']) * 1000
tt = np.array(dic['tmin'])

sh = sh#/(dummy) *-100

sh

tt = tt[pos]
pp = pp[pos]
qq = qq[pos]
sh = sh[pos]

sheardiff = shearbins[0:-1]+((shearbins[1::]-shearbins[0:-1])/2)[0]
qdiff= qbins[0:-1] + ((qbins[1::]-qbins[0:-1])/2)[0]

shlist = []
qlist = []

outprob = np.zeros((nbs,nbq))
outperc = np.zeros((nbs,nbq))
outval = np.zeros((nbs,nbq))
outshi = np.zeros((nbs,nbq))

corrlist = []
stdlist = []

for isq, qql in enumerate(qbins[0:-1]):

    possi = np.where((qq >= qql) & (qq < qbins[isq + 1]))
    r, p = u_stat.pcor(pp[possi], sh[possi], qq[possi])

    # if p> 0.05:
    #     r = np.nan
    corrlist.append(r)

    # plt.figure()
    # plt.hist(pp[possi],range=(1,100))
    # plt.ylim(0,110)
    # plt.show()

    for issh, shl in enumerate(shearbins[0:-1]):

        poss = np.where((sh >= shl) & (sh < shearbins[issh+1]) & (qq>=qql) & (qq < qbins[isq+1]))
        print('bigger than',shl )
        print('smaller than', shearbins[issh+1])

        pallt = np.concatenate(pall[poss])



        try:
            cmean = np.percentile(pp[poss], 90)
        except IndexError:
            cmean = np.nan

        prob = np.sum(pp[poss]>=60) / np.sum(pp[poss]>=1)


        print('Rainy pixel number', pallt.size)


        outprob[issh,isq] = prob
        outperc[issh,isq] = cmean
        outval[issh,isq] = len(poss[0])

    for issh, shl in enumerate(realshearbins[0:-1]):

        poss = np.where((shi >= shl) & (shi < realshearbins[issh + 1]) & (qq >= qql) & (qq < qbins[isq + 1]))

        try:
            cmean = np.percentile(pp[poss], 90)
        except IndexError:
            cmean = np.nan

        outshi[issh, isq] = cmean

# plt.figure()
# plt.imshow(outperc, cmap='viridis', vmin=30, vmax=80)
# plt.show()

# X, Y = np.meshgrid(shearbins,areabins/1000)
#
# Zm = ma.masked_where(np.isnan(outarea),outarea)
#
# cmapp = uplot.discrete_cmap(10, base_cmap='RdBu')
# ax5 = fig.add_subplot(33)
# mappable = ax5.pcolormesh(X, Y, Zm.T, cmap=cmapp, vmin=30, vmax=40) # viridis_r
#
# ax5.set_ylabel('Area (10$^3$km$^2$)')
# ax5.set_xlabel('Shear')
# cbar = fig.colorbar(mappable, ticks=np.linspace(30,45,11))
# # cbar.set_label('90th centile')

X, Y = np.meshgrid(shearbins,qbins)
cmapp = uplot.discrete_cmap(10, base_cmap='RdBu')
ax5 = fig.add_subplot(339)
#outperc[outval<30] = np.nan
Zm = ma.masked_where(np.isnan(outperc),outperc)
mappable = ax5.pcolormesh(X, Y, Zm.T, cmap=cmapp, vmin=40, vmax=70) # viridis_r

ax5.set_ylabel('Max. q925hPa')
ax5.set_xlabel('Max. u925hPa (equally populated)')
ax5.set_title('')
cbar = fig.colorbar(mappable, ticks=np.linspace(40,70,11)) # ticks=np.linspace(30,45,11)
cbar.set_label('90th centile max. rain')

# ax5 = fig.add_subplot(235)
# Zm = ma.masked_where(np.isnan(outval),outval)
# mappable = ax5.pcolormesh(X, Y, Zm.T, cmap=cmapp)#, vmin=30, vmax=40) # viridis_r
#
# ax5.set_ylabel('Area (10$^3$km$^2$)')
# ax5.set_xlabel('Shear')
# cbar = fig.colorbar(mappable)#, ticks=np.linspace(30,45,11))
# cbar.set_label('90th centile')


xtick = qbins[0:-1]
xtickwidth= (qbins[1::]-qbins[0:-1])

ax5 = fig.add_subplot(338)
#ax5.plot(qbins[0:-1],np.array(corrlist), '-o')
ax5.bar(xtick, np.array(corrlist), xtickwidth, align='edge', ec='black')
ax5.set_ylabel('p-corr (Max. rain/u925hPa | q removed')
ax5.set_xlabel('Max. q925hPa (g/kg) | 0.5 bins')
ax5.set_title('')

plt.show()
