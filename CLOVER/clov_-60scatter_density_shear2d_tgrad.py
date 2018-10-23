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
from utils import constants_lappi as cnst


# In[2]:

dic = pkl.load( open (cnst.CLOVER_SAVES + 'bulk_-50_zeroRain_gt1k_shear_CP4.p', 'rb')) #MSG_TRMM_temp_pcp_300px2004-2013_new.p', 'rb'))

cc=15

pp = np.array(dic['pmax'])
sh = np.array(dic['shearmin']) * (-1)
umin = np.array(dic['umin_mid']) * (-1)
umax = np.array(dic['umax_srfc'])
qq = np.array(dic['qmax']) * 1000
tt = np.array(dic['tmin'])
month = np.array(dic['month'])
area = np.array(dic['area'])*(4.4**2)
lat = np.array(dic['clat'])

pos = np.where( (pp >= 1) & ((month>=9) | (month<=5)))# np.where((pp >= 3) & (sh >= 8) &  (sh <= 30) &  (area<=700000) & ((month<=5) | (month>=10)) & (lat>=5) )   # 5 + 10 look nicest

tt = tt[pos]
pp = pp[pos]
qq = qq[pos]
sh = sh[pos]
umin = umin[pos]
umax = umax[pos]
lats = lat[pos]
area = area[pos]

#new_array = np.vstack([tt,pp,qq,shi,sh,umin,umax]).T
struct = { 'qq' : qq,  'umax': umax, 'umin': umin, 'pp' : pp}
pd_frame = pd.DataFrame(data=struct) #,

parray = u_stat.multi_partial_correlation(pd_frame)


struct = {'pp' : pp, 'qq' : qq, 'sh': sh}
pd_frame = pd.DataFrame(data=struct) #,

psarray = u_stat.multi_partial_correlation(pd_frame)

#new_array = np.vstack([tt,pp,qq,shi,sh,umin,umax]).T
struct = { 'qq' : qq,  'umax': umax, 'umin': umin, 'tt' : tt}
pd_frame = pd.DataFrame(data=struct) #,

tarray = u_stat.multi_partial_correlation(pd_frame)


struct = {'tt' : tt, 'qq' : qq, 'sh': sh}
pd_frame = pd.DataFrame(data=struct) #,

tsarray = u_stat.multi_partial_correlation(pd_frame)

print(parray)
print(psarray)
print(tarray)
print(tsarray)

pall = (np.array(dic['p']))[pos]

rainy_area = np.log10(area)

nbs= 7
nbq= 7
nba = 7

shearbins = np.percentile(sh[(sh>=8) & (sh<=30)], np.linspace(0,100,nbs))

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

fig = plt.figure(figsize=(10, 6), dpi=200)

ax1 = fig.add_subplot(221)

xy = np.vstack([sh, tt])
z = gaussian_kde(xy)(xy)
test = z / (z.max() - z.min())

r = u_stat.pcor(sh,tt, qq)

mappable = ax1.scatter(sh, tt, c=test, edgecolor='', cmap='viridis_r', s=20) # viridis_r
ax1.set_ylabel('Min. T ($^{\circ}$C)')
ax1.set_xlabel('Max. zonal wind shear (m s$^{-1}$)')
ax1.set_title('P-corr. shear/T | q removed: '+str(np.round(r[0], decimals=2)), fontsize=cc)
ax1.tick_params(direction='in')
cbar = fig.colorbar(mappable)
cbar.set_label('Density')

####################################################################################
ax2 = fig.add_subplot(222)

xy = np.vstack([qq,tt])
z = gaussian_kde(xy)(xy)
test = z / (z.max() - z.min())

r = u_stat.pcor(qq,tt, sh)
print(r)
mappable = ax2.scatter(qq, tt, c=test, edgecolor='', cmap='viridis_r', s=15) # viridis_r
ax2.set_xlim(14, 20)
ax2.set_ylabel('Min. T ($^{\circ}$C)')
ax2.set_xlabel('Specific humidity (g kg$^{-1}$)')
ax2.set_title('P-corr. q/T | shear removed: '+str(np.round(r[0], decimals=2)), fontsize=cc)
ax2.tick_params(direction='in')
cbar = fig.colorbar(mappable)
cbar.set_label('Density')



####################################################################################
# ax4 = fig.add_subplot(223)
#
# ax4.bar(xtick, pprob, xtickwidth, align='edge', ec='black')
# ax4.set_xlabel('Zonal shear, bins (equally populated)')
# ax4.set_ylabel('Probability Rainfall > 60mm h-1')
# ax4.set_title('')


#####################################################################################

pp = np.array(dic['pmax'])
sh = np.array(dic['shearmin'])* (-1)
qq = np.array(dic['qmax']) * 1000
tt = np.array(dic['tmin'])

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
outt = np.zeros((nbs,nbq))
outval = np.zeros((nbs,nbq))
outshi = np.zeros((nbs,nbq))


for isq, qql in enumerate(qbins[0:-1]):

    for issh, shl in enumerate(shearbins[0:-1]):

        poss = np.where((sh >= shl) & (sh < shearbins[issh+1]) & (qq>=qql) & (qq < qbins[isq+1]))
        print('bigger than',shl )
        print('smaller than', shearbins[issh+1])

        try:
            cmean = np.percentile(pp[poss], 90)
        except IndexError:
            cmean = np.nan

        try:
            tmean = np.percentile(tt[poss], 10)
        except IndexError:
            tmean = np.nan

        prob = np.sum(pp[poss]>=60) / np.sum(pp[poss]>=1)



        outprob[issh,isq] = prob
        outperc[issh,isq] = cmean
        outval[issh,isq] = len(poss[0])
        outt[issh,isq]=tmean


# plt.figure()
# plt.imshow(outperc, cmap='viridis', vmin=30, vmax=80)
# plt.show()
#

X, Y = np.meshgrid(shearbins,qbins)
cmapp = uplot.discrete_cmap(10, base_cmap='RdBu_r')
ax5 = fig.add_subplot(223)
#outperc[outval<30] = np.nan
Zm = ma.masked_where(np.isnan(outt),outt)
mappable = ax5.pcolormesh(X, Y, Zm.T, cmap=cmapp, vmin=-80, vmax=-70) # viridis_r

ax5.set_ylabel('Max. q925hPa')
ax5.set_xlabel('Max. zonal shear (equally populated)')
ax5.set_title('')
cbar = fig.colorbar(mappable, ticks=np.linspace(-80,-70,11)) # ticks=np.linspace(30,45,11)
cbar.set_label('10th centile min. T')


X, Y = np.meshgrid(shearbins,qbins)
cmapp = uplot.discrete_cmap(10, base_cmap='RdBu')
ax6 = fig.add_subplot(224)
#outperc[outval<30] = np.nan
Zm = ma.masked_where(np.isnan(outperc),outperc)
mappable = ax6.pcolormesh(X, Y, Zm.T, cmap=cmapp, vmin=40, vmax=70) # viridis_r

ax6.set_ylabel('Max. q925hPa')
ax6.set_xlabel('Max. zonal shear (equally populated)')
ax6.set_title('')
cbar = fig.colorbar(mappable, ticks=np.linspace(40,70,11)) # ticks=np.linspace(30,45,11)
cbar.set_label('90th centile max. rain')

plt.tight_layout()
plt.show()
