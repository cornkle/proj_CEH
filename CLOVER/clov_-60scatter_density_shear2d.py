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

dic = pkl.load( open ('/users/global/cornkle/data/CLOVER/saves/bulk_-40_zeroRain_gt1k_shear_CP4_JJASNORTH.p', 'rb')) #MSG_TRMM_temp_pcp_300px2004-2013_new.p', 'rb'))

cc=0.8

pp = np.array(dic['pmax'])
sh = np.array(dic['shearmin']) * (-1)
qq = np.array(dic['qmax']) * 1000
tt = np.array(dic['tmin'])
month = np.array(dic['month'])
area = np.array(dic['area'])*(4.4**2)
rainy_area = []
for ip, it in zip(dic['t'], dic['p']):
        rainy_area.append(np.sum(r>=5))

rainy_area = np.array(rainy_area)*(4.4**2)
pos = np.where((pp >= 5) & (sh >= 7) &  (sh <= 30) & (month>=9) & (rainy_area<=400000)  & (area>=1000) )   # 5 + 10 look nicest

tt = tt[pos]
pp = pp[pos]
qq = qq[pos]
sh = sh[pos]

pall = (np.array(dic['p']))[pos]
pure_rainy = rainy_area[pos]
rainy_area = np.log10(area[pos])

nbs=7
nbq=7
nba = 6

shearbins = np.percentile(sh[(sh>=5) & (sh<=30)], np.linspace(0,100,nbs)) #np.percentile(p[(p>=7) & (p<=29)], np.arange(0,101,10)) #np.percentile(p[p>=8], np.arange(0,101,10)) #np.linspace(p[p>=8].min(), p[p>=8].max(),nbshear)
qbins = np.percentile(qq[(qq>=16) & (qq<=19)], np.linspace(0,100,nbq))
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

fig = plt.figure(figsize=(12, 6))

ax1 = fig.add_subplot(231)

xy = np.vstack([sh, pp])
z = gaussian_kde(xy)(xy)
test = z / (z.max() - z.min())

r = stats.pearsonr(sh,pp)

mappable = ax1.scatter(sh, pp, c=test, edgecolor='', cmap='viridis_r', s=20) # viridis_r
ax1.set_ylabel('Max. rainfall (mm h$^{-1}$)')
ax1.set_xlabel('Max. u-shear')
ax1.set_title('Pearsonr: '+str(np.round(r[0], decimals=2)))
cbar = fig.colorbar(mappable)
cbar.set_label('Density')
####################################################################################
ax2 = fig.add_subplot(232)

xy = np.vstack([pp,qq])
z = gaussian_kde(xy)(xy)
test = z / (z.max() - z.min())

r = stats.pearsonr(qq,pp)
print(r)
mappable = ax2.scatter(pp, qq, c=test, edgecolor='', cmap='viridis_r', s=15) # viridis_r
ax2.set_ylim(12, 22)
ax2.set_ylabel('max. q')
ax2.set_xlabel('Max. p')
ax2.set_ylim(12, 22)
ax2.set_title('Pearsonr: '+str(np.round(r[0], decimals=2)))
cbar = fig.colorbar(mappable)
cbar.set_label('Density')


###################################################################################

ax3 = fig.add_subplot(233)

xy = np.vstack([qq,sh])
z = gaussian_kde(xy)(xy)
test = z / (z.max() - z.min())

r = stats.pearsonr(sh,qq)
mappable = ax3.scatter(sh, qq, c=test, edgecolor='', cmap='viridis_r', s=15) # viridis_r

ax3.set_ylabel('max. q')
ax3.set_xlabel('Max. u-shear')
ax3.set_ylim(12, 22)
ax3.set_title('Pearsonr: '+str(np.round(r[0], decimals=2)))
cbar = fig.colorbar(mappable)
cbar.set_label('Density')
####################################################################################
ax4 = fig.add_subplot(234)

ax4.bar(xtick, pprob, xtickwidth, align='edge', ec='black')
ax4.set_xlabel('u-shear bins (deciles)')
ax4.set_ylabel('Probability Rainfall > 60mm h-1')


#####################################################################################

pp = np.array(dic['pmax'])
sh = np.array(dic['shearmin']) * (-1)
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
outval = np.zeros((nbs,nbq))

for issh, shl in enumerate(shearbins[0:-1]):
    for isq, qql in enumerate(qbins[0:-1]):

        poss = np.where((sh >= shl) & (sh < shearbins[issh+1]) & (qq>=qql) & (qq < qbins[isq+1]))
        print('bigger than',shl )
        print('smaller than', shearbins[issh+1])

        pallt = np.concatenate(pall[poss])

        try:
            cmean = np.mean(pp[poss]) #np.percentile(pallt[pallt>=10], 95) # np.mean(pp[poss])#
        except IndexError:
            cmean = np.nan

        prob = np.sum(pp[poss]>=60) / np.sum(pp[poss]>=1)


        print('Rainy pixel number', pallt.size)


        outprob[issh,isq] = prob
        outperc[issh,isq] = cmean
        outval[issh,isq] = len(poss[0])


outarea = np.zeros((nbs,nba))

for issh, shl in enumerate(shearbins[0:-1]):

    for isq, qql in enumerate(areabins[0:-1]):

        poss = np.where((sh >= shl) & (sh < shearbins[issh+1]) & (rainy_area>=qql) & (rainy_area < areabins[isq+1]))
        #print('bigger than',shl )
        #print('smaller than', shearbins[issh+1])
        print('bigger than area',qql )
        print('smaller than area', areabins[isq+1])

        try:
            cmean = np.mean(pp[poss])#percentile(pp[poss], 90)
        except IndexError:
            cmean = np.nan

        areaonly =  np.where((rainy_area>=qql) & (rainy_area < areabins[isq+1]))


        outarea[issh,isq] = cmean

plt.figure()
plt.imshow(outperc, cmap='viridis', vmin=30, vmax=80)
plt.show()

X, Y = np.meshgrid(shearbins,areabins/1000)

Zm = ma.masked_where(np.isnan(outarea),outarea)

cmapp = uplot.discrete_cmap(10, base_cmap='RdBu')
ax5 = fig.add_subplot(235)
mappable = ax5.pcolormesh(X, Y, Zm.T, cmap=cmapp, vmin=30, vmax=40) # viridis_r

ax5.set_ylabel('Area (10$^3$km$^2$)')
ax5.set_xlabel('Shear')
cbar = fig.colorbar(mappable, ticks=np.linspace(30,45,11))
cbar.set_label('90th centile')

X, Y = np.meshgrid(shearbins,qbins)
cmapp = uplot.discrete_cmap(10, base_cmap='RdBu')
ax5 = fig.add_subplot(236)
Zm = ma.masked_where(np.isnan(outperc),outperc)
mappable = ax5.pcolormesh(X, Y, Zm.T, cmap=cmapp, vmin=30, vmax=40) # viridis_r

ax5.set_ylabel('Q')
ax5.set_xlabel('Shear')
cbar = fig.colorbar(mappable, ticks=np.linspace(30,45,11))
cbar.set_label('90th centile')


