# import seaborn as sns
import pickle as pkl
# pal = sns.color_palette('Blues')
# sns.set_context("paper", font_scale=1.5)
# sns.set_style("ticks")
from scipy.stats import gaussian_kde
import numpy as np
import matplotlib.pyplot as plt
from utils import u_statistics as ustat
import pdb
import itertools
from utils import u_plot as uplot, constants as cnst
import scipy.stats as stats
import numpy.ma as ma
from utils import u_statistics as u_stat
import pandas as pd
import xarray as xr

def get_ERA(era, indic):

    dic = {}
    dic['u925'] = []
    dic['u650'] = []
    dic['q925'] = []
    dic['q700'] = []
    for id, date in enumerate(indic.date):

        getera =np.where((era['time.day']==(indic['date'])[id].day) & (era['time.month']==indic.month[id]) & (era['time.year']==indic.year[id]))

        try:
            era_day = era.isel(time=int(getera[0]))
        except TypeError:
            print('Era missing')
            dic['u925'].append(np.nan)
            dic['u650'].append(np.nan)
            dic['q925'].append(np.nan)
            dic['q700'].append(np.nan)

            continue

        dic['u925'].append(float(era_day['u'].sel(latitude=elat, longitude=elon, level=925, method='nearest').values))
        dic['u650'].append(float(era_day['u'].sel(latitude=elat, longitude=elon, level=650, method='nearest').values))
        dic['q925'].append(float(era_day['q'].sel(latitude=elat, longitude=elon, level=925, method='nearest').values))
        dic['q700'].append(float(era_day['q'].sel(latitude=elat, longitude=elon, level=700, method='nearest').values))

    return dic


#dic = pkl.load( open (cnst.CLOVER_SAVES + 'bulk_-40_zeroRain_gt5k_-40thresh_OBSera.p', 'rb')) #MSG_TRMM_temp_pcp_300px2004-2013_new.p', 'rb'))
pdf = pkl.load(open (cnst.CLOVER_SAVES + 'StormLoc_-50_5000km_WA.p', 'rb'))
pdf = pdf.where((pdf.clat>=4.5) & (pdf.clat<=8.5) & (pdf.clon>=-12) & (pdf.clon<=12))
pdf = pdf.dropna()

tt = np.array(pdf.tmin, dtype=float)
tmean = tt = np.array(pdf.tmean, dtype=float)
month = pdf.month.values
area = pdf.area.values
area70 = np.array(pdf['70area'], dtype=int)
lat = pdf.lat.values
lon = lat = pdf.lon.values

era = xr.open_dataset(cnst.ERA_DAILY_PL12UTC)

era_out = get_ERA(era, pdf)

dic.keys()

cc=11

sh = np.array(dic['u925']) #np.array(dic['shear']) * (-1)
umin = np.array(dic['u650']) * (-1)
umax = np.array(dic['u925'])
qq = np.array(dic['q925']) * 1000
tt = np.array(dic['tmin'])
month = np.array(dic['month'])
area = np.array(dic['area'])*(4.4**2)
lat = np.array(dic['clat'])
lon = np.array(dic['clon'])

pos = np.where( (pp >= 0.1) & (month>=9) & (umax>0) & (umin>0))# np.where((pp >= 3) & (sh >= 8) &  (sh <= 30) &  (area<=700000) & ((month<=5) | (month>=10)) & (lat>=5) )   # 5 + 10 look nicest

tt = tt[pos]
pp = tt
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


nbs= 7
nbq= 7
nba = 7

shearbins = np.percentile(sh[(sh>=5) & (sh<=10)], np.linspace(0,100,nbs))

#np.percentile(p[(p>=7) & (p<=29)], np.arange(0,101,10)) #np.percentile(p[p>=8], np.arange(0,101,10)) #np.linspace(p[p>=8].min(), p[p>=8].max(),nbshear)
qbins = np.linspace(14,18.5, nbq) #np.linspace(16,19, nbq)#np.linspace(16,19, nbq) #np.percentile(qq[(qq>=16) & (qq<=19)], np.linspace(0,100,nbq))

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

    prob = np.sum(pp[posi]>=20) / np.sum(pp[posi]>=1)
    pprob.append(prob)
    plist.append(cmean)
    shearlist.append(((shearbins[id+1])-c)/2)


xtick = shearbins[0:-1]
xtickwidth= (shearbins[1::]-shearbins[0:-1])

fig = plt.figure(figsize=(25, 45), dpi=70)

ax1 = fig.add_subplot(221)

xy = np.vstack([sh, pp])
z = gaussian_kde(xy)(xy)
test = z / (z.max() - z.min())

r = u_stat.pcor(sh,pp, qq)
r=stats.pearsonr(sh,pp)

mappable = ax1.scatter(sh, pp, c=test, edgecolor='', cmap='viridis_r', s=20) # viridis_r
ax1.set_ylabel('Max. rainfall (mm h$^{-1}$)')
ax1.set_xlabel('u925hPa')
ax1.set_title('P-corr. u925hPa/rain | q removed: '+str(np.round(r[0], decimals=2)), fontsize=cc)
ax1.tick_params(direction='in')
cbar = fig.colorbar(mappable)
cbar.set_label('Density')

print('Partial correlation umin_mid:', u_stat.pcor(umin,pp, qq) )
####################################################################################
ax2 = fig.add_subplot(222)

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



####################################################################################
ax4 = fig.add_subplot(223)

ax4.bar(xtick, pprob, xtickwidth, align='edge', ec='black')
ax4.set_xlabel('Max. u925hPa bins (equally populated)')
ax4.set_ylabel('Probability Rainfall > 60mm h-1')
ax4.set_title('')


#####################################################################################

pp = np.array(dic['pmax'])
sh = np.array(dic['u925'])/np.abs(np.array(dic['u650'])-np.array(dic['u925'])) * (100)
qq = np.array(dic['q925']) * 1000
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


for isq, qql in enumerate(qbins[0:-1]):


    for issh, shl in enumerate(shearbins[0:-1]):

        poss = np.where((sh >= shl) & (sh < shearbins[issh+1]) & (qq>=qql) & (qq < qbins[isq+1]))

        print('bigger than',shl )
        print('smaller than', shearbins[issh+1])

        try:
            cmean = np.percentile(pp[poss], 90)

        except IndexError:
            cmean = np.nan

        prob = np.sum(pp[poss]>=30) / np.sum(pp[poss]>=1)


        outprob[issh,isq] = prob
        outperc[issh,isq] = cmean
        outval[issh,isq] = len(poss[0])

X, Y = np.meshgrid(shearbins,qbins)
cmapp = uplot.discrete_cmap(10, base_cmap='RdBu')
ax5 = fig.add_subplot(224)
#outperc[outval<30] = np.nan
Zm = ma.masked_where(np.isnan(outperc),outperc)
mappable = ax5.pcolormesh(X, Y, Zm.T, cmap=cmapp) # viridis_rvmin=20, vmax=40

ax5.set_ylabel('Max. q925hPa')
ax5.set_xlabel('Max. u925hPa (equally populated)')
ax5.set_title('')
cbar = fig.colorbar(mappable) # ticks=np.linspace(30,45,11) , ticks=np.linspace(20,40,11)
cbar.set_label('90th centile max. rain')
plt.show()
