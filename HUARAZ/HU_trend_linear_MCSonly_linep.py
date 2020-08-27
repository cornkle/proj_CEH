import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from utils import u_darrays
import ipdb
from utils import constants as cnst, u_met
import salem
from utils import u_statistics as us
from scipy import stats
import numpy.ma as ma
import pickle as pkl
import shapely.geometry as shpg
import seaborn
import pandas as pd


def trend_all():

    #mcs = cnst.GRIDSAT_PERU + 'aggs/gridsat_WA_-40_allClouds_monthly.nc'
    mcs = cnst.GRIDSAT_PERU + 'aggs/gridsat_WA_count_-50_allClouds_monthly.nc'
    chirps = '/media/ck/Elements/SouthAmerica/CHIRPS/chirps-v2.0.monthly.nc'
    enso = '/home/ck/DIR/mymachine/ENSO/ONI.csv'#'/home/ck/DIR/mymachine/ENSO/meiv2.data'
    fpath = cnst.network_data + 'figs/HUARAZ/'

    fname = '/home/ck/DIR/cornkle/data/HUARAZ/shapes/riosan_sel_one.shp'
    isbuffer = [-79, -74, -12, -7]

    sdf = salem.read_shapefile(fname)
    sdf = salem.transform_geopandas(sdf, to_crs=salem.wgs84)

    da3 = xr.open_dataarray(mcs).sel(lon=slice(isbuffer[0], isbuffer[1]), lat=slice(isbuffer[2], isbuffer[3]))
    ca = xr.open_dataarray(chirps).sel(longitude=slice(isbuffer[0], isbuffer[1]), latitude=slice(isbuffer[2], isbuffer[3]))
    # This masks out the data which is not in the region

    ens = pd.read_csv(enso, sep=',', engine='python', names=np.arange(0, 13),index_col=0)

    ca[0,:,:].salem.roi(shape=sdf).plot.pcolormesh()

    da3 = da3.salem.roi(shape=sdf).mean(['lat', 'lon'])*100
    ca = ca.salem.roi(shape=sdf).mean(['latitude', 'longitude'])
    months= [1,2,3,4,5,6,7,8,9,10,11,12]#,4,5,6,9,10,11#,4,5,6,9,10,11,(3,5), (9,11)]#, 10,5,9]#[(12,2)]#[1,2,3,4,5,6,7,8,9,10,11,12]# #,2,3,11,12]#[(12,2)]#[1,2,3,4,5,6,7,8,9,10,11,12]# #,2,3,11,12]


    dicm = {}
    dicmean = {}




    f = plt.figure(figsize=(15,7.5), dpi=300)


    #sdf = salem.transform_geopandas(sdf, to_crs=salem.wgs84)

    for ids, m in enumerate(months):
        method = 'mk'

        if type(m)==int:
            m = [m]

        ensmonth = ens[m[0]]

        eens = ensmonth.loc['1985':'2019']

        sig = True

        da = da3[(da3['time.month'] == m[0]) & (da3['time.year'] >= 1985) & (
                da3['time.year'] <= 2019)]

        ch = ca[(ca['time.month'] == m[0]) & (ca['time.year'] >= 1985) & (
                ca['time.year'] <= 2019)]

        da = da - np.mean(da.values)
        ch = ch - np.mean(ch.values)

        sslope, sint, srval, spval, serr = stats.linregress(np.arange(len(da.values)), da.values)
        print('linear regression for shear', m, sslope, srval, spval)

        cslope, cint, crval, cpval, cerr = stats.linregress(np.arange(len(ch.values)), ch.values)
        print('linear regression for shear', m, cslope, crval, cpval)

        #ipdb.set_trace()


        dslope = sslope*10
        cdslope = cslope*10

        if len(m) == 1:
            fp = fpath + 'MCS_only_trendmap_Allmonths_count-50C_lines_Huaraz_1985-2018_'+str(m[0]).zfill(2)+'.png'
        else:
            fp = fpath + 'MCS_only_trendmap_1985-2018_' + str(m[0]).zfill(2) +'-'+ str(m[1]).zfill(2) + '.png'


        ax = f.add_subplot(3,4,ids+1)

        x = np.arange(0, len(ch['time.year']))

        ax.plot(ch['time.year'], ch, marker='o', markersize=3, label='Trend: '+str(np.round(cdslope,2))+ 'mm / month decade | '+ 'p='+str(np.round(cpval,2)), color='blue')
        ax.plot(ch['time.year'], cint + cslope * np.arange(0, len(ch['time.year'])), linestyle='dashed', color='blue')
        cc = []
        for enb in eens.values:
            if enb < 0:
                cc.append('lightblue')
            else:
                cc. append('red')

        ax.bar(eens.index.values, eens.values*30, color=cc)
        plt.title('Month: '+ str(m[0]))
        ax.set_ylabel('mm month$^{-1}$ (blue)')
        ax.set_ylim(-90,90)
        plt.legend(fontsize=6)
        ax1 = ax.twinx()
        ax1.plot(da['time.year'], da, color='orange', marker='o', markersize=3)
        ax1.plot(da['time.year'], sint + sslope * np.arange(0, len(da['time.year'])), color='orange', linestyle='dashed')
        ax1.set_ylabel('% cloud cover/month (orange)')
        ax1.set_ylim(-20, 20)


    plt.tight_layout()
    plt.savefig(fp)
    plt.close('all')
