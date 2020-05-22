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


def calc_trend(data, month, hour=None, method=None, sig=False, wilks=False):

    if method is None:
        'Please provide trend calc method: polyfit or mk (mann kendall)'
    if hour is not None:

        if len(month)>1:

            data = data[((data['time.month'] >= month[0]) & (data['time.month'] <= month[1])) & (data['time.hour'] == hour) & (data['time.year'] >= 1983) & (data['time.year'] <= 2017)]
        else:

            data = data[(data['time.month'] == month[0]) & (data['time.hour'] == hour) & (data['time.year'] >= 1983) & (data['time.year'] <= 2017)]
    else:
        if len(month)>1:
            data = data[((data['time.month'] >= month[0]) & (data['time.month'] <= month[1]))& (data['time.year'] >= 1983) & (data['time.year'] <= 2017)]
        else:
            data = data[(data['time.month'] == month[0]) & (data['time.year'] >= 1983) & (data['time.year'] <= 2017)]

    if len(data.time)==0:
        print('Data does not seem to have picked month or hour. Please check input data')

    mean_years = data.groupby('time.year').mean('time')

    # stack lat and lon into a single dimension called allpoints
    datastacked = mean_years.stack(allpoints=['latitude', 'longitude'])


    # apply the function over allpoints to calculate the trend at each point
    print('Entering trend calc')

    alpha = 0.05
    # NaNs means there is not enough data, slope = 0 means there is no significant trend.
    if method=='mk':
        dtrend = datastacked.groupby('allpoints').apply(u_darrays.linear_trend_mk, alpha=alpha, eps=0.01,nb_missing=10)
        dtrend = dtrend.unstack('allpoints')
        if sig:
            (dtrend['slope'].values)[dtrend['ind'].values==0] = 0

    # NaNs means there is not enough data, slope = 0 means there is no significant trend.
    if method=='polyfit':
        dtrend = datastacked.groupby('allpoints').apply(u_darrays.linear_trend_lingress,nb_missing=10)
        dtrend = dtrend.unstack('allpoints')

        if sig:
            (dtrend['slope'].values)[dtrend['pval'].values > alpha] = 0

    ddtrend = dtrend['slope']

    if wilks and sig:
        try:
            pthresh = us.fdr_threshold(dtrend['pval'].values[np.isfinite(dtrend['pval'].values)], alpha=alpha)
            ddtrend.values[(dtrend['pval'].values > pthresh) | np.isnan(dtrend['pval'].values)] = np.nan
        except ValueError:
            ddtrend.values = ddtrend.values * np.nan
            pthresh = np.nan
        print('p value threshold', pthresh)
    #ipdb.set_trace()
    # unstack back to lat lon coordinates
    return ddtrend, mean_years



def trend_all():

    #mcs = cnst.GRIDSAT_PERU + 'aggs/gridsat_WA_-40_allClouds_monthly.nc'
    mcs = cnst.GRIDSAT_PERU + 'aggs/gridsat_WA_count_-55_allClouds_monthly.nc'
    fpath = cnst.network_data + 'figs/HUARAZ/'

    #box = [-79, -74, -12, -7]  # small
    #box=[-79,-65,-17,-3]#  [-18,40,0,25] #
    #box = [-80, -53, -30, -1]
    box = [-78, -75, -10.5, -8]

    da3 = xr.open_dataarray(mcs)#/100
    da3 = da3.sel(lon=slice(box[0], box[1]), lat=slice(box[2],box[3]))

    grid = da3.salem.grid.regrid(factor=1)

    tir = grid.lookup_transform(da3, method=np.nanmean)  #t2d.salem.lookup_transform(da3['tir']) #

    grid = grid.to_dataset()
    tir = xr.DataArray(tir, coords=[da3['time'],  grid['y'], grid['x']], dims=['time',  'latitude','longitude'])

    months= [11]#[3,4,5,6,9,10,11]#,4,5,6,9,10,11#,4,5,6,9,10,11,(3,5), (9,11)]#, 10,5,9]#[(12,2)]#[1,2,3,4,5,6,7,8,9,10,11,12]# #,2,3,11,12]#[(12,2)]#[1,2,3,4,5,6,7,8,9,10,11,12]# #,2,3,11,12]

    dicm = {}
    dicmean = {}


    f = plt.figure(figsize=(9,5), dpi=300)

    fname = '/home/ck/DIR/cornkle/data/HUARAZ/shapes/riosan_sel_one.shp'

    sdf = salem.read_shapefile(fname)
    #sdf = salem.transform_geopandas(sdf, to_crs=salem.wgs84)

    for ids, m in enumerate(months):
        method = 'mk'

        if type(m)==int:
            m = [m]

        sig = True


        tirtrend, tirmean = calc_trend(tir, m, method=method, sig=sig, wilks=False)

        tirm_mean = tirmean.mean('year')


        

        # plt.plot(np.arange(1985,2018), tirmean.mean(['latitude', 'longitude']))
        # return

        tirtrend_unstacked = ((tirtrend.values)*10./ tirm_mean.values) * 100.
        linemean = tirmean.where(tirtrend_unstacked > 20).mean(['latitude', 'longitude'])
        #ipdb.set_trace()
        tirtrend_out = xr.DataArray(tirtrend_unstacked, coords=[grid['y'], grid['x']], dims=['latitude','longitude'])
        tirtrend_out.name = 'tir'
        #tirmean_out = xr.DataArray(tirm_mean, coords=[grid['y'], grid['x']], dims=['latitude','longitude'])

        dicm[m[0]] = tirtrend_out
        dicmean[m[0]] = tirm_mean

        ti_da = tirtrend_out


        if len(m) == 1:
            fp = fpath + 'MCS_only_trendmap_Allmonths_count-50C_lines_small'+str(m[0]).zfill(2)+'.png'
        else:
            fp = fpath + 'MCS_only_trendmap_' + str(m[0]).zfill(2) +'-'+ str(m[1]).zfill(2) + '.png'


        ax1 = f.add_subplot(1,1,1)

        ax1.plot(linemean['year'], linemean, label=str(m))


    plt.legend()
    plt.tight_layout()
    plt.savefig(fp)
    plt.close('all')
