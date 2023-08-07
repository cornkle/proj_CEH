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

    y1 = 2000
    y2 = 2019

    if method is None:
        'Please provide trend calc method: polyfit or mk (mann kendall)'
    if hour is not None:

        if len(month)>1:

            data = data[((data['time.month'] >= month[0]) & (data['time.month'] <= month[1])) & (data['time.hour'] == hour) & (data['time.year'] >= y1) & (data['time.year'] <= y2)]
        else:

            data = data[(data['time.month'] == month[0]) & (data['time.hour'] == hour) & (data['time.year'] >= y1) & (data['time.year'] <= y2)]
    else:
        if len(month)>1:
            data = data[((data['time.month'] >= month[0]) & (data['time.month'] <= month[1]))& (data['time.year'] >= y1) & (data['time.year'] <= y2)]
        else:
            data = data[(data['time.month'] == month[0]) & (data['time.year'] >= y1) & (data['time.year'] <= y2)]

    if len(data.time)==0:
        print('Data does not seem to have picked month or hour. Please check input data')

    mean_years = data.groupby('time.year').mean('time')

    # stack lat and lon into a single dimension called allpoints
    datastacked = mean_years.stack(allpoints=['lat', 'lon'])


    # apply the function over allpoints to calculate the trend at each point
    print('Entering trend calc')

    alpha = 0.05
    # NaNs means there is not enough data, slope = 0 means there is no significant trend.
    if method=='mk':
        dtrend = datastacked.groupby('allpoints').apply(u_darrays.linear_trend_mk, alpha=alpha, eps=0.01,nb_missing=5)
        dtrend = dtrend.unstack('allpoints')
        if sig:
            (dtrend['slope'].values)[dtrend['ind'].values==0] = 0

    # NaNs means there is not enough data, slope = 0 means there is no significant trend.
    if method=='polyfit':
        dtrend = datastacked.groupby('allpoints').apply(u_darrays.linear_trend_lingress,nb_missing=5)
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
    mcs = '/home/ck/DIR/cornkle/data/HUARAZ/Lorenz_NDVI/RioSanta_NDVI_monthly.nc'

    fpath = cnst.network_data + 'figs/HUARAZ/'


    da3 = xr.open_dataarray(mcs)#/100

    #ipdb.set_trace()

    months= [1,2,3,4,5,6,7,8,9,10,11,12]#[3,4,5,6,9,10,11]#,4,5,6,9,10,11#,4,5,6,9,10,11,(3,5), (9,11)]#, 10,5,9]#[(12,2)]#[1,2,3,4,5,6,7,8,9,10,11,12]# #,2,3,11,12]#[(12,2)]#[1,2,3,4,5,6,7,8,9,10,11,12]# #,2,3,11,12]


    f = plt.figure(figsize=(14, 9), dpi=300)


    for ids, m in enumerate(months[0:2]):

        print('Doing', m)
        method = 'mk'

        if type(m)==int:
            m = [m]

        sig = False

        tirtrend, tirmean = calc_trend(da3, m, method=method, sig=sig, wilks=False)

        tirm_mean = tirmean.mean('year')

        #ipdb.set_trace()
        tirtrend_out = xr.DataArray(tirtrend, coords=[da3['lat'], da3['lon']], dims=['latitude','longitude'])
        tirtrend_out.name = 'ndvi'
        tirmean_out = xr.DataArray(tirm_mean, coords=[da3['lat'], da3['lon']], dims=['latitude','longitude'])

        ti_da = tirtrend_out

        if len(m) == 1:
            fp = fpath + 'NDVI_nosig_'+str(m[0]).zfill(2)+'.png'
        else:
            fp = fpath + 'NDVI_nosig_' + str(m[0]).zfill(2) +'-'+ str(m[1]).zfill(2) + '.png'

        map = ti_da.salem.get_map()

        ax1 = f.add_subplot(3, 4, ids + 1)
        map.set_contour((tirm_mean) * 100, interp='linear', levels=[10, 20, 40, 60, 80], colors='k', linewidths=0.5)
        # map.set_contour((tirm_mean), interp='linear', levels=[-55,-50,-45], colors='k', linewidths=0.5)
        # .values).astype(np.float64)

        ti_da.values[ti_da.values == 0] = np.nan
        map.set_data(tirtrend)  #
        # map.set_data(tirm_mean)

        # map.set_geometry(sdf, linewidth=2)
        # coord = [18, 25, -28, -20]
        # geom = shpg.box(coord[0], coord[2], coord[1], coord[3])
        # map.set_geometry(geom, zorder=99, color='darkorange', linewidth=3, linestyle='--', alpha=0.3)

        # map.set_plot_params(cmap='RdBu_r', extend='both', levels=np.arange(-1.5,1.6,0.25)) #)  #, levels=np.arange(-7,7,25)  # levels=np.arange(20,101,20)  #np.arange(20,101,20)
        map.set_plot_params(cmap='RdBu', extend='both', levels=np.arange(-30, 31, 10))
        dic = map.visualize(ax=ax1, title=str(m) + ': -50C frequency change', cbar_title='% decade-1', addcbar=True)
        contours = dic['contour'][0]
        plt.clabel(contours, inline=True, fontsize=7, fmt='%1.1f')

    plt.tight_layout()
    plt.savefig(fp)
    plt.close('all')
