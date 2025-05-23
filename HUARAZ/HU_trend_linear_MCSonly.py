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

    y1 = 1984
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
    mcs = cnst.GRIDSAT_PERU + 'aggs/gridsat_WA_count_-40_allClouds_monthly.nc'

    fpath = cnst.network_data + 'figs/HUARAZ/'

    box = [-81, -33, -55, 11] #[-79, -74, -12, -8]  # small latitude=slice(-25,0), longitude=slice(-81,-65) [-81,-65,-25,0]
    #box=[-79,-65,-17,-3]#  [-18,40,0,25] #
    #box = [-80, -53, -30, -1]
    #box = [-79,-69,-17,-7] #[-79, -75, -10.5, -8]

    da3 = xr.open_dataarray(mcs)#/100
    da3 = da3.sel(lon=slice(box[0], box[1]), lat=slice(box[2],box[3]))
    print('Doing regrid')
    grid = da3.salem.grid.regrid(factor=0.25)

    tir = grid.lookup_transform(da3, method=np.nanmean)  #t2d.salem.lookup_transform(da3['tir']) #

    grid = grid.to_dataset()
    tir = xr.DataArray(tir, coords=[da3['time'],  grid['y'], grid['x']], dims=['time',  'latitude','longitude'])
    print('Created array')
    tir.to_netcdf(cnst.GRIDSAT_PERU + 'aggs/gridsat_WA_count_-40_allClouds_monthly_regrid025.nc')
    ipdb.set_trace()
    months= [1,2,3,4,5,6,7,8,9,10,11,12]#[3,4,5,6,9,10,11]#,4,5,6,9,10,11#,4,5,6,9,10,11,(3,5), (9,11)]#, 10,5,9]#[(12,2)]#[1,2,3,4,5,6,7,8,9,10,11,12]# #,2,3,11,12]#[(12,2)]#[1,2,3,4,5,6,7,8,9,10,11,12]# #,2,3,11,12]

    dicm = {}
    dicmean = {}


    f = plt.figure(figsize=(14, 9), dpi=300)

    fname = '/home/ck/DIR/cornkle/data/HUARAZ/shapes/riosan_sel_one.shp'

    sdf = salem.read_shapefile(fname)
    #sdf = salem.transform_geopandas(sdf, to_crs=salem.wgs84)

    for ids, m in enumerate(months):

        print('Doing', m)
        method = 'mk'

        if type(m)==int:
            m = [m]

        sig = True


        tirtrend, tirmean = calc_trend(tir, m, method=method, sig=sig, wilks=False)

        tirm_mean = tirmean.mean('year')

        #linemean = tirmean.mean(['year','lat','lon'])


        # plt.plot(np.arange(1985,2018), tirmean.mean(['latitude', 'longitude']))
        # return

        tirtrend_unstacked = ((tirtrend.values)*10./ tirm_mean.values) * 100. #/ tirm_mean.values
        #ipdb.set_trace()
        tirtrend_out = xr.DataArray(tirtrend_unstacked, coords=[grid['y'], grid['x']], dims=['latitude','longitude'])
        tirtrend_out.name = 'tir'
        tirmean_out = xr.DataArray(tirm_mean, coords=[grid['y'], grid['x']], dims=['latitude','longitude'])

        # dicm[m[0]] = tirtrend_out
        # dicmean[m[0]] = tirm_mean
        #
        ti_da = tirtrend_out


        if len(m) == 1:
            fp = fpath + 'all_coldCloud_-40_1985-2019_'+str(m[0]).zfill(2)+'.png'
        else:
            fp = fpath + 'all_coldCloud_-40_1985-2019_' + str(m[0]).zfill(2) +'-'+ str(m[1]).zfill(2) + '.png'

        map = ti_da.salem.get_map()

        ax1 = f.add_subplot(3,4,ids+1)
        map.set_contour((tirm_mean)*100, interp='linear', levels=[10, 20, 40,60,80], colors='k', linewidths=0.5)
        #map.set_contour((tirm_mean), interp='linear', levels=[-55,-50,-45], colors='k', linewidths=0.5)
        #.values).astype(np.float64)

        ti_da.values[ti_da.values==0] = np.nan
        map.set_data(ti_da)  #
        #map.set_data(tirm_mean)
        map.set_shapefile(sdf,facecolor='white', color='red', linewidth=1.5, alpha=0.5)
        #map.set_geometry(sdf, linewidth=2)
        # coord = [18, 25, -28, -20]
        # geom = shpg.box(coord[0], coord[2], coord[1], coord[3])
        #map.set_geometry(geom, zorder=99, color='darkorange', linewidth=3, linestyle='--', alpha=0.3)

        #map.set_plot_params(cmap='RdBu_r', extend='both', levels=np.arange(-1.5,1.6,0.25)) #)  #, levels=np.arange(-7,7,25)  # levels=np.arange(20,101,20)  #np.arange(20,101,20)
        map.set_plot_params(cmap='RdBu', extend='both', levels=np.arange(-30, 31, 10))
        dic = map.visualize(ax=ax1, title=str(m)+': -40C frequency change', cbar_title='% decade-1',addcbar=True)
        contours = dic['contour'][0]
        plt.clabel(contours, inline=True, fontsize=7, fmt='%1.1f')

    plt.tight_layout()
    plt.savefig(fp)
    plt.close('all')
