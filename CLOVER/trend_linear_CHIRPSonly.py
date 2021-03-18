import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from utils import u_darrays
import ipdb
from utils import constants as cnst, u_met
import salem
from utils import u_statistics as us, u_darrays as uda
from scipy import stats
import numpy.ma as ma
import pickle as pkl
import shapely.geometry as shpg
import seaborn


def calc_trend(data, month, method=None, sig=False, wilks=False):

    y1 = 1983
    y2 = 2019

    if method is None:
        'Please provide trend calc method: polyfit or mk (mann kendall)'

    data = data[(data['year'] >= y1) & (data['year'] <= y2)]


    mean_years = data.mean('year')

    # stack lat and lon into a single dimension called allpoints
    datastacked = mean_years.stack(allpoints=['latitude', 'longitude'])


    # apply the function over allpoints to calculate the trend at each point
    print('Entering trend calc')

    alpha = 0.05
    # NaNs means there is not enough data, slope = 0 means there is no significant trend.
    if method=='mk':
        dtrend = datastacked.groupby('allpoints').apply(u_darrays.linear_trend_mk, alpha=alpha, eps=0.01,nb_valid=10)
        dtrend = dtrend.unstack('allpoints')
        if sig:
            (dtrend['slope'].values)[dtrend['ind'].values==0] = 0

    # NaNs means there is not enough data, slope = 0 means there is no significant trend.
    if method=='polyfit':
        dtrend = datastacked.groupby('allpoints').apply(u_darrays.linear_trend_lingress,nb_valid=10)
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
    #ipdb.set_trace()
    ddout = ddtrend #uda.flip_lat(ddtrend)

    return ddout, mean_years



def trend_all():

    #mcs = cnst.GRIDSAT_PERU + 'aggs/gridsat_WA_-40_allClouds_monthly.nc'
    #mcs = cnst.GRIDSAT_PERU + 'aggs/gridsat_WA_count_-50_allClouds_monthly.nc'

    mcs = '/media/ck/Elements/Africa/CHIRPS/daily/*.nc' #'/media/ck/Elements/SouthAmerica/CHIRPS/chirps-v2.0.monthly.nc'

    fpath = cnst.network_data + 'figs/CLOVER/'

    box =  [-14,14,3.5,15] #[-78, -77, -10.2, -8.63]#[-80, -65, -20, -1] #[-79, -76, -11, -8] # [-79, -74, -12, -8]  # small
    #box=[-79,-65,-17,-3]#  [-18,40,0,25] #
    #box = [-80, -53, -30, -1]
    #box = [-79,-69,-17,-7] #[-79, -75, -10.5, -8]

    da3_out = xr.open_mfdataset(mcs)#/100
    da3_out = da3_out['precip'].sel(longitude=slice(box[0], box[1]), latitude=slice(box[2],box[3]))#.sel(latitude=slice(-10.2,-8.63), longitude=slice(-78, -77))#.sel(longitude=slice(box[0], box[1]), latitude=slice(box[2],box[3]))

    #grid = da3.salem.grid.regrid(factor=1)

    #tir = grid.lookup_transform(da3, method=np.nanmean)  #t2d.salem.lookup_transform(da3['tir']) #

    #grid = grid.to_dataset()

    months= [1,2,3,4,5,6,7,8,9,10,11,12]#[3,4,5,6,9,10,11]#,4,5,6,9,10,11#,4,5,6,9,10,11,(3,5), (9,11)]#, 10,5,9]#[(12,2)]#[1,2,3,4,5,6,7,8,9,10,11,12]# #,2,3,11,12]#[(12,2)]#[1,2,3,4,5,6,7,8,9,10,11,12]# #,2,3,11,12]

    dicm = {}
    dicmean = {}

    md = {1:31, 2:28, 3:31, 4:30, 5:31, 6:30,7:31, 8:31, 9:30, 10:31, 11:30, 12:31}


    f = plt.figure(figsize=(14, 8), dpi=300)


    for ids, m in enumerate(months):

        tir = da3_out[da3_out['time.month']==m].load()

        tir.values[tir.values<5] = np.nan

        tir = tir.groupby('time.year').count('time')/md[m]

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

        tirtrend_unstacked = ((tirtrend.values)*10.)#/ tirm_mean.values) * 100. #/ tirm_mean.values
        #ipdb.set_trace()
        tirtrend_out = xr.DataArray(tirtrend_unstacked, coords=[tir.latitude, tir.longitude], dims=['latitude','longitude'])
        tirtrend_out.name = 'tir'
        tirmean_out = xr.DataArray(tirm_mean, coords=[tir.latitude, tir.longitude], dims=['latitude','longitude'])

        # dicm[m[0]] = tirtrend_out
        # dicmean[m[0]] = tirm_mean
        #
        ti_da = tirtrend_out


        if len(m) == 1:
            fp = fpath + 'CHIRPS_only_trendmap_Allmonths_SIG_inPerc_mediumDomain_1983-2018_'+str(m[0]).zfill(2)+'.png'
        else:
            fp = fpath + 'CHIRPS_only_trendmap_' + str(m[0]).zfill(2) +'-'+ str(m[1]).zfill(2) + '.png'

        map = tir.salem.get_map()

        ax1 = f.add_subplot(3,4,ids+1)
        map.set_contour((tirm_mean), interp='linear', levels=[10, 20, 40,60,80,100,200,300], colors='k', linewidths=0.5)

        ti_da.values[ti_da.values==0] = np.nan

        map.set_data(ti_da)

        map.set_plot_params(cmap='RdBu_r', extend='both', levels = np.linspace(-20, 20, 8)) #)  #, levels=np.arange(-7,7,25)  # levels=np.arange(20,101,20)  #np.arange(20,101,20)
        #map.set_plot_params(cmap='RdBu', extend='both', levels=[-100, -80,-60,-40,-20,-10,-5,-2,0,2,5,10,20,40,60,80,100])
        dic = map.visualize(ax=ax1, title=str(m)+': Monthly precipitation change', cbar_title='% decade-1',addcbar=True)
        contours = dic['contour'][0]
        plt.clabel(contours, inline=True, fontsize=7, fmt='%1.1f')

    plt.tight_layout()
    plt.savefig(fp)
    plt.close('all')
