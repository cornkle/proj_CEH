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
from metpy import calc
from metpy.units import units


def calc_trend(data, month, hour=None, method=None, sig=False, wilks=False):

    y0 = 1980
    if method is None:
        'Please provide trend calc method: polyfit or mk (mann kendall)'
    if hour is not None:

        if len(month)>1:

            data = data[((data['time.month'] >= month[0]) & (data['time.month'] <= month[1])) & (data['time.hour'] == hour) & (data['time.year'] >= y0) & (data['time.year'] <= 2018)]
        else:

            data = data[(data['time.month'] == month[0]) & (data['time.hour'] == hour) & (data['time.year'] >= y0) & (data['time.year'] <= 2018)]
    else:
        if len(month)>1:
            data = data[((data['time.month'] >= month[0]) & (data['time.month'] <= month[1]))& (data['time.year'] >= y0) & (data['time.year'] <= 2018)]
        else:
            data = data[(data['time.month'] == month[0]) & (data['time.year'] >= y0) & (data['time.year'] <= 2018)]

    if len(data.time)==0:
        print('Data does not seem to have picked month or hour. Please check input data')


    mean_years = data.groupby('time.year').mean(dim='time').squeeze().load()

    # stack lat and lon into a single dimension called allpoints
    datastacked = mean_years.stack(allpoints=['lat', 'lon'])

    # apply the function over allpoints to calculate the trend at each point
    print('Entering trend calc')

    alpha = 0.05
    # NaNs means there is not enough data, slope = 0 means there is no significant trend.
    if method=='mk':
        dtrend = datastacked.groupby('allpoints').apply(u_darrays.linear_trend_mk, alpha=alpha, eps=0.0001,nb_missing=10)
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

    # unstack back to lat lon coordinates
    return ddtrend, mean_years




def trend_all():

    pl = '/media/ck/Elements/Africa/WestAfrica/MERRA2/MERRA2_*instM_3d_*.nc' #cnst.ERA_MONTHLY_PL_SYNOP
    tl = '/media/ck/Elements/Africa/WestAfrica/MERRA2/MERRA2_*instM_2d_*.nc'
    fpath = cnst.network_data + 'figs/CLOVER/months/ERA5_WA/'

    box=[-18,30,0,25]#[-18,30,0,25]#  [-18,40,0,25] #

    da = xr.open_mfdataset(pl) #xr.open_dataset(pl)
    da2 = xr.open_mfdataset(tl)
    # da = u_darrays.flip_lat(da)
    da = da.sel(lon=slice(box[0], box[1]), lat=slice(box[2],box[3]))
    da2 = da2.sel(lon=slice(box[0], box[1]), lat=slice(box[2], box[3]))
    low_press = 950
    up_press = 650
    mid_press = 700

    q = da['QV'].sel(lev=925)
    theta_e = da['RH'].sel(lev=600)*100#da2['cape']
    tt = da['T'].sel(lev=925)
    tcw = da2['TQV']
    q.values = q.values*1000

    months= [3,10]#[3,4,5,6,9,10,11]#,4,5,6,9,10,11#,4,5,6,9,10,11,(3,5), (9,11)]#, 10,5,9]#[(12,2)]#[1,2,3,4,5,6,7,8,9,10,11,12]# #,2,3,11,12]#[(12,2)]#[1,2,3,4,5,6,7,8,9,10,11,12]# #,2,3,11,12]
    mstr = ['Mar: ', 'Oct: ']
    fp = fpath + 'use/MERRA_-70_use_nosig_RHQ.png'
    map = q.salem.get_map(countries=False)
    # Change the country borders
    map.set_shapefile(countries=True, color='grey', linewidths=0.5)

    map.set_lonlat_contours(add_ytick_labels=True, interval=5, linewidths=0.01,
                            linestyles='-', colors='white')

    pllist = []

    for ids, m in enumerate(months):
        method = 'mk'

        if type(m)==int:
            m = [m]

        mmstr = mstr[ids]

        sig = False

        qtrend, qmean = calc_trend(q, m, method=method, sig=sig, wilks=False) #hour=12,
        q_mean = qmean.mean(axis=0)

        thetatrend, thetamean = calc_trend(theta_e, m, method=method, sig=sig, wilks=False) #hour=12,
        theta_mean = thetamean.mean(axis=0)
        tttrend, ttmean = calc_trend(tt, m, method=method, sig=sig, wilks=False) #hour=12,
        tt_mean = ttmean.mean(axis=0)
        tcwtrend, tcwmean = calc_trend(tcw, m, method=method, sig=sig, wilks=False) #hour=12,
        tcw_mean = tcwmean.mean(axis=0)

        q_da = qtrend * 10.  # warming over decade
        theta_da = thetatrend * 10
        tt_da = tttrend * 10
        tcw_da = tcwtrend * 10

        pllist.append([q_da,theta_da, tt_da,tcw_da])


    f = plt.figure(figsize=(12,7), dpi=300)
    ax1 = f.add_subplot(2, 2, 1)
    map.set_shapefile(countries=True, linewidths=1.2, color='grey')
    map.set_data(pllist[0][3].values,interp='linear')  # interp='linear'
    map.set_contour((pllist[0][2].values).astype(np.float64),interp='linear', colors='k',  linewidths=1.8, levels=[-0.5,-0.4,-0.3,0.3,0.4,0.5]) #[6,8,10,12,14,16] #levels=[-0.6,-0.4,-0.2,0.2,0.4, 0.6],
    map.set_plot_params(levels=[-1.5,-1,-0.8,-0.6,-0.4,-0.2,0.2, 0.4,0.6,0.8,1,1.5], cmap='RdBu', extend='both')  # levels=np.arange(-0.5,0.51,0.1), [-0.6,-0.4,-0.2,0.2,0.4,0.6]
    dic = map.visualize(ax=ax1, title=r'Mar: TCWV (shading) | 925hPa T (contours) ', cbar_title=r'kg m$^{-2}$ decade$^{-1}$')
    contours = dic['contour'][0]
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.1f')

    ax2 = f.add_subplot(2, 2, 2)
    map.set_data(pllist[0][1].values,interp='linear')  # interp='linear'
    map.set_contour((pllist[0][0].values).astype(np.float64),interp='linear', colors='k',  linewidths=1.8, levels=[-0.6,-0.4,-0.2,0.2,0.4, 0.6]) #[6,8,10,12,14,16] #levels=[-0.6,-0.4,-0.2,0.2,0.4, 0.6],
    map.set_plot_params(levels=np.array([-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3, 0.4])*10, cmap='RdBu', extend='both')  # levels=np.arange(-0.5,0.51,0.1), [-0.6,-0.4,-0.2,0.2,0.4,0.6]
    dic = map.visualize(ax=ax2, title=r'Mar: 650hPa RH (shading) | 925hPa q (contours)', cbar_title=r'% decade$^{-1}$')
    contours = dic['contour'][0]
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.1f')

    ax3 = f.add_subplot(2,2,3)
    map.set_shapefile(countries=True, linewidths=1.2, color='grey')
    map.set_data(pllist[1][3].values,interp='linear')  # interp='linear'
    map.set_contour((pllist[1][2].values).astype(np.float64),interp='linear', colors='k',  linewidths=1.8, levels=[-0.5,-0.4,-0.3,0.3,0.4,0.5]) #[6,8,10,12,14,16] #levels=[-0.6,-0.4,-0.2,0.2,0.4, 0.6],
    map.set_plot_params(levels=[-1.5,-1,-0.8,-0.6,-0.4,-0.2,0.2, 0.4,0.6,0.8,1,1.5], cmap='RdBu', extend='both')  # levels=np.arange(-0.5,0.51,0.1), [-0.6,-0.4,-0.2,0.2,0.4,0.6]
    dic = map.visualize(ax=ax3, title=r'Oct: TCWV (shading) | 925hPa T (contours)', cbar_title=r'kg m$^{-2}$ decade$^{-1}$')
    contours = dic['contour'][0]
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.1f')

    ax4 = f.add_subplot(2, 2,4)
    map.set_data(pllist[1][1].values,interp='linear')  # interp='linear'
    map.set_contour((pllist[1][0].values).astype(np.float64),interp='linear', colors='k',  linewidths=1.8, levels=[-0.6,-0.4,-0.2,0.2,0.4, 0.6]) #[6,8,10,12,14,16] #levels=[-0.6,-0.4,-0.2,0.2,0.4, 0.6],
    map.set_plot_params(levels=np.array([-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3, 0.4])*10, cmap='RdBu', extend='both')  # levels=np.arange(-0.5,0.51,0.1), [-0.6,-0.4,-0.2,0.2,0.4,0.6]
    dic = map.visualize(ax=ax4, title=r'Oct: 650hPa RH (shading) | 925hPa q (contours)', cbar_title=r'% decade$^{-1}$')
    contours = dic['contour'][0]
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.1f')

    plt.tight_layout()

    plt.annotate('a)', xy=(0.02, 0.96), xytext=(0, 4), size=15, xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points')
    plt.annotate('b)', xy=(0.49, 0.96), xytext=(0, 4), size=15, xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points')
    plt.annotate('c)', xy=(0.02, 0.48), xytext=(0, 4), size=15, xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points')
    plt.annotate('d)', xy=(0.49, 0.48), xytext=(0, 4), size=15, xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points')

    plt.savefig(fp)
    plt.close('all')
