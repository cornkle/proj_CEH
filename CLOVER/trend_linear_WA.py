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


    mean_years = data.groupby('time.year').mean(axis=0)

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

    # unstack back to lat lon coordinates
    return ddtrend, mean_years




def trend_all():

    srfc = cnst.ERA5_MONTHLY_SRFC_SYNOP #cnst.ERA_MONTHLY_SRFC_SYNOP
    pl = cnst.ERA5_MONTHLY_PL_SYNOP #cnst.ERA_MONTHLY_PL_SYNOP
    mcs = cnst.GRIDSAT + 'aggs/gridsat_WA_-70_monthly_mean_5000km2.nc'

    fpath = cnst.network_data + 'figs/CLOVER/months/ERA5_WA/'

    box=[-18,30,0,25]#  [-18,40,0,25] #

    da = xr.open_dataset(pl) #xr.open_dataset(pl)
    #da = xr.decode_cf(da)
    da = u_darrays.flip_lat(da)
    da = da.sel(longitude=slice(box[0], box[1]), latitude=slice(box[2],box[3]))
    da2 = xr.open_dataset(srfc) #xr.open_dataset(srfc)
    #da2 = xr.decode_cf(da2)
    da2 = u_darrays.flip_lat(da2)
    da2 = da2.sel(longitude=slice(box[0], box[1]), latitude=slice(box[2],box[3]))
    da3 = xr.open_dataarray(mcs)*100
    da3 = da3.sel(lon=slice(box[0], box[1]), lat=slice(box[2],box[3]))

    lons = da.longitude
    lats = da.latitude


    press = da2['tcwv']
    press = press[press['time.hour'] == 12]
    #press.values = press.values#*1000
    low_press = 925
    up_press = 650
    mid_press = 700

    q = da['q'].sel(level=slice(low_press-20, low_press)).mean('level')
    q = q[q['time.hour']==12]
    t2d = da2['t2m']#['t2m']
    #t2d = da['t'].sel(level=slice(800, 850)).mean('level')
    t2d = t2d[t2d['time.hour']==12]

    theta_low = u_met.theta_e(da.level.values, da['t'].sel(level=low_press), da['q'].sel(level=low_press))
    theta_high = u_met.theta_e(da.level.values, da['t'].sel(level=mid_press), da['q'].sel(level=mid_press))

    theta_e = theta_low - theta_high

    u600 = da['u'].sel(level=slice(up_press-20, up_press)).mean('level')
    u600 = u600[u600['time.hour']==12]
    v600 = da['v'].sel(level=slice(up_press-20, up_press)).mean('level')
    v600 = v600[v600['time.hour']==12]

    ws600 = u_met.u_v_to_ws_wd(u600, v600)

    u800 = da['u'].sel(level=slice(low_press-20, low_press)).mean('level')
    u800 = u800[u800['time.hour']==12]

    v800 = da['v'].sel(level=slice(low_press-20, low_press)).mean('level')
    v800 = v800[v800['time.hour']==12]

    shear_u = u600-u800 #u600-
    shear_v = v600-v800 # v600-

    ws_shear = u_met.u_v_to_ws_wd(shear_u.values, shear_v.values)

    ws_600 = t2d.copy(deep=True)
    ws_600.name = 'ws'

    ws_600.values = ws600[0]

    shear = t2d.copy(deep=True)
    shear.name = 'shear'
    shear.values = ws_shear[0]

    u6 = shear_u#u800
    v6 = shear_v#v800

    q.values = q.values*1000

    grid = t2d.salem.grid.regrid(factor=0.5)
    t2 = t2d # grid.lookup_transform(t2d)
    tir = grid.lookup_transform(da3)  #t2d.salem.lookup_transform(da3['tir']) #

    grid = grid.to_dataset()
    tir = xr.DataArray(tir, coords=[da3['time'],  grid['y'], grid['x']], dims=['time',  'latitude','longitude'])

    months= [2,3,4,5,6,7,8,9,10,11]#[3,4,5,6,9,10,11]#,4,5,6,9,10,11#,4,5,6,9,10,11,(3,5), (9,11)]#, 10,5,9]#[(12,2)]#[1,2,3,4,5,6,7,8,9,10,11,12]# #,2,3,11,12]#[(12,2)]#[1,2,3,4,5,6,7,8,9,10,11,12]# #,2,3,11,12]

    dicm = {}
    dicmean = {}

    for m in months:
        method = 'polyfit'

        if type(m)==int:
            m = [m]

        sig = True

        t2trend, t2mean = calc_trend(t2, m,  method=method, sig=sig,hour=12, wilks=False) #hour=12,
        t2_mean = t2mean.mean(axis=0)

        tirtrend, tirmean = calc_trend(tir, m, method=method, sig=sig, wilks=False)

        tirm_mean = tirmean.mean(axis=0)

        qtrend, qmean = calc_trend(q, m, method=method, sig=sig,hour=12, wilks=False) #hour=12,
        q_mean = qmean.mean(axis=0)

        sheartrend, shearmean = calc_trend(shear, m, method=method, sig=sig,hour=12, wilks=False) #hour=12,
        shear_mean = shearmean.mean(axis=0)

        presstrend, pressmean = calc_trend(press, m, method=method, sig=sig,hour=12, wilks=False) #hour=12,
        press_mean = pressmean.mean(axis=0)

        u6trend, u6mean = calc_trend(u6, m,  method=method, sig=sig, hour=12,wilks=False) #hour=12,
        u6_mean = u6mean.mean(axis=0)
        v6trend, v6mean = calc_trend(v6, m, method=method, sig=sig, hour=12,wilks=False) #hour=12,
        v6_mean = v6mean.mean(axis=0)

        thetatrend, thetamean = calc_trend(theta_e, m, method=method, sig=sig, hour=12,wilks=False) #hour=12,
        theta_mean = thetamean.mean(axis=0)

        t2trend_unstacked = t2trend*10. # warming over decade
        qtrend_unstacked = qtrend * 10.  # warming over decade
        sheartrend_unstacked = sheartrend * 10.  # warming over decade
        u6trend_unstacked = u6trend * 10
        v6trend_unstacked = v6trend * 10
        presstrend_unstacked = presstrend * 10
        thetatrend_unstacked = thetatrend * 10

        tirtrend_unstacked = ((tirtrend.values)*10. / tirm_mean.values) * 100.
        #ipdb.set_trace()
        tirtrend_out = xr.DataArray(tirtrend_unstacked, coords=[grid['y'], grid['x']], dims=['latitude','longitude'])
        tirtrend_out.name = 'tir'
        #tirmean_out = xr.DataArray(tirm_mean, coords=[grid['y'], grid['x']], dims=['latitude','longitude'])

        dicm[m[0]] = tirtrend_out
        dicmean[m[0]] = tirm_mean

        t_da = t2trend_unstacked
        q_da = qtrend_unstacked
        s_da = sheartrend_unstacked
        ti_da = tirtrend_out
        tcwv_da = presstrend_unstacked
        theta_da  = thetatrend_unstacked

        if len(m) == 1:
            fp = fpath + 'use/ERA5_trend_synop_WA_sig_poly_tcwv_'+str(m[0]).zfill(2)+'.png'
        else:
            fp = fpath + 'use/ERA5_trend_synop_WA_sig_poly_tcwv_' + str(m[0]).zfill(2) +'-'+ str(m[1]).zfill(2) + '.png'
        map = shear.salem.get_map()
        ti_da = t2d.salem.transform(ti_da)

        f = plt.figure(figsize=(15,8), dpi=300)

        # transform their coordinates to the map reference system and plot the arrows
        xx, yy = map.grid.transform(shear.longitude.values, shear.latitude.values,
                                    crs=shear.salem.grid.proj)

        xx, yy = np.meshgrid(xx, yy)
        #Quiver only every 7th grid point
        u = u6trend_unstacked.values[1::2, 1::2]
        v = v6trend_unstacked.values[1::2, 1::2]

        #Quiver only every 7th grid point
        uu = u6_mean.values[1::2, 1::2]
        vv = v6_mean.values[1::2, 1::2]

        xx = xx[1::2, 1::2]
        yy = yy[1::2, 1::2]

        pdic = {
            'tlin' : (t2_mean.values-273.15).astype(np.float64),
            'tmean' : (t2_mean.values-273.15).astype(np.float64),
            'qmean' : (q_mean.values).astype(np.float64),
            'qlin'  : q_da.values,
            'shearlin' : s_da.values,
            'u' : u,
            'v' : v,
            'xx' : xx,
            'yy' : yy,
            'tirmean' : tirm_mean,


        }

        pkl.dump(dicm,
                 open(cnst.network_data + 'data/CLOVER/saves/storm_frac_synop12UTC_WA.p',
                      'wb'))

        ax1 = f.add_subplot(221)
        map.set_data(t_da.values, interp='linear')  # interp='linear'

        map.set_contour((t2_mean.values-273.15).astype(np.float64), interp='linear', colors='k', linewidths=0.5, levels=[20,23,26,29,32,35])
        map.set_plot_params(levels=[-0.5,-0.4,-0.3,-0.2,0.2,0.3,0.4,0.5], cmap='RdBu_r', extend='both')  # levels=np.arange(-0.5,0.51,0.1),

        #map.set_contour((t2_mean.values).astype(np.float64), interp='linear', colors='k', linewidths=0.5, levels=np.linspace(800,925,8))
        #map.set_plot_params(levels=[-0.5,-0.4,-0.3,-0.2,-0.1,-0.05,-0.02, 0.02,0.05,0.1,0.2,0.3,0.4,0.5], cmap='RdBu_r', extend='both')  # levels=np.arange(-0.5,0.51,0.1),

        dic = map.visualize(ax=ax1, title='2m temperature trend | contours: mean T', cbar_title='K decade-1')
        contours = dic['contour'][0]
        plt.clabel(contours, inline=True, fontsize=7, fmt='%1.1f')

        ax2 = f.add_subplot(222)
        map.set_data(tcwv_da.values,interp='linear')  # interp='linear'
        map.set_contour((press_mean.values).astype(np.float64),interp='linear', colors='k', levels=[20,30,40,50,60], linewidths=0.5) #[6,8,10,12,14,16]
        map.set_plot_params(levels=[-0.8,-0.6,-0.4,-0.2,0.2,0.4, 0.6,0.8], cmap='RdBu', extend='both')  # levels=np.arange(-0.5,0.51,0.1), [-0.6,-0.4,-0.2,0.2,0.4,0.6]

        dic = map.visualize(ax=ax2, title='925hPa Spec. humidity trend | contours: mean q', cbar_title='g kg-1 decade-1')
        contours = dic['contour'][0]
        plt.clabel(contours, inline=True, fontsize=7, fmt='%1.1f')


        ax3 = f.add_subplot(223)
        map.set_data(s_da.values, interp='linear')  # interp='linear'
        map.set_contour(s_da.values, interp='linear', levels=np.arange(-15,-8,8), colors='k')

        map.set_plot_params(levels=[-0.8,-0.6,-0.4,-0.2,0.2,0.4, 0.6,0.8], cmap='RdBu_r', extend='both')  # levels=np.arange(-0.5,0.51,0.1)
        map.visualize(ax=ax3, title='925-650hPa wind shear trend, mean 650hPa wind vectors', cbar_title='m s-1 decade-1')
        qu = ax3.quiver(xx, yy, u, v, scale=60, width=0.002)

        qk = plt.quiverkey(qu, 0.4, 0.03, 1, '1 m s$^{-1}$',
                           labelpos='E', coordinates='figure')

        ax4 = f.add_subplot(224)
        map.set_contour((tirm_mean), interp='linear', levels=[0.1,0.5,1,2.5], colors='k', linewidths=0.5)
        #.values).astype(np.float64)

        ti_da.values[ti_da.values==0] = np.nan
        map.set_data(ti_da)  #
        coord = [18, 25, -28, -20]
        geom = shpg.box(coord[0], coord[2], coord[1], coord[3])
        #map.set_geometry(geom, zorder=99, color='darkorange', linewidth=3, linestyle='--', alpha=0.3)

        map.set_plot_params(cmap='viridis', extend='both', levels=np.arange(10,51,10))  # levels=np.arange(20,101,20)  #np.arange(20,101,20)
        dic = map.visualize(ax=ax4, title='-70C cloud cover change | >5000km2', cbar_title='$\%$ decade-1')
        contours = dic['contour'][0]
        plt.clabel(contours, inline=True, fontsize=7, fmt='%1.1f')

        plt.tight_layout()
        plt.savefig(fp)
        plt.close('all')

    pkl.dump(dicm,
             open(cnst.network_data + 'data/CLOVER/saves/storm_frac_synop12UTC_WA.p',
                  'wb'))

    pkl.dump(dicmean,
                 open(cnst.network_data + 'data/CLOVER/saves/storm_frac_mean_synop12UTC_WA.p',
                      'wb'))
