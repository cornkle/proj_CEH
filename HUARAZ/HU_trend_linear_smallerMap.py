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

    #ipdb.set_trace()
    mean_years = data.groupby('time.year').mean('time')

    highpos = 0


    if mean_years.name == 'z':
        highpos = []
        three = mean_years.coarsen(year=11, boundary='trim').mean()
        for enu, my in enumerate(three):

            if enu>0:
                pos = np.unravel_index(np.argmax(my.values), my.values.shape)
                my.values[pos]=0
                pos = np.unravel_index(np.argmax(my.values), my.values.shape)
                my.values[pos] = 0
                pos = np.unravel_index(np.argmax(my.values), my.values.shape)
                my.values[pos] = 0
                pos = np.unravel_index(np.argmax(my.values), my.values.shape)
                my.values[pos] = 0


            pos = np.unravel_index(np.argmax(my.values), my.values.shape)
            #ipdb.set_trace()
            mlat = my.latitude.values[pos[0]]
            mlon = my.longitude.values[pos[1]]
            highpos.append((mlon,mlat))

    if mean_years.name == 'u200':

        three = mean_years.coarsen(year=11, boundary='trim').mean()
        three.values = np.abs(three.values)
       # ipdb.set_trace()
        highpos = three.argmin('latitude')



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
    return ddtrend, mean_years, highpos


def get_trend(var, m, sig=False, wilks=False, method='polyfit'):
    vartrend, varmean, isnt = calc_trend(var, m, method=method, sig=sig, hour=12, wilks=False)  # hour=12,
    varmean = varmean.mean(axis=0)
    vartrend = vartrend * 10.

    return (vartrend, varmean, isnt)




def trend_all():


    pl = cnst.ERA5_MONTHLY_PL_SYNOP_HU #cnst.ERA_MONTHLY_PL_SYNOP
    #mcs = cnst.GRIDSAT + 'aggs/gridsat_WA_-70_monthly_mean_5000km2.nc'

    fpath = cnst.network_data + 'figs/HUARAZ/monthly/'

    box=[-82,-40,-28,4]#  [-18,40,0,25] #

    topo = xr.open_dataset('/media/ck/Elements/SouthAmerica/ERA5/monthly/ERA5_static_synop_0.7deg.nc')
    topo = u_darrays.flip_lat(topo)
    z = topo['z'].isel(number=0, time=0)
    z = z.sel(longitude=slice(box[0], box[1]), latitude=slice(box[2],box[3])).values

    da = xr.open_mfdataset(pl+'/*.nc') #xr.open_dataset(pl)
    #da = xr.decode_cf(da)
    da = u_darrays.flip_lat(da)
    da = da.sel(longitude=slice(box[0], box[1]), latitude=slice(box[2],box[3]), time=(da['time.hour']==12)).load()
    #ipdb.set_trace()
    #da3 = xr.open_dataarray(mcs)*100
    #da3 = da3.sel(lon=slice(box[0], box[1]), lat=slice(box[2],box[3]))



    lons = da.longitude
    lats = da.latitude

    #ipdb.set_trace()


    low_press = 850
    up_press = 200
    mid_press = 550

    gp = da['z'].mean('time') / 9.81

    gp_high = da['z'].sel(level=up_press)/9.81
    w_mid = da['w'].sel(level=mid_press)

    low_z = gp.sel(level=low_press) > 1400
    mid_z = gp.sel(level=mid_press) > 5500

    #ipdb.set_trace()

    tlow = da['t'].sel(level=low_press).where(low_z)-273.15
    qlow = da['q'].sel(level=low_press).where(low_z)*1000

    tmid = da['t'].sel(level=mid_press)-273.15#.where(mid_z)-273.15
    qmid = da['q'].sel(level=mid_press)*1000#.where(mid_z)*1000

    theta_low = u_met.theta_e(low_press,tlow, qlow)
    theta_high = u_met.theta_e(mid_press, tmid, qmid)

    theta_e = theta_low - theta_high

    u600 = da['u'].sel(level=up_press)#.where(mid_z)
    v600 = da['v'].sel(level=up_press)#.where(mid_z)

    u600.name = 'u200'

    ws600 = u_met.u_v_to_ws_wd(u600, v600)

    u800 = da['u'].sel(level=mid_press) # 200-500 shear
    v800 = da['v'].sel(level=mid_press)

    shear_u = u600-u800
    shear_v = v600-v800

    ws_shear = u_met.u_v_to_ws_wd(shear_u.values, shear_v.values)

    ws_600 = u600.copy(deep=True)
    ws_600.name = 'ws'

    ws_600.values = ws600[0]

    shear = u600.copy(deep=True)
    shear.name = 'shear'
    shear.values = ws_shear[0]


    vars = ['t600', 'q600', 'shear', 'q550', 'u200', 'v200', 'u550', 'v550', 'gp', 'w550']
    data = [tmid, qmid, shear, qmid, u600, v600, u800, v800, gp_high, w_mid]

    months= [11]#,6]#,2,3,4,5,6,7,8,9,10,11,12]#[3,4,5,6,9,10,11]#,4,5,6,9,10,11#,4,5,6,9,10,11,(3,5), (9,11)]#, 10,5,9]#[(12,2)]#[1,2,3,4,5,6,7,8,9,10,11,12]# #,2,3,11,12]#[(12,2)]#[1,2,3,4,5,6,7,8,9,10,11,12]# #,2,3,11,12]

    #ipdb.set_trace()
    for m in months:

        dic = {}
        for v in vars:
            dic[v] = 0

        if type(m)==int:
            m = [m]

        for v, dat in zip(vars,data):
            print('Doing ', v)
            dic[v] = get_trend(dat, m, sig=True, wilks=False, method='mk')

        if len(m) == 1:
            fp = fpath + 'low_ERA5_trend_synop_HU_poly_quatro_'+str(m[0]).zfill(2)+'.png'
        else:
            fp = fpath + 'low_ERA5_trend_synop_HU_poly_quatro_' + str(m[0]).zfill(2) +'-'+ str(m[1]).zfill(2) + '.png'
        #ipdb.set_trace()

        map = shear.salem.get_map()
        # ti_da = t2d.salem.transform(ti_da)

        f = plt.figure(figsize=(12,5), dpi=300)

        # transform their coordinates to the map reference system and plot the arrows
        xo, yo = map.grid.transform(shear.longitude.values, shear.latitude.values,
                                    crs=shear.salem.grid.proj)

        xH, yH = map.grid.transform(-77.52, -9.52,
                                    crs=shear.salem.grid.proj)

        xx, yy = np.meshgrid(xo, yo)

        ss = 3
        #Quiver only every 7th grid point
        u = (dic['u200'][0]).values[1::ss, 1::ss] # 200hpa
        v = (dic['v200'][0]).values[1::ss, 1::ss]

        #Quiver only every 7th grid point
        uu = (dic['u200'][1]).values[1::ss, 1::ss] # 200mean
        vv = (dic['v200'][1]).values[1::ss, 1::ss]

        u500 = (dic['u550'][0]).values[1::ss, 1::ss]
        v500 = (dic['v550'][0]).values[1::ss, 1::ss]

        #Quiver only every 7th grid point
        uu500 = (dic['u550'][1]).values[1::ss, 1::ss]
        vv500 = (dic['v550'][1]).values[1::ss, 1::ss]

        xx = xx[1::ss, 1::ss]
        yy = yy[1::ss, 1::ss]

        ax1 = f.add_subplot(121)
        map.set_data((dic['t600'][0]).values, interp='linear')  # interp='linear'

        map.set_contour(((dic['t600'][1]).values).astype(np.float64), interp='linear', colors='k', linewidths=0.5)
        map.set_plot_params(cmap='RdBu_r', extend='both', levels=[-0.5,-0.4,-0.3,-0.2,0.2,0.3,0.4,0.5])  #levels=[-0.5,-0.4,-0.3,-0.2,0.2,0.3,0.4,0.5], ,  levels=np.arange(-0.5,0.51,0.1),
        ax1.plot(xH, yH, 'ko', markersize=7)
        ax1.plot(xH, yH, 'wo', markersize=5)
        ax1.text(xH, yH+10, 'Huaraz', fontsize=11)
        #map.set_contour((t2_mean.values).astype(np.float64), interp='linear', colors='k', linewidths=0.5, levels=np.linspace(800,925,8))
        #map.set_plot_params(levels=[-0.5,-0.4,-0.3,-0.2,-0.1,-0.05,-0.02, 0.02,0.05,0.1,0.2,0.3,0.4,0.5], cmap='RdBu_r', extend='both')  # levels=np.arange(-0.5,0.51,0.1),
        cc = ['b', 'r', 'k']
        tags = ['1985-1994', '1995-2004', '2005-2014']
        for ids, pp in enumerate((dic['gp'][2])):
            xh, yh = map.grid.transform(pp[0], pp[1],
                                        crs=shear.salem.grid.proj)
            if ids < 2:
                ax1.plot(xh-80, yh-130, color=cc[ids], marker='x', markersize=10, label=tags[ids])
            else:
                ax1.plot(xh, yh, color=cc[ids], marker='x', markersize=10, label=tags[ids])
            #ax1.text(xh + 10, yh, tags[ids], fontsize=9)
        plt.legend(loc='upper left')
        pdic = map.visualize(ax=ax1, title='550hpa T trend (shading) and mean (contour)', cbar_title='K decade-1')
        contours = pdic['contour'][0]
        plt.clabel(contours, inline=True, fontsize=7, fmt='%1.1f')

        ax3 = f.add_subplot(122)
        map.set_data((dic['u200'][0]).values, interp='linear')  # interp='linear'
        #map.set_contour(u6_mean.values, interp='linear', colors='k')
        map.set_contour()
        ax3.plot(xH, yH, 'ko', markersize=7)
        ax3.plot(xH, yH, 'wo', markersize=5)
        map.set_plot_params(levels=[-3,-2.5,-2,-1.5,-1,-0.5,0.5,1,1.5,2,2.5,3], cmap='RdBu_r', extend='both')  # levels=np.arange(-0.5,0.51,0.1)
        map.visualize(ax=ax3, title='200hPa wind trend (shading) and mean (vectors)', cbar_title='m s-1 decade-1')
        qu = ax3.quiver(xx, yy, uu, vv, scale=100, width=0.002)

        cc = ['b', 'r', 'k']
        tags = ['1985-1994', '1995-2004', '2005-2014']
        for ii in range(3):
            use = (dic['u200'][2]).values[ii,:]
            for id, xpos in enumerate(xo):
                #ipdb.set_trace()
                ax3.plot(xpos, yo[use[id]], color=cc[ii], marker='o', markersize=4, label=tags[ii])
                # if id == 15:
                #     ax3.text(xpos + 10, yo[use[id]], tags[ii], fontsize=9)
        #plt.legend(loc='upper right')

        qk = plt.quiverkey(qu, 0.9, 0.03, 4, '4 m s$^{-1}$',
                           labelpos='E', coordinates='figure')

        ax1.text(14,-60, '(x) "Bolivian High centre" (max. 200hPa geopotential)')

        ax3.text(14,-60, '(o) Area of easterlies')

        plt.tight_layout()
        plt.savefig(fp)
        plt.close('all')
