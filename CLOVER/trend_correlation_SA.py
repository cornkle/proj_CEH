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


c_box = [15,32,-20,-15]#Middle[15,32,-20,-15]# DJF[17, 25, -28, -22]

def corr_box():
    srfc = cnst.ERA_MONTHLY_SRFC_SYNOP
    pl = cnst.ERA_MONTHLY_PL_SYNOP
    mcs = cnst.GRIDSAT + 'aggs/gridsat_WA_-65_monthly_count_-40base_15-21UTC_1000km2.nc'  # -70count

    fpath = cnst.network_data + 'figs/CLOVER/months/'

    dicm = pkl.load(open(cnst.network_data + 'data/CLOVER/saves/storm_frac_synop12UTC_SA.p', 'rb'))
    dicmean = pkl.load(open(cnst.network_data + 'data/CLOVER/saves/storm_frac_mean_synop12UTC_SA.p', 'rb'))

    # mcsbox = cnst.GRIDSAT + 'aggs/SAboxWest_meanT-40_1000km2.nc' # box_13W-13E-4-8N_meanT-50_from5000km2.nc'
    # mcs_temp = xr.open_dataset(mcsbox)
    # mcs_temp = mcs_temp['tir']

    box=[5,55,-36,0]#[-18,55,-35,35]#[-10,55,-35,0]

    da = xr.open_dataset(pl)
    da = u_darrays.flip_lat(da)
    da = da.sel(longitude=slice(box[0], box[1]), latitude=slice(box[2], box[3]))  # latitude=slice(36, -37))
    da2 = xr.open_dataset(srfc)
    da2 = u_darrays.flip_lat(da2)
    da2 = da2.sel(longitude=slice(box[0], box[1]), latitude=slice(box[2], box[3]))  # latitude=slice(36, -37))
    da3 = xr.open_dataset(mcs)
    da3 = da3.sel(lon=slice(box[0], box[1]), lat=slice(box[2], box[3]))

    lons = da.longitude
    lats = da.latitude

    press = da2['sp']
    press = press[press['time.hour'] == 12]
    press.values = press.values*1000
    low_press = 850
    up_press = 550

    q = da['q'].sel(level=slice(low_press-50, low_press)).mean('level')
    q = q[q['time.hour']==12]
    t2d = da2['t2m']#['t2m']
    #t2d = da['t'].sel(level=slice(800, 850)).mean('level')
    t2d = t2d[t2d['time.hour']==12]

    u600 = da['u'].sel(level=slice(up_press-50, up_press)).mean('level')
    u600 = u600[u600['time.hour']==12]
    v600 = da['v'].sel(level=slice(up_press-50, up_press)).mean('level')
    v600 = v600[v600['time.hour']==12]

    ws600 = u_met.u_v_to_ws_wd(u600, v600)

    u800 = da['u'].sel(level=slice(low_press-50, low_press)).mean('level')
    u800 = u800[u800['time.hour']==12]

    v800 = da['v'].sel(level=slice(low_press-50, low_press)).mean('level')
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

    q.values = q.values * 1000

    tir = da3['tir']
    tir = t2d.salem.lookup_transform(tir)

    def array_juggling(data, month, hour=None):

        m = month

        if hour is not None:
            if len(month) > 1:

                data = data[((data['time.month'] >= month[0]) | (data['time.month'] <= month[1])) & (
                            data['time.hour'] == hour) & (data['time.year'] >= 1983) & (data['time.year'] <= 2017)]
            else:

                data = data[
                    (data['time.month'] == month[0]) & (data['time.hour'] == hour) & (data['time.year'] >= 1983) & (
                                data['time.year'] <= 2017)]
        else:
            if len(month) > 1:
                data = data[((data['time.month'] >= month[0]) | (data['time.month'] <= month[1])) & (
                            data['time.year'] >= 1983) & (data['time.year'] <= 2017)]
            else:
                data = data[
                    (data['time.month'] == month[0]) & (data['time.year'] >= 1983) & (data['time.year'] <= 2017)]

        data_years = data.groupby('time.year').mean(axis=0)

        data_mean = data.mean(axis=0)

        # diff = xr.DataArray(data_years.values[1::, :, :] - data_years.values[0:-1, :, :],
        #                     coords=[data_years.year[1::], data.latitude, data.longitude], dims=['year','latitude', 'longitude'] )
        diff = xr.DataArray(data_years.values, coords=[data_years.year, data.latitude, data.longitude],
                            dims=['year', 'latitude', 'longitude'])
        # unstack back to lat lon coordinates
        return diff, data_mean


    def corr(a, b, bsingle=None, c_box=None):
        ds = xr.Dataset()
        ds['pval'] = a.copy(deep=True).sum('year') * np.nan
        ds['r'] = a.copy(deep=True).sum('year') * np.nan
        ds['slope'] = a.copy(deep=True).sum('year') * np.nan

        corr_box = c_box

        if bsingle:
            bb = b
        else:
            bb = b.sel(latitude=slice(corr_box[2], corr_box[3]), longitude=slice(corr_box[0], corr_box[1])).mean(dim=['latitude', 'longitude'])

        for lat in a.latitude.values:
            for lon in a.longitude.values:
                aa = a.sel(latitude=lat, longitude=lon)
                if bsingle:
                    r, p = stats.pearsonr(aa.values, bb)

                    pf = np.polyfit(aa.values, bb, 1)
                else:
                    r, p = stats.pearsonr(aa.values, bb.values)
                    pf = np.polyfit(aa.values, bb.values, 1)


                slope = pf[0]

                if (np.nansum(aa.values == 0) >= 10):
                    p = np.nan
                    r = np.nan

                ds['r'].loc[{'latitude': lat, 'longitude': lon}] = r
                ds['pval'].loc[{'latitude': lat, 'longitude': lon}] = p
                ds['slope'].loc[{'latitude': lat, 'longitude': lon}] = slope

        return ds



    box_dic = {

        1 : [16,30,-25,-20],
        2 : [12,28,-10,3],
        3: [16,25,-23,-18],
        4: [16,25,-23,-18],
        10 : [14,28,-15,-5],
        11 : [16,30,-20,-10],
        12 : [16,30,-22,-12],
        (11,1) : [18,30,-23,-18]
    }


    months = [(11,1),2,3,10]

    for m in months:

        c_box = box_dic[m]

        if type(m)==int:
            m = [m]

        t2diff, t2year = array_juggling(t2d, m) #
        qdiff, qyear = array_juggling(q, m) #, hour=12
        shdiff, sheyear = array_juggling(shear, m) #, hour=12
        vdiff, vyear = array_juggling(v800, m)  # , hour=12
        udiff, uyear = array_juggling(u800, m)  # , hour=12
        tirdiff, tiryear = array_juggling(tir, m)  # average frequency change

        #mcs_month = mcs_temp[mcs_temp['time.month'] == m] # meanT box average change

        #tirdiff = mcs_month.values[1::]-mcs_month.values[0:-1]

        bs = False
        try:
            qcorr = corr(qdiff, tirdiff, bsingle=bs, c_box=c_box)
        except:
            continue
        shearcorr = corr(shdiff, tirdiff, bsingle=bs, c_box=c_box)
        tcorr = corr(t2diff, tirdiff, bsingle=bs, c_box=c_box)

        dicm[m[0]].values[dicm[m[0]].values==0] = np.nan

        print('plot')

        if len(m) == 1:
            fp = fpath + 'corr_mid_-70C_synop_-50base_linear_850T'+str(m[0]).zfill(2)+'.png'
        else:
            fp = fpath + 'corr_mid_-70C_synop_-50base_linear_850T' + str(m[0]).zfill(2) +'-'+ str(m[1]).zfill(2) + '.png'


        map = shear.salem.get_map()

        xx, yy = map.grid.transform(shear.longitude.values, shear.latitude.values,
                                    crs=shear.salem.grid.proj)

        xx, yy = np.meshgrid(xx, yy)

        # Quiver only every 7th grid point
        u = uyear.values[1::2, 1::2]
        v = vyear.values[1::2, 1::2]

        xx = xx[1::2, 1::2]
        yy = yy[1::2, 1::2]

        #ipdb.set_trace()
        f = plt.figure(figsize=(15,8), dpi=350)
        ax1 = f.add_subplot(221)

        map.set_data(tcorr['r'], interp='linear')  # interp='linear'
        map.set_contour(t2year-273.15, interp='linear', levels=np.arange(24,37,4),colors='k', linewidths=0.5)


        map.set_plot_params(cmap='RdBu_r', extend='both',levels=[-0.7,-0.6,-0.5,-0.4,-0.3,0.3,0.4,0.5,0.6,0.7])  # levels=np.arange(-0.5,0.51,0.1),
        dic = map.visualize(ax=ax1, title='2m temperature corr. | contours: mean T', cbar_title='K decade-1')
        contours = dic['contour'][0]
        plt.clabel(contours, inline=True, fontsize=7, fmt='%1.1f')

        ax2 = f.add_subplot(222)
        map.set_data(qcorr['r'],interp='linear')  # interp='linear'
        map.set_contour(qyear,interp='linear', levels=np.arange(5,19,3), colors='k', linewidths=0.5)

        map.set_plot_params(cmap='RdBu_r', extend='both',levels=[-0.7,-0.6,-0.5,-0.4,-0.3,0.3,0.4,0.5,0.6,0.7])  # levels=np.arange(-0.5,0.51,0.1),
        dic = map.visualize(ax=ax2, title='800hPa Spec. humidity corr. | contours: mean q', cbar_title='g kg-1 decade-1')
        contours = dic['contour'][0]
        plt.clabel(contours, inline=True, fontsize=7, fmt='%1.1f')

        ax3 = f.add_subplot(223)
        # map.set_data(shearcorr['r'], interp='linear')  # interp='linear'
        # map.set_plot_params(cmap='RdBu_r', extend='both', levels=[-0.7,-0.6,-0.5,-0.4, -0.3,0.3,0.4,0.5,0.6,0.7])
        #map.set_contour(sheyear, interp='linear', levels=np.arange(-10,1,6), colors='k', linewidths=0.5)
        map.set_data(tcorr['r'], interp='linear')  # interp='linear'
        map.set_plot_params(cmap='RdBu_r', extend='both',levels=[-0.7,-0.6,-0.5,-0.4,-0.3,0.3,0.4,0.5,0.6,0.7])

        qu = ax3.quiver(xx, yy, u, v, scale=50, width=0.002)
        qk = plt.quiverkey(qu, 0.4, 0.03, 4, '4 m s$^{-1}$',
                           labelpos='E', coordinates='figure')

          # levels=np.arange(-0.5,0.51,0.1)
        dic = map.visualize(ax=ax3, title='800-500hPa wind shear corr., mean 500hPa wind vectors', cbar_title='m s-1 decade-1')
        contours = dic['contour'][0]
        plt.clabel(contours, inline=True, fontsize=7, fmt='%1.1f')

        ax4 = f.add_subplot(224)
        map.set_contour(dicmean[m[0]], interp='linear', levels=[0.1,0.5,1,2.5], colors='k', linewidths=0.5)

        map.set_data(dicm[m[0]])  #
        #ax4.axhspan(-26,18)  #[15,25,-26,-18]
        coord = c_box#[17, 25, -28, -22]
        geom = shpg.box(coord[0], coord[2], coord[1], coord[3])
        map.set_geometry(geom, zorder=99, color='darkorange', linewidth=3, linestyle='--', alpha=0.3)
        # ax4.axvline(x=25, ymin=-26, ymax=-18)
        # ax4.axhline(y=-26, xmin=15, xmax=25)
        # ax4.axhline(y=-18, xmin=15, xmax=25)

        map.set_plot_params(cmap='viridis', extend='both', levels=np.arange(10,41,10))  # levels=np.arange(20,101,20)  #np.arange(20,101,20)
        dic = map.visualize(ax=ax4, title='-65C cloud cover change | >1000km2 -40C', cbar_title='$\%$ decade-1')
        contours = dic['contour'][0]
        plt.clabel(contours, inline=True, fontsize=7, fmt='%1.1f')

        plt.tight_layout()
        plt.savefig(fp)
        plt.close('all')
