import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from utils import u_darrays
import ipdb
from utils import constants as cnst
import salem
from utils import u_statistics as us
from scipy import stats
import numpy.ma as ma
import pickle as pkl
import shapely.geometry as shpg


def corr_box():
    srfc = cnst.ERA_MONTHLY_SRFC_SYNOP
    pl = cnst.ERA_MONTHLY_PL_SYNOP # _SYNOP
    mcs = cnst.GRIDSAT + 'aggs/gridsat_WA_-65_monthly_count_-40base_1000km2.nc'  # -70count

    fpath = cnst.network_data + 'figs/CLOVER/months/'


    box=[-10,55,-35,0]#[-18,55,-35,35]#[-10,55,-35,0]


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


    q = da['q'].sel(level=slice(800)).mean('level')
    #q = q[q['time.hour']==12]
    t2 = da2['t2m']#.sel(level=slice(800)).mean('level')
    #t2 = t2[t2['time.hour']==12]
    u925 = da['u'].sel(level=slice(800)).mean('level')
    #u925 = u925[u925['time.hour']==12]
    u600 = da['u'].sel(level=slice(500,550)).mean('level')
    v600 = da['v'].sel(level=slice(500,550)).mean('level')

    shear = u925 # u600-

    q.values = q.values * 1000

    tir = da3['tir']
    tir = t2.salem.lookup_transform(tir)


    months = np.arange(1, 13)
    months = [1,12]

    def array_juggling(data, month, hour=None):

        m = month

        if hour != None:
            data = data[(data['time.month'] == m) & (data['time.hour'] == hour) & (data['time.year'] >= 1983)& (
            data['time.year'] <= 2017)]
        else:
            data = data[(data['time.month'] == m) & (data['time.year'] >= 1983) & (data['time.year'] <= 2017)]
        data_years = data.groupby('time.year').mean(axis=0)

        data_mean = data.mean(axis=0)

        # diff = xr.DataArray(data_years.values[1::, :, :] - data_years.values[0:-1, :, :],
        #                     coords=[data_years.year[1::], data.latitude, data.longitude], dims=['year','latitude', 'longitude'] )
        diff = xr.DataArray(data_years.values, coords=[data_years.year, data.latitude, data.longitude],
                            dims=['year', 'latitude', 'longitude'])
        # unstack back to lat lon coordinates
        return diff, data_mean


    def corr(a, b, bsingle=None):
        ds = xr.Dataset()
        ds['pval'] = a.copy(deep=True).sum('year') * np.nan
        ds['r'] = a.copy(deep=True).sum('year') * np.nan
        ds['slope'] = a.copy(deep=True).sum('year') * np.nan

        #corr_box = [-10,11,4.5,8]
        corr_box = [16,24,-24,-18]#West: [13,25,-23,-10]

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

                # if (np.nansum(aa.values == 0) >= 10):
                #     p = np.nan
                #     r = np.nan

                ds['r'].loc[{'latitude': lat, 'longitude': lon}] = r
                ds['pval'].loc[{'latitude': lat, 'longitude': lon}] = p
                ds['slope'].loc[{'latitude': lat, 'longitude': lon}] = slope

        return ds

    for m in months:

        t2diff, t2year = array_juggling(t2, m, hour=12) #
        qdiff, qyear = array_juggling(q, m, hour=12) #, hour=12
        shdiff, sheyear = array_juggling(shear, m, hour=12) #, hour=12
        vdiff, vyear = array_juggling(v600, m, hour=12)  # , hour=12
        udiff, uyear = array_juggling(u600, m, hour=12)  # , hour=12
        tirdiff, tiryear = array_juggling(tir, m)  # average frequency change

        mcs_month = mcs_temp[mcs_temp['time.month'] == m] # meanT box average change

        #tirdiff = mcs_month.values[1::]-mcs_month.values[0:-1]

        bs = False
        try:
            qcorr = corr(qdiff, tirdiff, bsingle=bs)
        except:
            continue
        shearcorr = corr(shdiff, tirdiff, bsingle=bs)
        tcorr = corr(t2diff, tirdiff, bsingle=bs)


        # pthresh = us.fdr_threshold(qcorr['pval'].values[np.isfinite(qcorr['pval'].values)], alpha=0.05)
        # print(pthresh)
        pthresh = 0.05
        #cloud['slope'].values[cloud['pval'].values > pthresh] = np.nan
        #qcorr['r'].values[qcorr['pval'].values > pthresh] = 0

        # pthresh = us.fdr_threshold(shearcorr['pval'].values[np.isfinite(shearcorr['pval'].values)], alpha=0.05)
        # print(pthresh)
        #shearcorr['r'].values[shearcorr['pval'].values > pthresh] = 0

        # pthresh = us.fdr_threshold(tcorr['pval'].values[np.isfinite(tcorr['pval'].values)], alpha=0.05)
        # print(pthresh)
        #tcorr['r'].values[tcorr['pval'].values > pthresh] = 0

        dicm[m].values[dicm[m].values==0] = np.nan

        print('plot')
        fp = fpath + 'corr_box_SYNOP_SAWest_-50base_' + str(m).zfill(2) + '.png'
        map = shear.salem.get_map()
        #xx, yy = map.grid.xy_coordinates

        # transform their coordinates to the map reference system and plot the arrows
        xx, yy = map.grid.transform(shear.longitude.values, shear.latitude.values,
                                     crs=shear.salem.grid.proj)

        xx, yy = np.meshgrid(xx,yy)

        #ipdb.set_trace()
        f = plt.figure(figsize=(15,8), dpi=200)
        ax1 = f.add_subplot(221)
        # map.set_shapefile(rivers=True)
        # bla = ma.masked_invalid(tcorr['r'].values)

        map.set_data(tcorr['r'], interp='linear')  # interp='linear'
        contours = map.set_contour(t2year-273.15, interp='linear', levels=np.arange(24,37,4), cmap='inferno')

        #plt.clabel(contours, inline=True, fontsize=7, fmt='%1.1f')

        map.set_plot_params(cmap='RdBu_r', extend='both',levels=[-0.7,-0.6,-0.5,-0.4,0.4,0.5,0.6,0.7])  # levels=np.arange(-0.5,0.51,0.1),
        map.visualize(ax=ax1, title='2m temperature')

        ax2 = f.add_subplot(222)
        map.set_data(qcorr['r'],interp='linear')  # interp='linear'
        map.set_contour(qyear,interp='linear', levels=np.arange(5,19,3), cmap='inferno')

        map.set_plot_params(cmap='RdBu_r', extend='both',levels=[-0.7,-0.6,-0.5,-0.4,0.4,0.5,0.6,0.7])  # levels=np.arange(-0.5,0.51,0.1),
        map.visualize(ax=ax2, title='800hPa Spec. humidity')

        ax3 = f.add_subplot(223)
        map.set_data(shearcorr['r'], interp='linear')  # interp='linear'
        map.set_plot_params(cmap='RdBu_r', extend='both', levels=[-0.7,-0.6,-0.5,-0.4, 0.4,0.5,0.6,0.7])
        map.set_contour(sheyear, interp='linear', levels=np.arange(-10,1,6), cmap='inferno')
        cs = ax3.contour(sheyear, interp='linear', levels=np.arange(-10,1,6), cmap='inferno')
        plt.clabel(cs, inline=1, fontsize=10)

        # Quiver only every 7th grid point
        u = uyear[4::7, 4::7]
        v = vyear[4::7, 4::7]

        xx = xx[4::7, 4::7]
        yy = yy[4::7, 4::7]


        qu = ax3.quiver(xx, yy, u.values, v.values, scale=50)

          # levels=np.arange(-0.5,0.51,0.1)
        map.visualize(ax=ax3, title='500-800hPa Zonal wind shear')

        ax4 = f.add_subplot(224)
        map.set_contour(dicmean[m], interp='linear', levels=[0.1,0.5,1,2.5], cmap='inferno')

        map.set_data(dicm[m])  #
        #ax4.axhspan(-26,18)  #[15,25,-26,-18]
        coord = [16,24,-24,-18]
        geom = shpg.box(coord[0], coord[2], coord[1], coord[3])
        map.set_geometry(geom, zorder=99, color='darkorange', linewidth=3, linestyle='--', alpha=0.3)
        # ax4.axvline(x=25, ymin=-26, ymax=-18)
        # ax4.axhline(y=-26, xmin=15, xmax=25)
        # ax4.axhline(y=-18, xmin=15, xmax=25)

        map.set_plot_params(cmap='viridis', extend='both', levels=np.arange(10,51,10))  # levels=np.arange(20,101,20)  #np.arange(20,101,20)
        map.visualize(ax=ax4, title='-65C cloud cover change', cbar_title='$\%$ decade-1')

        plt.tight_layout()
        plt.savefig(fp)
        plt.close('all')
