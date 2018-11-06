import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from utils import u_darrays
import pdb
from utils import constants as cnst
import salem
from utils import u_statistics as us
from scipy import stats
import numpy.ma as ma
import pickle as pkl


def corr_box():
    srfc = cnst.ERA_MONTHLY_SRFC_SYNOP
    pl = cnst.ERA_MONTHLY_PL_SYNOP
    mcs = cnst.GRIDSAT + 'aggs/gridsat_WA_-70_monthly_count.nc'

    fpath = cnst.network_data + 'figs/CLOVER/months/'

    dicm = pkl.load(open(cnst.network_data + 'data/CLOVER/saves/storm_frac_synop12UTC.p', 'rb'))
    dicmean = pkl.load(open(cnst.network_data + 'data/CLOVER/saves/storm_frac_mean_synop12UTC.p', 'rb'))

    mcsbox = cnst.GRIDSAT + 'aggs/box_13W-13E-4-8N_meanT-50_from5000km2.nc'
    mcs_temp = xr.open_dataset(mcsbox)
    mcs_temp = mcs_temp['tir']

    da = xr.open_dataset(pl)
    da = u_darrays.flip_lat(da)
    da = da.sel(longitude=slice(-18, 40), latitude=slice(0, 25))  # latitude=slice(36, -37))
    da2 = xr.open_dataset(srfc)
    da2 = u_darrays.flip_lat(da2)
    da2 = da2.sel(longitude=slice(-18, 40), latitude=slice(0, 25))  # latitude=slice(36, -37))
    da3 = xr.open_dataset(mcs)
    da3 = da3.sel(lon=slice(-18, 40), lat=slice(0, 25))

    lons = da.longitude
    lats = da.latitude

    q = da['q'].sel(level=slice(850,925)).mean('level')
    t2 = da['t'].sel(level=slice(850, 925)).mean('level')
    u925 = da['u'].sel(level=slice(850, 925)).mean('level')
    u600 = da['u'].sel(level=slice(550,700)).mean('level')

    shear = u600-u925 # u600-

    q.values = q.values * 1000

    tir = da3['tir']
    tir = t2.salem.lookup_transform(tir)

    months = np.arange(1, 13)
    months = [3]

    def array_juggling(data, month, hour=None):

        m = month

        if hour != None:
            data = data[(data['time.month'] == m) & (data['time.hour'] == hour) & (data['time.year'] >= 1983)& (
            data['time.year'] <= 2017)]
        else:
            data = data[(data['time.month'] == m) & (data['time.year'] >= 1983) & (data['time.year'] <= 2017)]
        data_years = data.groupby('time.year').mean(axis=0)

        data_mean = data.mean(axis=0)

        diff = xr.DataArray(data_years.values[1::, :, :] - data_years.values[0:-1, :, :],
                            coords=[data_years.year[1::], data.latitude, data.longitude], dims=['year','latitude', 'longitude'] )
        # diff = xr.DataArray(data_years.values, coords=[data_years.year, data.latitude, data.longitude],
        #                     dims=['year', 'latitude', 'longitude'])
        # unstack back to lat lon coordinates
        return diff, data_mean


    def corr(a, b, bsingle=None):
        ds = xr.Dataset()
        ds['pval'] = a.copy(deep=True).sum('year') * np.nan
        ds['r'] = a.copy(deep=True).sum('year') * np.nan
        ds['slope'] = a.copy(deep=True).sum('year') * np.nan

        if bsingle:
            bb = b
        else:
            bb = b.sel(latitude=slice(4, 8), longitude=slice(-10, 11)).mean(dim=['latitude', 'longitude'])

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
        tirdiff, tiryear = array_juggling(tir, m)

        mcs_month = mcs_temp[mcs_temp['time.month'] == m]

        #tirdiff = mcs_month.values#[1::]-mcs_month.values[0:-1]

        qcorr = corr(qdiff, tirdiff, bsingle=False)
        shearcorr = corr(shdiff, tirdiff, bsingle=False)
        tcorr = corr(t2diff, tirdiff, bsingle=False)

        # pthresh = us.fdr_threshold(qcorr['pval'].values[np.isfinite(qcorr['pval'].values)], alpha=0.05)
        # print(pthresh)
        pthresh = 0.05
        #cloud['slope'].values[cloud['pval'].values > pthresh] = np.nan
        qcorr['r'].values[qcorr['pval'].values > pthresh] = 0

        # pthresh = us.fdr_threshold(shearcorr['pval'].values[np.isfinite(shearcorr['pval'].values)], alpha=0.05)
        # print(pthresh)
        shearcorr['r'].values[shearcorr['pval'].values > pthresh] = 0

        # pthresh = us.fdr_threshold(tcorr['pval'].values[np.isfinite(tcorr['pval'].values)], alpha=0.05)
        # print(pthresh)
        tcorr['r'].values[tcorr['pval'].values > pthresh] = 0


        fp = fpath + 'corr_SYNOP_box_detrended' + str(m).zfill(2) + '.png'
        map = shear.salem.get_map()

        f = plt.figure(figsize=(13,7), dpi=300)
        ax1 = f.add_subplot(221)
        # map.set_shapefile(rivers=True)
        # bla = ma.masked_invalid(tcorr['r'].values)

        map.set_data(tcorr['r'].values, interp='linear')  # interp='linear'
        contours = map.set_contour(t2year-273.15, interp='linear', levels=np.arange(24,37,4), cmap='inferno')

        #plt.clabel(contours, inline=True, fontsize=7, fmt='%1.1f')

        map.set_plot_params(cmap='RdBu', extend='both',levels=[-0.7,-0.6,-0.5,-0.4,0.4,0.5,0.6,0.7])  # levels=np.arange(-0.5,0.51,0.1),
        map.visualize(ax=ax1, title='925hP temperature')

        ax2 = f.add_subplot(222)
        map.set_data(qcorr['r'],interp='linear')  # interp='linear'
        map.set_contour(qyear,interp='linear', levels=np.arange(5,19,3), cmap='inferno')

        map.set_plot_params(cmap='RdBu', extend='both',levels=[-0.7,-0.6,-0.5,-0.4,0.4,0.5,0.6,0.7])  # levels=np.arange(-0.5,0.51,0.1),
        map.visualize(ax=ax2, title='925hPa Spec. humidity')

        ax3 = f.add_subplot(223)
        map.set_data(shearcorr['r'], interp='linear')  # interp='linear'
        map.set_contour(sheyear, interp='linear', levels=np.arange(-10,1,3), cmap='inferno')
        map.set_plot_params(cmap='RdBu_r', extend='both', levels=[-0.7,-0.6,-0.5,-0.4,0.4,0.5,0.6,0.7])  # levels=np.arange(-0.5,0.51,0.1)
        map.visualize(ax=ax3, title='600-925hPa Zonal wind shear')

        ax4 = f.add_subplot(224)
        map.set_contour(dicmean[m], interp='linear', levels=[0.1,0.5,1,2.5], cmap='inferno')
        map.set_data(dicm[m])  #


        map.set_plot_params(cmap='viridis', extend='both', levels=np.arange(10,51,10))  # levels=np.arange(20,101,20)  #np.arange(20,101,20)
        map.visualize(ax=ax4, title='-70C cloud cover change', cbar_title='$\%$ decade-1')

        plt.tight_layout()
        plt.savefig(fp)
        plt.close('all')