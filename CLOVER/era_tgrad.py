import pandas as pd
import numpy as np
import xarray as xr
import pdb
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy
import salem
from utils import u_plot as up, u_mann_kendall as mk, constants_lappi as cnst
from utils import u_darrays
from scipy import stats
from matplotlib import lines
from utils import u_statistics as ustats


def t_trend():
    #file = '/users/global/cornkle/data/ERA-I monthly/ERA-WA-Monthly-2mTemp.nc'
    file = '/localscratch/wllf030/cornkle/ERA-I/monthly/monthly_1979-2017_srfc.nc'

    fpath = '/users/global/cornkle/figs/CLOVER/months/'


    dam = xr.open_dataset(file)
    dam = dam['t2m']
    months = np.arange(1,13)

    for m in months:

        da = dam[(dam['time.month']==m)]
        da = da.sel(longitude=slice(-18,51), latitude=slice(36, -37))
        da = da.groupby('time.year').mean(axis=0)

        lons = da.longitude
        lats = np.flip(da.latitude.values, axis=0)

        # stack lat and lon into a single dimension called allpoints
        stacked = da.stack(allpoints=['latitude','longitude'])
        # apply the function over allpoints to calculate the trend at each point
        trend = stacked.groupby('allpoints').apply(u_darrays.linear_trend)
        # unstack back to lat lon coordinates
        trend_unstacked = trend.unstack('allpoints')

        trend_unstacked = trend_unstacked*10. # warming over decade
        da2 = xr.DataArray(trend_unstacked, coords=[lats, lons], dims=['latitude', 'longitude'])

        fp = fpath + 'ttrend_WA'+str(m).zfill(2)+'.png'

        up.quick_map_salem(da2, levels=np.arange(-0.5,0.5,0.1), cmap='RdBu_r', save=fp)  #

        plt.close('all')


def t_trend_slice():
    #file = '/users/global/cornkle/data/ERA-I monthly/ERA-WA-Monthly-2mTemp.nc'
    file = '/localscratch/wllf030/cornkle/ERA-I/monthly/old/ERA-Int-Monthly-2mTemp.nc'

    fpath = '/users/global/cornkle/figs/CLOVER/months/'


    dam = xr.open_dataarray(file)
    lower = 9
    higher = 11

    da = dam[(dam['time.month']>=lower) & (dam['time.month']<=higher)]
    da = da.sel(longitude=slice(-18,51), latitude=slice(36, -37))
    da = da.groupby('time.year').mean(axis=0)

    lons = da.longitude
    lats = np.flip(da.latitude.values, axis=0)

    # define a function to compute a linear trend of a timeseries
    def linear_trend(x):

        #pf = np.polyfit(np.arange(len(x)), x, 1)
        pf, slope, int, p, ind = mk.test(np.arange(len(x)),x.squeeze().values, eps=0.001, alpha=0.01, Ha='upordown')

        # we need to return a dataarray or else xarray's groupby won't be happy

        if ind == 1:
            issig = slope
        else:
            issig = np.nan

        return xr.DataArray(issig, )

    # stack lat and lon into a single dimension called allpoints
    stacked = da.stack(allpoints=['latitude','longitude'])
    # apply the function over allpoints to calculate the trend at each point
    trend = stacked.groupby('allpoints').apply(linear_trend)
    # unstack back to lat lon coordinates
    trend_unstacked = trend.unstack('allpoints')

    trend_unstacked = trend_unstacked*10. # warming over decade
    da2 = xr.DataArray(trend_unstacked, coords=[lats, lons], dims=['latitude', 'longitude'])

    fp = fpath + 'ttrend_'+str(lower).zfill(2)+'-'+str(higher).zfill(2)+'.png'

    up.quick_map_salem(da2, vmin=-0.4, vmax=0.4, cmap='RdBu_r', save=fp)  #

    plt.close('all')




def t_trend_polyfit():
    #file = '/users/global/cornkle/data/ERA-I monthly/ERA-WA-Monthly-2mTemp.nc'
    file = '/localscratch/wllf030/cornkle/ERA-I/monthly/old/ERA-Int-Monthly-2mTemp.nc'

    fpath = '/users/global/cornkle/figs/CLOVER/months/'


    dam = xr.open_dataarray(file)
    lower = 9
    higher = 11

    da = dam[(dam['time.month']>=lower) & (dam['time.month']<=higher)]
    da = da.sel(longitude=slice(-18,51), latitude=slice(36, -37))
    da = da.groupby('time.year').mean(axis=0)

    lons = da.longitude
    lats = np.flip(da.latitude.values, axis=0)

    # define a function to compute a linear trend of a timeseries
    def linear_trend(x):

        #pf = np.polyfit(np.arange(len(x)), x, 1)
        pf, slope, int, p, ind = mk.test(np.arange(len(x)),x.squeeze().values, eps=0.001, alpha=0.01, Ha='upordown')

        # we need to return a dataarray or else xarray's groupby won't be happy

        if ind == 1:
            issig = slope
        else:
            issig = np.nan

        return xr.DataArray(issig, )

    # stack lat and lon into a single dimension called allpoints
    stacked = da.stack(allpoints=['latitude','longitude'])
    # apply the function over allpoints to calculate the trend at each point
    trend = stacked.groupby('allpoints').apply(linear_trend)
    # unstack back to lat lon coordinates
    trend_unstacked = trend.unstack('allpoints')

    trend_unstacked = trend_unstacked*10. # warming over decade
    da2 = xr.DataArray(trend_unstacked, coords=[lats, lons], dims=['latitude', 'longitude'])

    fp = fpath + 'ttrend_'+str(lower).zfill(2)+'-'+str(higher).zfill(2)+'.png'

    up.quick_map_salem(da2, vmin=-0.4, vmax=0.4, cmap='RdBu_r', save=fp)  #

    plt.close('all')




def t_mean():
    # file = '/users/global/cornkle/data/ERA-I monthly/ERA-WA-Monthly-2mTemp.nc'
    file = '/users/global/cornkle/data/ERA-I monthly/ERA-Int-Monthly-2mTemp.nc'

    fpath = '/users/global/cornkle/figs/gap_filling_Tgrad/months/'


    dam = xr.open_dataarray(file)
    months = np.arange(1, 13)

    for m in months:
        da = dam[(dam['time.month'] == m)]
        da = da.sel(longitude=slice(-18, 51), latitude=slice(36, -37))
        da = da.mean(axis=0)-273.15

        fp = fpath + 'tmean_' + str(m).zfill(2) + '.png'

        up.quick_map_salem(da, levels=np.arange(20,41,2), cmap='jet', save=fp)

def tgrad_shear_trend():

    srfc = cnst.ERA_MONTHLY_SRFC
    pl = cnst.ERA_MONTHLY_PL
    mcs = cnst.GRIDSAT + 'aggs/box_13W-13E-4-8N_meanT-50_from5000km2.nc'
    out = cnst.network_data + 'figs/CLOVER/'


    box = [-10,10,5.5,9]
    TNORTH = [-10,10,11,18]
    TSOUTH = [-10,10,5,9]

    dam = xr.open_dataset(srfc)
    dam = u_darrays.flip_lat(dam)
    dam = dam['t2m']
    tsouth = dam.sel(longitude=slice(TSOUTH[0], TSOUTH[1]), latitude=slice(TSOUTH[2],TSOUTH[3]))
    tnorth = dam.sel(longitude=slice(TNORTH[0], TNORTH[1]), latitude=slice(TNORTH[2],TNORTH[3]))

    da = xr.open_dataset(pl)
    da = u_darrays.flip_lat(da)
    da = da.sel(longitude=slice(box[0], box[1]), latitude=slice(box[2],box[3]))#latitude=slice(36, -37))

    u925 = da['u'].sel(level=925)#).mean(dim='level')  #slice(850
    u600 = da['u'].sel(level=650)
    qq = da['q'].sel(level=925)

    mcs_temp = xr.open_dataset(mcs)
    mcs_temp = mcs_temp['tir']
    months = [3,4,5,9,10]
    mnames = {3 : 'March', 4 : 'April', 5 : 'May', 9 : 'September', 10 : 'October'}

    for m in months:

        south = tsouth[(tsouth['time.month']==m)]
        north = tnorth[(tnorth['time.month']==m)]
        ulow = u925[(u925['time.month']==m)]
        uhigh = u600[(u600['time.month']==m)]
        mcs_month = mcs_temp[mcs_temp['time.month']==m]
        qmonth = qq[(qq['time.month']==m)]
        ##da = da.sel(longitude=slice(-18,51), latitude=slice(36, -37))

        south_peryear = south.groupby('time.year').mean('longitude').min('latitude')
        south_peryear5 = south_peryear.rolling(time=5, center=True).mean()
        north_peryear = north.groupby('time.year').mean('longitude').max('latitude')
        north_peryear5 = north_peryear.rolling(time=5, center=True).mean()

        u925_peryear = ulow.groupby('time.year').mean('longitude').max('latitude') #ulow.groupby('time.year').mean()
        u925_peryear5 = u925_peryear.rolling(time=5, center=True).mean()

        u600_peryear = uhigh.groupby('time.year').mean('longitude').min('latitude')#.mean() # ('latitude').min()
        u600_peryear5 = u600_peryear.rolling(time=5, center=True).mean()

        q_peryear = qmonth.groupby('time.year').mean('longitude').max('latitude')#.mean() # ('latitude').min()
        q_peryear5 = q_peryear.rolling(time=5, center=True).mean()


        tgrad = ((north_peryear-south_peryear)[4::])
        shear = (u600_peryear-u925_peryear)[4::] # -q_peryear[4::]#

        r = stats.pearsonr(shear.values[1::]-shear.values[0:-1],mcs_month.values[1::]-mcs_month.values[0:-1])
        tshear_cor = stats.pearsonr(shear.values[1::]-shear.values[0:-1],tgrad.values[1::]-tgrad.values[0:-1])

        sslope, sint = ustats.linear_trend(shear)
        mslope, mint = ustats.linear_trend(mcs_month)

        x = np.arange(0,len(shear))
        rr = r[0]
        f=plt.figure(figsize=(6,3))
        ax = f.add_subplot(111)
        ax.plot(tgrad.year, shear, 'x-', label='Zonal wind shear 600-925hPa', color='k')
        ax.plot(tgrad.year, sint + x*sslope, '--', color='k')
        ax.set_ylim(-14,-7)
        ax1 = ax.twinx()
        ax1.plot(mcs_month['time.year'],mcs_month, 'o-', label='Mean MCS temp.', color='r')
        ax1.plot(tgrad.year, mint + x*mslope, '--', color='r')
        mcsline = lines.Line2D([],[], color='r', label='Mean MCS temp.', linestyle='solid', marker='o')
        shearline = lines.Line2D([],[], color='k', label='Zonal wind shear 600-925hPa', linestyle='solid', marker='x', markersize=5)
        ax1.set_ylabel('degC')
        ax.set_ylabel('m s-1')
        ax.set_title(mnames[m]+' | Corr.:'+ str(np.round(rr, decimals=2)) + '| Tgrad/Shear corr: ' + str(np.round(tshear_cor[0], decimals=2)))
        if m==3:
            ax.legend(handles=[mcsline, shearline])
        f.savefig(out + 'trend_timeseries_'+str(m)+'.png')

        plt.close('all')

        #fp = fpath + 'ttrend_WA'+str(m).zfill(2)+'.png'

        #plt.close('all')
