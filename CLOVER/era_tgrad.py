import pandas as pd
import numpy as np
import xarray as xr
import pdb
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy
import salem
from utils import u_plot as up, u_mann_kendall as mk, constants as cnst
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



def t_trend_polyfit_year():
    #file = '/users/global/cornkle/data/ERA-I monthly/ERA-WA-Monthly-2mTemp.nc'
    file = cnst.ERA_MONTHLY_SRFC_SYNOP

    fpath = cnst.network_data + '/figs/CLOVER/months/'

    # define a function to compute a linear trend of a timeseries
    def linear_trend(x):
        # pf = np.polyfit(np.arange(len(x)), x, 1)
        pf, slope, int, p, ind = mk.test(np.arange(len(x)), x.squeeze().values, eps=0.001, alpha=0.01, Ha='upordown')

        # we need to return a dataarray or else xarray's groupby won't be happy
        issig = slope

        # if ind == 1:
        #     issig = slope
        # else:
        #     issig = np.nan

        return xr.DataArray(issig, )


    ## a clean way of plotting - use matplotlib functions directly:

    def draw_map(ax, data, lon, lat, title=None, mask_sig=None, quiver=None, contour=None, cbar_label=None, toggle_cbar=None, **kwargs):

        mapp = ax.contourf(lon, lat, data, transform=ccrs.PlateCarree(), **kwargs)  # this is the actual plot

        ## mask for significance indicator
        if mask_sig is not None:
            plt.contourf(lon, lat, mask_sig, colors='none', hatches='.',
                         levels=[0.5, 1], linewidth=0.1)

        ## quiver list
        if quiver is not None:
            qu = ax.quiver(quiver['x'], quiver['y'], quiver['u'], quiver['v'], scale=quiver['scale'])
        ## additional contour on plot
        if contour is not None:
            ax.contour(contour['x'], contour['y'], contour['data'], levels=contour['levels'], cmap=contour['cmap'])

        ax.coastlines()  ## adds coastlines
        # Gridlines
        xl = ax.gridlines(draw_labels=True);  # adds latlon grid lines
        xl.xlabels_top = False  ## labels off
        xl.ylabels_right = False
        plt.title(title)
        # Countries
        ax.add_feature(cartopy.feature.BORDERS, linestyle='--');  # adds country borders
        if toggle_cbar is not None:
            cbar = plt.colorbar(mapp)  # adds colorbar
            cbar.set_label(cbar_label)
        plt.show()


    dam = xr.open_dataset(file)
    dam = dam['t2m']

    f = plt.figure(figsize=(9,11))

    m_dic = { 1 : 'Jan', 2 : 'Feb', 3 : 'Mar', 4 : 'Apr', 5 : 'May', 6 : 'Jun', 7 : 'Jul', 8 : 'Aug',
              9 : 'Sep', 10 : 'Oct', 11 : 'Nov', 12 : '12'}

    for mm in range(1,13):

        ax = f.add_subplot(4,3,mm, projection=ccrs.PlateCarree())  # this opens a new plot axis

        da = dam[(dam['time.month']==mm)]

        str_tit = m_dic[mm]

        da = da.sel(longitude=slice(-18,51), latitude=slice(36, -37))
        da = da.groupby('time.year').mean(axis=0)

        lons = da.longitude
        lats = np.flip(da.latitude.values, axis=0)

        # stack lat and lon into a single dimension called allpoints
        stacked = da.stack(allpoints=['latitude','longitude'])
        # apply the function over allpoints to calculate the trend at each point
        trend = stacked.groupby('allpoints').apply(linear_trend)
        # unstack back to lat lon coordinates
        trend_unstacked = trend.unstack('allpoints')

        trend_unstacked = trend_unstacked*10. # warming over decade
        da2 = xr.DataArray(trend_unstacked, coords=[lats, lons], dims=['latitude', 'longitude'])

        bla = True
        # if mm in [3,6,9,12]:
        #     bla = True

        draw_map(ax, da2.values, da2.longitude, da2.latitude, levels=[-0.4, -0.3, -0.2, 0.2, 0.3, 0.4],
                 extend='both', title=str_tit, mask_sig=None, quiver=None, contour=None, cbar_label='K decade$^{-1}$', cmap='RdBu_r', toggle_cbar=bla)


        print('Doing ', str(mm))

    plt.tight_layout()
    fp = fpath + 'ttrend_year.png'

    f.savefig(fp)
    plt.close('all')


def t_trend_polyfit_season():
    #file = '/users/global/cornkle/data/ERA-I monthly/ERA-WA-Monthly-2mTemp.nc'
    file = cnst.ERA_MONTHLY_SRFC_SYNOP

    fpath = cnst.network_data + '/figs/CLOVER/months/ERA5_WA/'

    # define a function to compute a linear trend of a timeseries
    def linear_trend(x):
        # pf = np.polyfit(np.arange(len(x)), x, 1)
        pf, slope, int, p, ind = mk.test(np.arange(len(x)), x.squeeze().values, eps=0.001, alpha=0.01, Ha='upordown')

        # we need to return a dataarray or else xarray's groupby won't be happy
        issig = slope
        #
        # if ind == 1:
        #     issig = slope
        # else:
        #     issig = np.nan

        return xr.DataArray(issig, )


    ## a clean way of plotting - use matplotlib functions directly:

    def draw_map(ax, data, lon, lat, title=None, mask_sig=None, quiver=None, contour=None, cbar_label=None, toggle_cbar=None, **kwargs):

        mapp = ax.contourf(lon, lat, data, transform=ccrs.PlateCarree(), **kwargs)  # this is the actual plot

        ## mask for significance indicator
        if mask_sig is not None:
            plt.contourf(lon, lat, mask_sig, colors='none', hatches='.',
                         levels=[0.5, 1], linewidth=0.1)

        ## quiver list
        if quiver is not None:
            qu = ax.quiver(quiver['x'], quiver['y'], quiver['u'], quiver['v'], scale=quiver['scale'])
        ## additional contour on plot
        if contour is not None:
            ax.contour(contour['x'], contour['y'], contour['data'], levels=contour['levels'], cmap=contour['cmap'])

        ax.coastlines()  ## adds coastlines
        # Gridlines
        xl = ax.gridlines(draw_labels=True);  # adds latlon grid lines
        xl.xlabels_top = False  ## labels off
        xl.ylabels_right = False
        plt.title(title)
        # Countries
        ax.add_feature(cartopy.feature.BORDERS, linestyle='--');  # adds country borders
        if toggle_cbar is not None:
            cbar = plt.colorbar(mapp)  # adds colorbar
            cbar.set_label(cbar_label)
        plt.show()


    dam = xr.open_dataset(file)
    dam = dam['t2m']

    f = plt.figure(figsize=(10,9), dpi=300)

    m_dic = { 1 : 'Jan', 2 : 'Feb', 3 : 'Mar', 4 : 'Apr', 5 : 'May', 6 : 'Jun', 7 : 'Jul', 8 : 'Aug',
              9 : 'Sep', 10 : 'Oct', 11 : 'Nov', 12 : '12'}

    seasons = [ [12,2, 'DJF',1], [3,5, 'MAM',2], [6,8, 'JJA',3],[9,11, 'SON',4]]

    for mm in seasons:

        ax = f.add_subplot(2,2,mm[3], projection=ccrs.PlateCarree())  # this opens a new plot axis
        if mm[0] == 12:
            da = dam[(dam['time.month']>=mm[0]) | (dam['time.month']<=mm[1])]
        else:
            da = dam[(dam['time.month'] >= mm[0]) & (dam['time.month'] <= mm[1])]

        str_tit = mm[2]

        da = da.sel(longitude=slice(-18,51), latitude=slice(36, -37))
        da = da.groupby('time.year').mean(axis=0)

        lons = da.longitude
        lats = np.flip(da.latitude.values, axis=0)

        # stack lat and lon into a single dimension called allpoints
        stacked = da.stack(allpoints=['latitude','longitude'])
        # apply the function over allpoints to calculate the trend at each point
        #u_darrays.linear_trend_mk, alpha=alpha, eps=0.0001,nb_missing=10
        trend = stacked.groupby('allpoints').apply(linear_trend)
        # unstack back to lat lon coordinates
        trend_unstacked = trend.unstack('allpoints')

        trend_unstacked = trend_unstacked*10 # warming over decade
        da2 = xr.DataArray(trend_unstacked, coords=[lats, lons], dims=['latitude', 'longitude'])

        bla = True
        # if mm in [3,6,9,12]:
        #     bla = True

        draw_map(ax, da2.values, da2.longitude, da2.latitude, levels=[-0.4, -0.3, -0.2, 0.2, 0.3, 0.4],
                 extend='both', title=str_tit, mask_sig=None, quiver=None, contour=None, cbar_label='K decade$^{-1}$', cmap='RdBu_r', toggle_cbar=bla)


        print('Doing ', str(mm))

    plt.tight_layout()
    fp = fpath + 'ttrend_year.png'

    f.savefig(fp)
    plt.close('all')

def t_trend_polyfit_month():
    #file = '/users/global/cornkle/data/ERA-I monthly/ERA-WA-Monthly-2mTemp.nc'
    file = cnst.ERA_MONTHLY_SRFC_SYNOP

    fpath = cnst.network_data + '/figs/CLOVER/months/'

    # define a function to compute a linear trend of a timeseries
    def linear_trend(x):
        # pf = np.polyfit(np.arange(len(x)), x, 1)
        pf, slope, int, p, ind = mk.test(np.arange(len(x)), x.squeeze().values, eps=0.001, alpha=0.01, Ha='upordown')

        # we need to return a dataarray or else xarray's groupby won't be happy
        issig = slope

        # if ind == 1:
        #     issig = slope
        # else:
        #     issig = np.nan

        return xr.DataArray(issig, )


    ## a clean way of plotting - use matplotlib functions directly:

    def draw_map(ax, data, lon, lat, title=None, mask_sig=None, quiver=None, contour=None, cbar_label=None, toggle_cbar=None, **kwargs):

        mapp = ax.contourf(lon, lat, data, transform=ccrs.PlateCarree(), **kwargs)  # this is the actual plot

        ## mask for significance indicator
        if mask_sig is not None:
            plt.contourf(lon, lat, mask_sig, colors='none', hatches='.',
                         levels=[0.5, 1], linewidth=0.1)

        ## quiver list
        if quiver is not None:
            qu = ax.quiver(quiver['x'], quiver['y'], quiver['u'], quiver['v'], scale=quiver['scale'])
        ## additional contour on plot
        if contour is not None:
            ax.contour(contour['x'], contour['y'], contour['data'], levels=contour['levels'], cmap=contour['cmap'])

        ax.coastlines()  ## adds coastlines
        # Gridlines
        xl = ax.gridlines(draw_labels=True);  # adds latlon grid lines
        xl.xlabels_top = False  ## labels off
        xl.ylabels_right = False
        plt.title(title)
        # Countries
        ax.add_feature(cartopy.feature.BORDERS, linestyle='--');  # adds country borders
        if toggle_cbar is not None:
            cbar = plt.colorbar(mapp)  # adds colorbar
            cbar.set_label(cbar_label)
        plt.show()


    dam = xr.open_dataset(file)
    dam = dam['t2m']

    f = plt.figure(figsize=(9,8), dpi=300)

    m_dic = { 1 : 'Jan', 2 : 'Feb', 3 : 'Mar', 4 : 'Apr', 5 : 'May', 6 : 'Jun', 7 : 'Jul', 8 : 'Aug',
              9 : 'Sep', 10 : 'Oct', 11 : 'Nov', 12 : '12'}

    seasons = [ [11,2, 'Nov',1], [2,5, 'Feb',2], [5,8, 'May',3],[8,11, 'Aug',4]]

    for mm in seasons:

        ax = f.add_subplot(2,2,mm[3], projection=ccrs.PlateCarree())  # this opens a new plot axis

        da = dam[(dam['time.month'] == mm[0])]

        str_tit = mm[2]

        da = da.sel(longitude=slice(-18,51), latitude=slice(36, -37))
        da = da.groupby('time.year').mean(axis=0)

        lons = da.longitude
        lats = np.flip(da.latitude.values, axis=0)

        # stack lat and lon into a single dimension called allpoints
        stacked = da.stack(allpoints=['latitude','longitude'])
        # apply the function over allpoints to calculate the trend at each point
        trend = stacked.groupby('allpoints').apply(linear_trend)
        # unstack back to lat lon coordinates
        trend_unstacked = trend.unstack('allpoints')

        trend_unstacked = trend_unstacked*10. # warming over decade
        da2 = xr.DataArray(trend_unstacked, coords=[lats, lons], dims=['latitude', 'longitude'])

        bla = True
        # if mm in [3,6,9,12]:
        #     bla = True

        draw_map(ax, da2.values, da2.longitude, da2.latitude, levels=[-0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4],
                 extend='both', title=str_tit, mask_sig=None, quiver=None, contour=None, cbar_label='K decade$^{-1}$', cmap='RdBu_r', toggle_cbar=bla)


        print('Doing ', str(mm))

    plt.tight_layout()
    fp = fpath + 'ttrend_year_singleM.png'

    f.savefig(fp)
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
    mcs = cnst.GRIDSAT + 'aggs/box_13W-13E-4-8N_meanT-50_from5000km2.nc' #box_25-33E-28-10S_meanT-50_from5000km2.nc'#box_13W-13E-4-8N_meanT-50_from5000km2.nc'
    out = cnst.network_data + 'figs/CLOVER/timeseries/'

    #
    box = [-10,10,5.5,8]
    tpick = [-10,10,6,25]
    Tlons = [-10,10]

    # box = [18,30,-22,-10]
    # tpick = [18,30,-30,-20]
    # Tlons = [18,30]

    dam = xr.open_dataset(srfc)
    dam = u_darrays.flip_lat(dam)
    dam = dam['t2m']

    tgrad = dam.sel(longitude=slice(tpick[0], tpick[1]), latitude=slice(tpick[2],tpick[3])).mean('longitude').groupby('time.month').mean('time')

    Tgrad_lat = []

    for tgrad_ts in tgrad:

        ttgrad = np.argmax(tgrad_ts.values[2::] - tgrad_ts.values[0:-2])

        lat_pos = tgrad.latitude[ttgrad+1]

        Tgrad_lat.append(float(lat_pos))

    da = xr.open_dataset(pl)
    da = u_darrays.flip_lat(da)
    da = da.sel(longitude=slice(box[0], box[1]), latitude=slice(box[2],box[3]))#latitude=slice(36, -37))

    u925 = da['u'].sel(level=925)#(level=slice(850, 925)).mean(dim='level')#).mean(dim='level')  #slice(850
    u600 = da['u'].sel(level=700)#(level=slice(550, 650)).mean(dim='level')

    qq925 = da['q'].sel(level=slice(850,925)).mean(dim='level')
    qq600 = da['q'].sel(level=slice(550, 650)).mean(dim='level')

    mcs_temp = xr.open_dataset(mcs)
    mcs_temp = mcs_temp['tir']
    months = [1,2,3,4,5,6,7,8,9,10,11,12]
    mnames = {1:'January',2:'February', 3 : 'March', 4 : 'April', 5 : 'May', 6:'June', 7:'July', 8:'August', 9 : 'September', 10 : 'October', 11:'November', 12:'December'}

    mcs_change = []
    shear_change = []
    tgrad_change = []

    for m in months:

        tsouthlice = slice(Tgrad_lat[m-1]-1.25, Tgrad_lat[m-1]-0.25)
        tnorthslice = slice(Tgrad_lat[m-1]+0.25, Tgrad_lat[m-1]+1.25)

        tsouth = dam.sel(longitude=slice(Tlons[0], Tlons[1]), latitude=slice(5,8))
        tnorth = dam.sel(longitude=slice(Tlons[0], Tlons[1]), latitude=slice(11,18))

        south = tsouth[(tsouth['time.month']==m)]
        north = tnorth[(tnorth['time.month']==m)]
        ulow = u925[(u925['time.month']==m)]
        uhigh = u600[(u600['time.month']==m)]
        mcs_month = mcs_temp[mcs_temp['time.month']==m]#*(-1)
        qlow = qq925[(qq925['time.month']==m)]
        qmid = qq600[(qq600['time.month'] == m)]
        ##da = da.sel(longitude=slice(-18,51), latitude=slice(36, -37))

        south_peryear = south.groupby('time.year').mean('longitude').min('latitude')
        north_peryear = north.groupby('time.year').mean('longitude').max('latitude')

        u925_peryear = ulow.groupby('time.year').mean('longitude').max('latitude') #ulow.groupby('time.year').mean()

        u600_peryear = uhigh.groupby('time.year').mean('longitude').min('latitude')#.mean() # ('latitude').min()

        qlow_peryear = qlow.groupby('time.year').mean('longitude').max('latitude')#.mean() # ('latitude').min()
        qmid_peryear = qmid.groupby('time.year').mean('longitude').max('latitude')


        tgrad = ((north_peryear-south_peryear)[4::])
        shear = (u600_peryear-u925_peryear)[4::] # -q_peryear[4::]#
        q = qlow_peryear[4::]*1000

        r = stats.pearsonr(shear.values[1::]-shear.values[0:-1],mcs_month.values[1::]-mcs_month.values[0:-1])
        rq = stats.pearsonr(q.values[1::]-q.values[0:-1],mcs_month.values[1::]-mcs_month.values[0:-1])
        tshear_cor = stats.pearsonr(shear.values[1::]-shear.values[0:-1],tgrad.values[1::]-tgrad.values[0:-1])

        #sslope, sint = ustats.linear_trend(shear)
        sslope, sint = ustats.linear_trend(shear)
        sslope, sint, sr, sp, sstd_err = stats.linregress(np.arange(len(shear)), shear)
        try:
            #mslope, mint = ustats.linear_trend(mcs_month)
            mslope, mint, mrr, mp, mstd_err = stats.linregress(np.arange(len(mcs_month)), mcs_month)

        except:
            continue

        mr = stats.pearsonr(np.arange(len(mcs_month)), mcs_month)
        print(m, mr, mrr, mslope)
        sig = 'SIG'
        if mr[1]>0.05:
            #continue
            sig = 'not_sig'

        mcs_change.append(mcs_month[-5::].mean()-mcs_month[0:5].mean())
        shear_change.append(shear[-5::].mean()-shear[0:5].mean())
        tgrad_change.append(tgrad[-5::].mean()-tgrad[0:5].mean())

        x = np.arange(0,len(shear))
        rr = r[0]
        f=plt.figure(figsize=(6,3))
        ax = f.add_subplot(111)
        ax.plot(tgrad.year, shear, 'x-', label='Zonal wind shear 600-925hPa', color='k')
        ax.plot(tgrad.year, sint + x*sslope, '--', color='k')
        ax.text(0.5,0.8, 'MCS: '+ sig +'_' + str(np.round(mslope,decimals=2)*10) + ' p-val: ' + str(np.round(mp, decimals=2)), transform=ax.transAxes)
        ax.set_ylim(-8,0)
        #ax.set_ylim(9.8, 10)
        ax1 = ax.twinx()
        ax1.plot(mcs_month['time.year'],mcs_month, 'o-', label='Mean MCS temp.', color='r')
        ax1.plot(tgrad.year, mint + x*mslope, '--', color='r')
        mcsline = lines.Line2D([],[], color='r', label='Mean MCS temp.', linestyle='solid', marker='o')
        shearline = lines.Line2D([],[], color='k', label='Zonal wind shear 600-925hPa', linestyle='solid', marker='x', markersize=5)
        ax1.set_ylabel('degC')
        ax.set_ylabel('m s-1')
        ax.set_title(mnames[m]+' | Corr.:'+ str(np.round(rr, decimals=2)) +'|' + str(np.round(rq[0], decimals=2)) + '| Tgrad/Shear corr: ' + str(np.round(tshear_cor[0], decimals=2)))
        if m==3:
            ax.legend(handles=[mcsline, shearline], loc='lower left')
        f.savefig(out + 'timeseries_WA'+str(m)+'.png')

        plt.close('all')

        #fp = fpath + 'ttrend_WA'+str(m).zfill(2)+'.png'

        #plt.close('all')

    # plt.figure()
    # plt.plot([2,3,4,5,6,7,8,9,10,11], tgrad_change, 'ro-', label='tgrad')
    # plt.plot([2,3,4,5,6,7,8,9,10,11], shear_change, 'ko-', label='u600')
    # plt.plot([2,3,4,5,6,7,8,9,10,11], mcs_change, 'bo-', label='mcs mean temp')
    # plt.legend()
