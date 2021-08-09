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

    y1 = 2000
    y2 = 2020

    if method is None:
        'Please provide trend calc method: polyfit or mk (mann kendall)'

    data = data[(data['year'] >= y1) & (data['year'] <= y2)]


    mean_years = data.mean('year')

    # stack lat and lon into a single dimension called allpoints
    datastacked = data.stack(allpoints=['lat', 'lon'])


    # apply the function over allpoints to calculate the trend at each point
    print('Entering trend calc')

    alpha = 0.05

    if method=='polyfit':
        dtrend = datastacked.groupby('allpoints').apply(u_darrays.linear_trend_lingress,nb_valid=5)
        dtrend = dtrend.unstack('allpoints')

        if sig:
            (dtrend['slope'].values)[dtrend['pval'].values > alpha] = 0

    ddtrend = dtrend['slope']

    #ipdb.set_trace()
    # unstack back to lat lon coordinates
    #ipdb.set_trace()
    ddout = ddtrend #uda.flip_lat(ddtrend)

    return ddout, mean_years



def trend_all():



    mcs = '/media/ck/Elements/global/GPM/*.nc4' #'/media/ck/Elements/SouthAmerica/CHIRPS/chirps-v2.0.monthly.nc'
    mcs = '/media/ck/Elements/global/CHIRPS/chirps-v2.0.monthly.nc'
    fpath = cnst.network_data + 'figs/CLOVER/CHIRPS/'


    da3_out = xr.open_dataset(mcs)#/100
    da3_out = da3_out.rename({'longitude':'lon', 'latitude':'lat'})


    months = [311]



    f = plt.figure(figsize=(9, 7), dpi=300)


    for ids, m in enumerate(months):
       # ipdb.set_trace()
        tir = da3_out.sel(time=(da3_out['time.month']>=3) &(da3_out['time.month']<=11))#.load()
        box = [-13, -6, 4, 10] #[-13, 10, 4, 10] #

        tir = tir['precip'].sel(lon=slice(box[0], box[1]), lat=slice(box[2], box[3])).groupby('time.year').mean('time')/30/24
        #ipdb.set_trace()
        print('Doing', m)
        method = 'polyfit'

        if type(m)==int:
            m = [m]

        sig = False


        tirtrend, tirmean = calc_trend(tir, m, method=method, sig=sig, wilks=False)

        tirm_mean = tirmean#.mean('time')


        tirtrend_unstacked = ((tirtrend)*20.)#/ tirm_mean.values) * 100. #/ tirm_mean.values

        ti_da = tirtrend_unstacked


        if len(m) == 1:
            fp = fpath + 'CHIRPS_only_trendmap_Allmonths_mediumDomain_1983-2018_'+str(m[0]).zfill(2)+'.png'
        else:
            fp = fpath + 'CHIRPS_only_trendmap_' + str(m[0]).zfill(2) +'-'+ str(m[1]).zfill(2) + '.png'

        map = tir.salem.get_map()

        ax1 = f.add_subplot(2,2,ids+1)
      #  ipdb.set_trace()
        map.set_contour((tirm_mean), interp='linear', levels=np.array([1,5, 10, 20, 40,80,100])/100, colors='k', linewidths=0.5)

        coord = [-13, -10, 7.3, 7.8]
        geom = shpg.box(coord[0], coord[2], coord[1], coord[3])
        map.set_geometry(geom, zorder=99, color='darkorange', linewidth=2, linestyle='--', alpha=0.3)

        coord = [-10.4, -5.5, 6.2, 6.8]
        geom2 = shpg.box(coord[0], coord[2], coord[1], coord[3])
        map.set_geometry(geom2, zorder=99, color='darkorange', linewidth=2, linestyle='--', alpha=0.3)

        ti_da.values[ti_da.values==0] = np.nan

        map.set_data(ti_da)
        # np.linspace(0, 8, 6)/100)
        map.set_plot_params(cmap='RdBu', extend='both', levels = np.linspace(-8, 8, 12)/100) #)  #, levels=np.arange(-7,7,25)  # levels=np.arange(20,101,20)  #np.arange(20,101,20)
        #map.set_plot_params(cmap='RdBu', extend='both', levels=[-100, -80,-60,-40,-20,-10,-5,-2,0,2,5,10,20,40,60,80,100])
        dic = map.visualize(ax=ax1, title='CHIRPS Mar-Nov 2000-2020 precip trend', cbar_title='mm/h',addcbar=True)
        contours = dic['contour'][0]
        plt.clabel(contours, inline=True, fontsize=7, fmt='%1.1f')

    ax1 = f.add_subplot(2, 2, 3)
    plt.plot(tirm_mean.sel(lon=slice(-13,-10)).lon, tirm_mean.sel(lon=slice(-13,-10), lat=slice(7.3,7.8)).mean('lat'), color='r', label='mean')
    plt.plot(tirm_mean.sel(lon=slice(-13,-10)).lon, ti_da.sel(lon=slice(-13,-10), lat=slice(7.3,7.8)).mean('lat'), label='trend')
    plt.title('Sierra Leone box')
    plt.axvline(-11.9, color='k')
    plt.axvline(-11.2, color='k')
    plt.axhline(0, color='k')

    plt.legend()

    ax1 = f.add_subplot(2, 2, 4)
    plt.plot(tirm_mean.sel(lon=slice(-10.4,-5.5)).lon, tirm_mean.sel(lon=slice(-10.4,-5.5), lat=slice(6.2,6.8)).mean('lat'), color='r', label='mean')
    plt.plot(tirm_mean.sel(lon=slice(-10.4,-5.5)).lon, ti_da.sel(lon=slice(-10.4,-5.5), lat=slice(6.2,6.8)).mean('lat'), label='trend')
    plt.axvline(-7.5, color='k')
    plt.axvline(-8.4, color='k')
    plt.axhline(0, color='k')
    plt.title("Cote d'Ivoire box")
    plt.legend()

    plt.tight_layout()
    plt.savefig(fp)
    plt.close('all')
