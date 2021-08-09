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


def calc_trend(data, method=None, sig=False):

    y1 = 2000
    y2 = 2020

    if method is None:
        'Please provide trend calc method: polyfit or mk (mann kendall)'
    #ipdb.set_trace()
    data = data.isel(year=(data['year'] >= y1) & (data['year'] <= y2))

    data = data.where(data['precipitationQualityIndex']>5).mean(['lat','lon'])

    mean_years = data['precipitation'].mean(['year'])

    alpha = 0.05

    if method=='polyfit':
        dtrend = u_darrays.linear_trend_lingress(data['precipitation'],nb_valid=5)
        if sig:
            (dtrend['slope'].values)[dtrend['pval'].values > alpha] = 0

    ddout = dtrend #uda.flip_lat(ddtrend)

    f = plt.figure()
    ax = f.add_subplot(121)
    ax.plot(data.year, data['precipitation'], label=str(mean_years))
    plt.legend()


    return ddout, mean_years



def trend_all():


    mcs = '/media/ck/Elements/global/GPM/*.nc4' #'/media/ck/Elements/SouthAmerica/CHIRPS/chirps-v2.0.monthly.nc'


    da3_out = xr.open_mfdataset(mcs)


    tir = da3_out.sel(time=(da3_out['time.month']>=3) &(da3_out['time.month']<=11)).load()
    box = [-13, 10, 4, 10] #[-13, 10, 4, 10] #
    #tir = tir.sel(lon=slice(box[0], box[1]), lat=slice(box[2], box[3])).groupby('time.year').mean('time')
    tir = tir.groupby('time.year').mean('time')

    method = 'polyfit'

    sig = False


    tirtrend, tirmean = calc_trend(tir, method=method, sig=sig)
    #ipdb.set_trace()
    print('Trend', tirtrend['slope'].values, 'pval', tirtrend['pval'].values)
    print('Mean', tirmean.values)
    print('Relative', tirtrend['slope'].values/tirmean.values*100)
