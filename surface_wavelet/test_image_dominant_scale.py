import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from wavelet import util
import pdb
from scipy import ndimage
from utils import u_arrays as ua, constants
import pandas as pd
import time
import matplotlib.gridspec as gridspec
from scipy.interpolate import griddata
from utils import u_grid
from utils import u_plot
import matplotlib
import salem as sm

#nice example 25.06.2006

#24.09.2008 and following day wet/dry place makes no sense
#12.08.2009 is cool
#2.8.2009
#29.07.2007

DATE = {'day' : 29,
        'month' : 6,
        'year' : 2006}

def run_waveletDry():


    file = '/users/global/cornkle/data/OBS/MSG_LSTA/lsta_netcdf/lsta_daily_'+str(DATE['year'])+str(DATE['month']).zfill(2)+str(DATE['day']).zfill(2)+'.nc'

    ds = xr.open_dataset(file)
    ds = ds.sel(lon=slice(-10,10), lat=slice(10,18))

    lsta = ds['LSTA'][0,:,:]
    cells = ds['cell'][0,:,:]

    lsta_raw = lsta.copy()

    points = np.where(np.isfinite(lsta.values))
    inter1 = np.where(np.isnan(lsta.values))

    # lsta.values[inter1] = griddata(points, np.ravel((lsta.values)[points]), inter1, method='linear')
    # inter = np.where(np.isnan(lsta))
    # lsta.values[inter] = griddata(points, np.ravel(lsta.values[points]), inter, method='nearest')
    lsta[inter1]=0
    # f = plt.figure()
    # plt.imshow(lsta.values)
    wav = util.waveletLSTA_domLocMax(lsta.values,3)

    wl = wav['dominant']
    power = wav['power']

    power[:, inter1[0], inter1[1]] = np.nan
    wl[inter1[0], inter1[1]] = np.nan

    scales = wav['scales']

    print(scales)

       # for id, s in enumerate(scales):
       #     wl[id, :, :][wl[id, :, :] <= s ** .5] = 0

    return scales, power, lsta, inter1, cells, lsta_raw, wl



def wav_checkDry():

    start_time = time.time()

    # smfile = '/users/global/cornkle/data/OBS/AMSRE/day_aqua/amsre_monthly_anomaly.nc'
    # sm = xr.open_dataset(smfile)
    # sm = sm.sel(lon=slice(-10, 10), lat=slice(10, 18))


    file = '/users/global/cornkle/data/OBS/MSG_LSTA/lsta_netcdf/lsta_daily_'+str(DATE['year'])+str(DATE['month']).zfill(2)+str(DATE['day']).zfill(2)+'.nc'

    ds = xr.open_dataset(file)

    ds = ds.sel(lon=slice(-10,10), lat=slice(10,18))

    lsta2 = ds['LSTA'][0,:,:]


    scales, power, lsta, inter, cells, lsta_raw , wll= run_waveletDry()

    print('Scales:', scales)

    latmin, latmax = (np.min(lsta['lat']), np.max(lsta['lat']))
    lonmin, lonmax = (np.min(lsta['lon']), np.max(lsta['lon']))

    daystring = str(DATE['year'])+'-'+str(DATE['month']).zfill(2)+'-'+str(DATE['day']).zfill(2)

    wllmean0=wll
    ppower_small= np.nansum(power[0:4,:,:], axis=0)
    ppower_big = np.nansum(power[4:7, :, :], axis=0)
    ppower_biggest = np.nansum(power[7::, :, :], axis=0)
    print(scales[0:4])
    print(scales[4:7])
    print(scales[7::])

    xv, yv = np.meshgrid(lsta['lon'], lsta['lat'])
    cmap = u_plot.discrete_cmap(24, base_cmap='gist_ncar')

    era = xr.open_dataset(constants.ERA_DAILY_SRFC_ANO)
    eday = era.sel(longitude=slice(-10, 10), latitude=slice(18, 10),
                   time=str(DATE['year']) + str(DATE['month']).zfill(2) + str(DATE['day']).zfill(2))

    eday['ws10'][0,:,:].plot.pcolormesh() #'p84.162'

    f= plt.figure()

    ax1=f.add_subplot(2,2,1)
    # cmap= plt.get_cmap('RdBu_r')
    # cmap.set_bad(color='gray', alpha=1.)
    # plt.pcolormesh(lsta['lon'], lsta['lat'], lsta_raw, vmin=-8, vmax=8, cmap=cmap)
    #
    # plt.colorbar()
    map = lsta.salem.get_map()
    xl, yl = map.grid.transform(xv, yv)
    map.set_data(lsta_raw)


    map.set_plot_params(cmap='RdBu_r', vmin=-8, vmax=8, extend='both')
    map.visualize(ax=ax1, addcbar=False )
    cax = ax1.scatter(xl, yl, c=cells, cmap=cmap, s=5)
    cbar = plt.colorbar(cax, ticks=np.arange(0,24))
    cbar.set_ticklabels(np.array(np.arange(0,24), dtype=str))


    plt.title(str(pd.to_datetime(lsta['time'].values))+' LSTA')
    #
    # f.add_subplot(2, 2, 2)
    # # plt.contourf(smarr['lon'], smarr['lat'], smarr.values[0,:,:], vmin=-8, vmax=8, nlevels=7, cmap='RdBu_r')
    # # plt.colorbar()
    # # plt.contour(night['lon'], night['lat'], night, cmap='prism')
    # arr = np.zeros_like(lsta2.values)
    # arr[lsta2.values<-0] = -4
    # arr[lsta2.values>0] = 4
    #
    # plt.pcolormesh(lsta['lon'], lsta['lat'], arr, vmin=-8, vmax=8, cmap='RdBu_r')
    # plt.colorbar()
    # plt.scatter(xv, yv, c=cells, cmap=cmap, s=5)
    # #
    # # plt.contourf(day_small['lon'], day_small['lat'], night, cmap='viridis')
    #
    # plt.title(str(pd.to_datetime(lsta['time'].values)) + ' Land surface temperature anomaly, +/-4')


    f.add_subplot(2, 2, 3)
    #plt.contourf(day_small['lon'], day_small['lat'], day_small)
    plt.pcolormesh(lsta['lon'], lsta['lat'], ppower_small, cmap='RdBu_r', vmin=-10, vmax=10)
    plt.colorbar()
    plt.scatter(xv, yv, c=cells, cmap=cmap, s=5)
    #plt.scatter(xv, yv, c=cells, cmap='jet', s=5)
    plt.title('Wavelet powers < 42km')


    f.add_subplot(2, 2, 4)
    #plt.contourf(night['lon'], night['lat'], night)
    plt.pcolormesh(lsta['lon'], lsta['lat'], ppower_big, cmap='RdBu_r', vmin=-10, vmax=10)
    plt.colorbar()
    plt.scatter(xv, yv, c=cells, cmap=cmap, s=5)
    #plt.scatter(xv, yv, c=cells, cmap='jet', s=5)
    #plt.gca().invert_yaxis()

    plt.title('Wavelet powers  35-51km')

    f.add_subplot(2, 2, 2)
    # plt.contourf(night['lon'], night['lat'], night)
    plt.pcolormesh(lsta['lon'], lsta['lat'], ppower_biggest, cmap='RdBu_r', vmin=-10, vmax=10)
    plt.colorbar()
    plt.scatter(xv, yv, c=cells, cmap=cmap, s=5)
    # plt.scatter(xv, yv, c=cells, cmap='jet', s=5)
    # plt.gca().invert_yaxis()

    plt.title('Wavelet powers  35-51km')




