import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from wavelet import util
import pdb
from scipy import ndimage
from utils import u_arrays as ua
import pandas as pd
import time
import matplotlib.gridspec as gridspec
from scipy.interpolate import griddata
from utils import u_grid
import matplotlib

#nice example 25.06.2006

#24.09.2008 and following day wet/dry place makes no sense
#12.08.2009 is cool
#2.8.2009
#29.07.2007

DATE = {'day' : 12,
        'month' : 8,
        'year' : 2009}

def run_waveletDry():


    file = '/users/global/cornkle/data/OBS/MSG_LSTA/lsta_new/lsta_daily_'+str(DATE['year'])+str(DATE['month']).zfill(2)+str(DATE['day']).zfill(2)+'.nc'
    file2 = '/users/global/cornkle/data/OBS/MSG_LSTA/lsta_netcdf/lsta_daily_' + str(DATE['year']) + str(
        DATE['month']).zfill(2) + str(DATE['day']).zfill(2) + '.nc'

    ds = xr.open_dataset(file)
    ds2 = xr.open_dataset(file2)
    ds = ds.sel(lon=slice(-10,10), lat=slice(10,18))
    ds2 = ds2.sel(lon=slice(-10, 10), lat=slice(10, 18))

    lsta = ds['LSTA'][0,:,:]
    cells = ds2['cell'][0,:,:]

    # f = plt.figure()
    # plt.imshow(lsta)

    lsta = lsta.where(lsta>-800)
    lsta_raw = lsta.copy()
    #lsta_raw.values[lsta_raw.values>0] =10
    #lsta_raw.values[lsta_raw.values < 0] = -10


    #lsta = lsta - lsta.mean()set

    # f = plt.figure()
    # plt.imshow(lsta.values)

    points = np.where(np.isfinite(lsta.values))
    inter1 = np.where(np.isnan(lsta.values))


    lsta.values[inter1] = griddata(points, np.ravel((lsta.values)[points]), inter1, method='linear')
    inter = np.where(np.isnan(lsta))
    lsta.values[inter] = griddata(points, np.ravel(lsta.values[points]), inter, method='nearest')
    #lsta[inter1]=0
    # f = plt.figure()
    # plt.imshow(lsta.values)
    wav = util.LSTA_LocalMax(lsta.values, 3)

    wl = wav['dominant']
    power = wav['power']

    # wl[inter[0], inter[1]] = np.nan
    # wl[inter1[0], inter1[1]] = np.nan
    # f = plt.figure()
    # plt.imshow(wl, cmap='RdBu', vmin=9, vmax=120)
    scales = wav['scales']

    print(scales)

  #  for id, s in enumerate(scales):
  #      wl[id, :, :][wl[id, :, :] <= s ** .5] = 0

    return wl, scales, power, lsta, inter1, cells, lsta_raw



def wav_checkDry():

    start_time = time.time()

    # smfile = '/users/global/cornkle/data/OBS/AMSRE/day_aqua/amsre_monthly_anomaly.nc'
    # sm = xr.open_dataset(smfile)
    # sm = sm.sel(lon=slice(-10, 10), lat=slice(10, 18))


    file = '/users/global/cornkle/data/OBS/MSG_LSTA/lsta_netcdf/lsta_daily_'+str(DATE['year'])+str(DATE['month']).zfill(2)+str(DATE['day']).zfill(2)+'.nc'

    ds = xr.open_dataset(file)

    ds = ds.sel(lon=slice(-10,10), lat=slice(10,18))

    lsta2 = ds['LSTA'][0,:,:]


    wll, scales, power, lsta, inter, cells, lsta_raw = run_waveletDry()

    # day = xr.open_dataarray(dayp)
    # night = xr.open_dataarray(nightp)
    #cells[np.isnan(cells)]=-1

    latmin, latmax = (np.min(lsta['lat']), np.max(lsta['lat']))
    lonmin, lonmax = (np.min(lsta['lon']), np.max(lsta['lon']))

    daystring = str(DATE['year'])+'-'+str(DATE['month']).zfill(2)+'-'+str(DATE['day']).zfill(2)
    #smarr = sm['SM'].sel(time=daystring)

    wllmean0=wll
    ppower_small= np.nansum(power[0:8,:,:], axis=0)
    ppower_big = np.nansum(power[-8::, :, :], axis=0)
    print(scales[0:8])

    xv, yv = np.meshgrid(lsta['lon'], lsta['lat'])

    f= plt.figure()
    f.add_subplot(2,2,1)
    cmap= plt.get_cmap('RdBu_r')
    cmap.set_bad(color='gray', alpha=1.)
    plt.pcolormesh(lsta['lon'], lsta['lat'], lsta_raw, vmin=-8, vmax=8, cmap=cmap)
    plt.scatter(xv, yv, c=cells, cmap='jet', s=5)
    plt.colorbar()
    # plt.contour(lsta['lon'], lsta['lat'], wllmean0, cmap='jet')
    # plt.colorbar()

    #plt.plot(day_pos[1], day_pos[0], 'ro')


    plt.title(str(pd.to_datetime(lsta['time'].values))+' DRY anomaly DAYTIME')

    f.add_subplot(2, 2, 2)
    # plt.contourf(smarr['lon'], smarr['lat'], smarr.values[0,:,:], vmin=-8, vmax=8, nlevels=7, cmap='RdBu_r')
    # plt.colorbar()
    # plt.contour(night['lon'], night['lat'], night, cmap='prism')
    plt.pcolormesh(lsta['lon'], lsta['lat'], lsta2, vmin=-8, vmax=8, cmap=cmap)
    plt.scatter(xv, yv, c=cells, cmap='jet', s=5)
    #
    # plt.contourf(day_small['lon'], day_small['lat'], night, cmap='viridis')

    plt.title(str(pd.to_datetime(lsta['time'].values)) + ' LSTA following day')
    plt.colorbar()

    f.add_subplot(2, 2, 3)
    #plt.contourf(day_small['lon'], day_small['lat'], day_small)
    plt.pcolormesh(lsta['lon'], lsta['lat'], ppower_small, cmap='Reds', vmax=30)
    plt.colorbar()
    plt.scatter(xv, yv, c=cells, cmap='jet', s=5)
    plt.title('Deep conv. frequency, <50C:, Reds: wavelet power')

    f.add_subplot(2, 2, 4)
    #plt.contourf(night['lon'], night['lat'], night)
    plt.pcolormesh(lsta['lon'], lsta['lat'], wllmean0, cmap='RdBu_r', vmin=-90, vmax=90)
    plt.colorbar()
    plt.scatter(xv, yv, c=cells, cmap='jet', s=5)
    #plt.gca().invert_yaxis()

    plt.title('Deep conv. frequency, <50C:, Reds: wavelet power')




