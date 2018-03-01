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


DATE = {'day' : 6,
        'month' : 6,
        'year' : 2006}

def run_waveletDry():


    file = '/users/global/cornkle/data/OBS/modis_LST/modis_netcdf/lsta_daily_'+str(DATE['year'])+str(DATE['month']).zfill(2)+str(DATE['day']).zfill(2)+'.nc'

    ds = xr.open_dataset(file)
    ds = ds.sel(lon=slice(-10,10), lat=slice(10,18))
    lsta = ds['LSTA'][0,:,:]

    # f = plt.figure()
    # plt.imshow(lsta)

    lsta[lsta<-800] = np.nan


    #lsta = lsta   - np.nanmean(lsta)

    points = np.where(np.isfinite(lsta.values))
    inter1 = np.where(np.isnan(lsta.values))

    lsta[inter1] = griddata(points, np.ravel(lsta.values[points]), inter1, method='linear')
    inter = np.where(np.isnan(lsta))
    lsta[inter] = griddata(points, np.ravel(lsta.values[points]), inter, method='nearest')
    #lsta[inter1]=0

    wav = util.waveletLSTA_dom(lsta.values,3)

    wl = wav['dominant']

    wl[inter[0], inter[1]] = np.nan
    wl[inter1[0], inter1[1]] = np.nan
    f = plt.figure()
    plt.imshow(wl, cmap='RdBu', vmin=9, vmax=120)
    scales = wav['scales']

    print(scales)

  #  for id, s in enumerate(scales):
  #      wl[id, :, :][wl[id, :, :] <= s ** .5] = 0

    return wl, scales, lsta, inter1



def wav_checkDry():

    start_time = time.time()
    file = '/users/global/cornkle/MCSfiles/blob_map_30km_-67_JJAS_points.nc'

    wll, scales, lsta, inter = run_waveletDry()

    # day = xr.open_dataarray(dayp)
    # night = xr.open_dataarray(nightp)
    day = xr.open_dataarray(file)

    latmin, latmax = (np.min(lsta['lat']), np.max(lsta['lat']))
    lonmin, lonmax = (np.min(lsta['lon']), np.max(lsta['lon']))

    daystring = str(DATE['year'])+'-'+str(DATE['month']).zfill(2)+'-'+str(DATE['day']).zfill(2)

    day_small = day.sel(time=slice(daystring+' 15:00',daystring+' 19:00'), lat=slice(latmin, latmax), lon=slice(lonmin,lonmax))
    day_small = day_small.where(day_small > -999)
    day_small.values = np.isfinite(day_small.values)
    day_small=day_small.sum(dim='time')
    nightstring = str(DATE['year'])+'-'+str(DATE['month']).zfill(2)+'-'+str(DATE['day']+1).zfill(2)
    night = day.sel(time=slice(nightstring+' 00:00',nightstring+' 05:00'), lat=slice(latmin, latmax), lon=slice(lonmin, lonmax))
    night = night.where(night>-999)
    night.values = np.isfinite(night.values)
    night = night.sum(dim='time')
    wllmean0=wll
    print(scales[0:5])

    bla=day_small.where(day_small>0)
    bla2 = night.where(night>0)


    f= plt.figure()
    f.add_subplot(2,2,1)
    plt.contourf(lsta['lon'], lsta['lat'], lsta, vmin=-8, vmax=8, nlevels=7, cmap='RdBu_r')
    plt.contour(day_small['lon'], day_small['lat'], day_small, cmap='prism')

    plt.contour(lsta['lon'], lsta['lat'], wllmean0, cmap='jet')
    plt.colorbar()

    plt.contourf(day_small['lon'], day_small['lat'], bla, cmap='viridis')


    plt.title(str(pd.to_datetime(lsta['time'].values))+' DRY anomaly DAYTIME')

    # f.add_subplot(2, 2, 2)
    # plt.contourf(lsta['lon'], lsta['lat'], wllmean0)
    # plt.colorbar()
    # plt.title('Sum of wavelet power - all scales')

    f.add_subplot(2, 2, 2)
    plt.contourf(lsta['lon'], lsta['lat'], lsta, vmin=-8, vmax=8, nlevels=7, cmap='RdBu_r')
    plt.colorbar()
    plt.contour(night['lon'], night['lat'], night, cmap='prism')
    plt.contour(lsta['lon'], lsta['lat'], wllmean0, cmap='jet')

    plt.contourf(day_small['lon'], day_small['lat'], bla2, cmap='viridis')

    plt.title(str(pd.to_datetime(lsta['time'].values)) + ' DRY anomaly NIGHTTIME')


    f.add_subplot(2, 2, 3)
    #plt.contourf(day_small['lon'], day_small['lat'], day_small)
    plt.contourf(lsta['lon'], lsta['lat'], wllmean0, cmap='jet')
    plt.colorbar()
    plt.title('Deep conv. frequency, <70C: '+daystring+', Reds: wavelet power')

    f.add_subplot(2, 2, 4)
    #plt.contourf(night['lon'], night['lat'], night)
    plt.contourf(lsta['lon'], lsta['lat'], wllmean0, cmap='jet')
    #plt.gca().invert_yaxis()
    plt.colorbar()
    plt.title('Deep conv. frequency, <70C: '+nightstring+' , Reds: wavelet power')




