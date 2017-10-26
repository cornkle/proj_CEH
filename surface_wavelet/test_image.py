import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from wavelet import util
import pdb
from scipy import ndimage
from utils import u_arrays as ua
import pandas as pd


DATE = {'day' : 2,
        'month' : 8,
        'year' : 2009}

def run_waveletDry():


    file = '/users/global/cornkle/data/OBS/modis_LST/modis_netcdf/lsta_daily_'+str(DATE['year'])+str(DATE['month']).zfill(2)+str(DATE['day']).zfill(2)+'.nc'

    ds = xr.open_dataset(file)

    lsta = ds['LSTA'][0,:,:]
    lsta[lsta<-800] = np.nan
    lsta = lsta - np.nanmean(lsta)

    wav = util.waveletLSTA(lsta.values,3, dry=True)

    wl = wav['power']
    scales = wav['scales']
    print(scales)

  #  for id, s in enumerate(scales):
  #      wl[id, :, :][wl[id, :, :] <= s ** .5] = 0

    return wl, scales, lsta

def run_waveletWet():

    file = '/users/global/cornkle/data/OBS/modis_LST/modis_netcdf/lsta_daily_'+str(DATE['year'])+str(DATE['month']).zfill(2)+str(DATE['day']).zfill(2)+'.nc'

    ds = xr.open_dataset(file)

    lsta = ds['LSTA'][0,:,:]
    lsta[lsta < -800] = np.nan
    lsta = lsta-np.nanmean(lsta)


    wav = util.waveletLSTA(lsta.values,3, wet=True)

    wl = wav['power']
    scales = wav['scales']

  #  for id, s in enumerate(scales):
  #      wl[id, :, :][wl[id, :, :] <= s ** .5] = 0

    return wl, scales, lsta



def dominantScaleDry():

    out = '/users/global/cornkle/data/OBS/modis_LST/modis_wav/'

    wll, scales, lsta = run_waveletDry()

    figure = wll[0,:,:].copy()*0

    yyy = []
    xxx = []
    scal = []

    maxoutt = (
        wll == ndimage.maximum_filter(wll, (10, 10, 10), mode='reflect', cval=np.amax(wll) + 1))  # (np.round(orig / 5))

    for nb in range(scales.size):

        wl = wll[nb,:,:]

        orig = scales[nb]
        print(np.round(orig))
        maxout = maxoutt[nb, :, :]

        try:
            yy, xx = np.where((maxout == 1) & (np.abs(lsta) >= 1.) & (wl >= np.percentile(wl[wl >= 0.5], 1)) )#& (wl >= np.percentile(wl[wl >= 0.5], 90)) & (wl > orig ** .5))
        except IndexError:

            continue

        yyy.append(yy)
        xxx.append(xx)



        for y, x in zip(yy, xx):

            ss = orig
            iscale = (np.ceil(ss / 2. / 3.)).astype(int)

            ycirc, xcirc = ua.draw_cut_circle(x, y, iscale, wl)

            figure[ycirc, xcirc] = 1#wl[ycirc, xcirc]
    #
    # f = plt.figure()
    # lsta.plot.contourf()
    # plt.plot((lsta['lon'])[xxx], (lsta['lat'])[yyy], 'bo', markersize=3)

    # da = xr.DataArray(figure, coords={'scale': scales,
    #                                          'lat': lsta['lat'],
    #                                          'lon': lsta['lon']},
    #                   dims=['scale', 'lat', 'lon'])  # .isel(time=0)
    #
    # ds = xr.Dataset({'wav': da})
    #
    # try:
    #     ds.to_netcdf(out+'test_wav.nc')
    # except OSError:
    #     print('Did not find ' + out)
    #     print('Out directory not found')
    # print('Wrote ' + out+'test_wav.nc')

    return xxx, yyy,  figure, lsta, wll

def dominantScaleWet():

    out = '/users/global/cornkle/data/OBS/modis_LST/modis_wav/'

    wll, scales, lsta = run_waveletWet()

    figure = wll[0,:,:].copy()*0

    yyy = []
    xxx = []
    scal = []

    maxoutt = (
        wll == ndimage.maximum_filter(wll, (10, 10, 10), mode='reflect', cval=np.amax(wll) + 1))  # (np.round(orig / 5))

    for nb in range(scales.size):

        wl = wll[nb,:,:]

        orig = scales[nb]
        print(np.round(orig))
        maxout = maxoutt[nb, :, :]

        try:
            yy, xx = np.where((maxout == 1) & (np.abs(lsta) >= 1.) & (wl >= np.percentile(wl[wl >= 0.5], 1)) )#& (wl >= np.percentile(wl[wl >= 0.5], 90)) & (wl > orig ** .5))
        except IndexError:

            continue

        yyy.append(yy)
        xxx.append(xx)



        for y, x in zip(yy, xx):

            ss = orig
            iscale = (np.ceil(ss / 2. / 3.)).astype(int)

            ycirc, xcirc = ua.draw_cut_circle(x, y, iscale, wl)

            figure[ycirc, xcirc] = 1#wl[ycirc, xcirc]
    #
    # f = plt.figure()
    # lsta.plot.contourf()
    # plt.plot((lsta['lon'])[xxx], (lsta['lat'])[yyy], 'bo', markersize=3)

    # da = xr.DataArray(figure, coords={'scale': scales,
    #                                          'lat': lsta['lat'],
    #                                          'lon': lsta['lon']},
    #                   dims=['scale', 'lat', 'lon'])  # .isel(time=0)
    #
    # ds = xr.Dataset({'wav': da})
    #
    # try:
    #     ds.to_netcdf(out+'test_wav.nc')
    # except OSError:
    #     print('Did not find ' + out)
    #     print('Out directory not found')
    # print('Wrote ' + out+'test_wav.nc')

    return xxx, yyy,  figure, lsta, wll


def wav_checkDry():

    nightp = '/users/global/cornkle/MCSfiles/blob_map_30km_-67_JJAS_0-3UTC_centrePoint.nc'
    dayp = '/users/global/cornkle/MCSfiles/blob_map_30km_-67_JJAS_17-19UTC_centrePoint.nc'

    nightp = '/localscratch/wllf030/cornkle/obs_data/blob_maps_MSG/blob_map_35km_-70_0-3UTC.nc'
    dayp = '/localscratch/wllf030/cornkle/obs_data/blob_maps_MSG/blob_map_35km_-70_15-18UTC.nc'


    xxx, yyy, figure, lsta, wll = dominantScaleDry()

    day = xr.open_dataarray(dayp)
    night = xr.open_dataarray(nightp)


    latmin, latmax = (np.min(lsta['lat']), np.max(lsta['lat']))
    lonmin, lonmax = (np.min(lsta['lon']), np.max(lsta['lon']))

    day = day.sel(time=slice(str(DATE['year'])+'-'+str(DATE['month']).zfill(2)+'-'+str(DATE['day']).zfill(2)+' 17:00',str(DATE['year'])+'-'+str(DATE['month']).zfill(2)+'-'+str(DATE['day']).zfill(2)+' 19:00'), lat=slice(latmin, latmax), lon=slice(lonmin,lonmax))
    day=day.sum(dim='time')

    night = night.sel(time=str(DATE['year'])+'-'+str(DATE['month']).zfill(2)+'-'+str(DATE['day']+1).zfill(2), lat=slice(latmin, latmax), lon=slice(lonmin, lonmax))
    night = night.sum(dim='time')
    wllmean0=np.sum(wll,axis=0)

    bla=day.where(day>0)
    bla2 = night.where(night>0)

    f= plt.figure()
    f.add_subplot(2,2,1)
    plt.contourf(lsta['lon'], lsta['lat'], lsta, vmin=-8, vmax=8, nlevels=7, cmap='RdBu_r')
    plt.colorbar()
    plt.contour(lsta['lon'], lsta['lat'], wllmean0, cmap='jet')
    plt.contourf(day['lon'], day['lat'], bla, cmap='viridis')

    plt.title(str(pd.to_datetime(lsta['time'].values))+' DRY anomaly')

    # f.add_subplot(2, 2, 2)
    # plt.contourf(lsta['lon'], lsta['lat'], wllmean0)
    # plt.colorbar()
    # plt.title('Sum of wavelet power - all scales')

    f.add_subplot(2, 2, 2)
    plt.contourf(lsta['lon'], lsta['lat'], lsta, vmin=-8, vmax=8, nlevels=7, cmap='RdBu_r')
    plt.colorbar()
    plt.contour(lsta['lon'], lsta['lat'], wllmean0, cmap='jet')
    plt.contourf(day['lon'], day['lat'], bla2, cmap='viridis')

    plt.title(str(pd.to_datetime(lsta['time'].values)) + ' DRY anomaly')


    f.add_subplot(2, 2, 3)
    plt.contourf(day['lon'], day['lat'], day)
    plt.contour(lsta['lon'], lsta['lat'], wllmean0, cmap='Reds')
    plt.colorbar()
    plt.title('Deep conv. frequency, <70C: Day (16-19UTC), Reds: wavelet power')

    f.add_subplot(2, 2, 4)
    plt.contourf(night['lon'], night['lat'], night)
    plt.contour(lsta['lon'], lsta['lat'], wllmean0, cmap='Reds')
    plt.colorbar()
    plt.title('Deep conv. frequency, <70C: Night (0-3UTC), Reds: wavelet power')




def wav_checkWet():

    nightp = '/users/global/cornkle/MCSfiles/blob_map_30km_-67_JJAS_0-3UTC_centrePoint.nc'
    dayp = '/users/global/cornkle/MCSfiles/blob_map_30km_-67_JJAS_17-19UTC_centrePoint.nc'

    xxx, yyy, figure, lsta, wll = dominantScaleWet()

    day = xr.open_dataarray(dayp)
    night = xr.open_dataarray(nightp)



    latmin, latmax = (np.min(lsta['lat']), np.max(lsta['lat']))
    lonmin, lonmax = (np.min(lsta['lon']), np.max(lsta['lon']))

    day = day.sel(time=slice('2006-06-19','2006-06-20'), lat=slice(latmin, latmax), lon=slice(lonmin,lonmax))
    day=day.sum(dim='time')

    night = night.sel(time=slice('2006-07-20', '2006-07-21'), lat=slice(latmin, latmax), lon=slice(lonmin, lonmax))
    night = night.sum(dim='time')
    wllmean0=np.sum(wll,axis=0)

    bla=day.where(day>0)

    f= plt.figure()
    f.add_subplot(2,2,1)
    plt.contourf(lsta['lon'], lsta['lat'], lsta, vmin=-8, vmax=8, nlevels=7, cmap='RdBu_r')
    plt.colorbar()
    plt.contour(lsta['lon'], lsta['lat'], wllmean0, cmap='jet')
    plt.contourf(day['lon'], day['lat'], bla, cmap='viridis')

    plt.title(str(pd.to_datetime(lsta['time'].values))+' WET anomaly')

    f.add_subplot(2, 2, 2)
    plt.contourf(lsta['lon'], lsta['lat'], wllmean0)
    plt.colorbar()
    plt.title('Sum of wavelet power - all scales')

    f.add_subplot(2, 2, 3)
    plt.contourf(day['lon'], day['lat'], day)
    plt.contour(lsta['lon'], lsta['lat'], wllmean0, cmap='Reds')
    plt.colorbar()
    plt.title('Deep conv. frequency, <70C: Day (16-19UTC), Reds: wavelet power')

    f.add_subplot(2, 2, 4)
    plt.contourf(night['lon'], night['lat'], night)
    plt.contour(lsta['lon'], lsta['lat'], wllmean0, cmap='Reds')
    plt.colorbar()
    plt.title('Deep conv. frequency, <70C: Night (0-3UTC), Reds: wavelet power')
