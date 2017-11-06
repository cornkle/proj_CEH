import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from wavelet import util
import pdb
from scipy import ndimage
from utils import u_arrays as ua
import pandas as pd
import time


DATE = {'day' : 22,
        'month' : 8,
        'year' : 2008}

def run_waveletDry():


    file = '/users/global/cornkle/data/OBS/modis_LST/modis_netcdf/lsta_daily_'+str(DATE['year'])+str(DATE['month']).zfill(2)+str(DATE['day']).zfill(2)+'.nc'

    ds = xr.open_dataset(file)
    ds = ds.sel(lon=slice(-10,10), lat=slice(10,20))
    lsta = ds['LSTA'][0,:,:]

    lsta[lsta<-800] = np.nan
    lsta = lsta #- np.nanmean(lsta)

    wav = util.waveletLSTA(lsta.values,3, method='dry')

    wl = wav['power']
    scales = wav['scales']
    print(scales)

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

    return xxx, yyy,  figure, lsta, wll, scales



def wav_checkDry():

    start_time = time.time()

    nightp = '/users/global/cornkle/MCSfiles/blob_map_30km_-67_JJAS_0-3UTC_centrePoint.nc'
    dayp = '/users/global/cornkle/MCSfiles/blob_map_30km_-67_JJAS_17-19UTC_centrePoint.nc'

    nightp = '/localscratch/wllf030/cornkle/obs_data/blob_maps_MSG/blob_map_35km_-70_0-3UTC.nc'
    dayp = '/localscratch/wllf030/cornkle/obs_data/blob_maps_MSG/blob_map_35km_-70_15-18UTC.nc'


    xxx, yyy, figure, lsta, wll, scales = dominantScaleDry()

    day = xr.open_dataarray(dayp)
    night = xr.open_dataarray(nightp)



    latmin, latmax = (np.min(lsta['lat']), np.max(lsta['lat']))
    lonmin, lonmax = (np.min(lsta['lon']), np.max(lsta['lon']))
    daystring = str(DATE['year'])+'-'+str(DATE['month']).zfill(2)+'-'+str(DATE['day']).zfill(2)
    day = day.sel(time=slice(daystring+' 17:00',daystring+' 19:00'), lat=slice(latmin, latmax), lon=slice(lonmin,lonmax))
    day=day.sum(dim='time')
    nightstring = str(DATE['year'])+'-'+str(DATE['month']).zfill(2)+'-'+str(DATE['day']+1).zfill(2)
    night = night.sel(time=nightstring, lat=slice(latmin, latmax), lon=slice(lonmin, lonmax))
    night = night.sum(dim='time')
    wllmean0=np.sum(wll[0:5],axis=0)
    print(scales[0:5])

    bla=day.where(day>0)
    bla2 = night.where(night>0)

    f= plt.figure()
    f.add_subplot(2,2,1)
    plt.contourf(lsta['lon'], lsta['lat'], lsta, vmin=-8, vmax=8, nlevels=7, cmap='RdBu_r')

   # plt.contour(lsta['lon'], lsta['lat'], wllmean0, cmap='RdBu')
    plt.colorbar()
    plt.contourf(day['lon'], day['lat'], bla, cmap='viridis')

    plt.title(str(pd.to_datetime(lsta['time'].values))+' DRY anomaly')

    # f.add_subplot(2, 2, 2)
    # plt.contourf(lsta['lon'], lsta['lat'], wllmean0)
    # plt.colorbar()
    # plt.title('Sum of wavelet power - all scales')

    f.add_subplot(2, 2, 2)
    plt.contourf(lsta['lon'], lsta['lat'], lsta, vmin=-8, vmax=8, nlevels=7, cmap='RdBu_r')
    plt.colorbar()
    plt.contour(lsta['lon'], lsta['lat'], wllmean0, cmap='RdBu')

    plt.contourf(day['lon'], day['lat'], bla2, cmap='viridis')

    plt.title(str(pd.to_datetime(lsta['time'].values)) + ' DRY anomaly')


    f.add_subplot(2, 2, 3)
    plt.contourf(day['lon'], day['lat'], day)
    plt.contourf(lsta['lon'], lsta['lat'], wllmean0, cmap='RdBu', levels=[-20, -15, -10, -5, 5, 10, 15,20])
    plt.colorbar()
    plt.title('Deep conv. frequency, <70C: '+daystring+', Reds: wavelet power')

    f.add_subplot(2, 2, 4)
    plt.contourf(night['lon'], night['lat'], night)
    plt.imshow(wllmean0, cmap='RdBu_r', vmin=-20, vmax=20)
    plt.gca().invert_yaxis()
    plt.colorbar()
    plt.title('Deep conv. frequency, <70C: '+nightstring+' , Reds: wavelet power')

    print('Elapsed time: ', time.time()-start_time)


