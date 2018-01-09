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


DATE = {'day' : 22,
        'month' : 8,
        'year' : 2008}

def run_waveletDry():


    file = '/users/global/cornkle/data/OBS/modis_LST/modis_netcdf/lsta_daily_'+str(DATE['year'])+str(DATE['month']).zfill(2)+str(DATE['day']).zfill(2)+'.nc'

    ds = xr.open_dataset(file)
    ds = ds.sel(lon=slice(-10,10), lat=slice(10,20))
    lsta = ds['LSTA'][0,:,:]

    lsta[lsta<-800] = np.nan
    points = np.where(np.isfinite(lsta))
    inter = np.where(np.isnan(lsta))

    # interpolate over sea from land points
    lsta[inter] = griddata(points, np.ravel(lsta.values[points]), inter, method='nearest')

    lsta = lsta #- np.nanmean(lsta)

    f = plt.figure()
    plt.imshow(lsta)

    wav = util.waveletLSTA(lsta.values,3, method='dry')

    wl = wav['power']
    f = plt.figure()
    plt.imshow(np.nansum(wl,axis=0))
    #wl[:, inter[0], inter[1]] = np.nan
    scales = wav['scales']
    dom = wav['dominant']
    dom[inter]=np.nan
    print(scales)

  #  for id, s in enumerate(scales):
  #      wl[id, :, :][wl[id, :, :] <= s ** .5] = 0

    return wl, scales, lsta, inter, dom



def wav_checkDry():

    start_time = time.time()

    nightp = '/users/global/cornkle/MCSfiles/blob_map_30km_-67_JJAS_0-3UTC_centrePoint.nc'
    dayp = '/users/global/cornkle/MCSfiles/blob_map_30km_-67_JJAS_17-19UTC_centrePoint.nc'

    nightp = '/localscratch/wllf030/cornkle/obs_data/blob_maps_MSG/blob_map_35km_-70_0-3UTC.nc'
    dayp = '/localscratch/wllf030/cornkle/obs_data/blob_maps_MSG/blob_map_35km_-70_15-18UTC.nc'

    wll, scales, lsta, inter, dom = run_waveletDry()

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
    posi = 250  # 116 ## 118

    plt.figure()
    plt.imshow(dom)
    plt.colorbar()
    plt.gca().invert_yaxis()

    f= plt.figure()
    f.add_subplot(2,2,1)
    plt.contourf(lsta['lon'], lsta['lat'], lsta, vmin=-8, vmax=8, nlevels=7, cmap='RdBu_r')
    plt.axhline(lsta['lat'][posi], linestyle='--', linewidth=2,
             color='black')
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

    f = plt.figure(figsize=(6.5, 11), dpi=300)
    #
    gridspec.GridSpec(3, 1)

    ax1 = plt.subplot2grid((3, 1), (0, 0), rowspan=2)
    ax3 = plt.subplot2grid((3, 1), (2, 0))
    #
  #  mt = ax1.contourf(np.arange(wll.shape[2]) * 3, np.arange(wll.shape[1]) * 5, lsta, cmap='Greys')
    ax1.plot(np.arange(wll.shape[2]) * 3, lsta[posi, :], linewidth=2,
             color='black')
    ax1.set_xlim(0, wll.shape[2]*3)
    mp = ax3.contourf(np.arange(wll.shape[2]) * 3, scales, wll[:, posi, :], cmap='viridis', levels=[0,0.25, 0.5 , 1, 2, 4, 6, 8])
   # plt.colorbar(mp)
    # levels=[0,1, 2,5,10,20,40,80,100, 130, 150, 180, 200, 300,400]

    # for p1, p2 in zip(ppos[1], ppos[0]):
    #    ax3.errorbar((np.arange(wll.shape[2])*5)[p1], arr[p2], xerr=arr[p2]/2, fmt='o', ecolor='white', color='white', capthick=3, ms=3, elinewidth=0.7)
    # ax3.set_xlim(100,700)
    ax3.set_ylim(15, 180)
    ax3.set_xlabel('Spatial extent (km)')
    ax3.set_ylabel('Length scale (km)')

    plt.tight_layout()

    f.subplots_adjust(right=0.86)


    plt.show()

    print('Elapsed time: ', time.time()-start_time)


