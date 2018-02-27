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


DATE = {'day' : 28,
        'month' : 9,
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


    wav = util.waveletLSTA(lsta.values,3, method='dry')

    wl = wav['power']

    wl[:,inter[0], inter[1]] = np.nan
    wl[:, inter1[0], inter1[1]] = np.nan
    f = plt.figure()
    plt.imshow(lsta, cmap='RdBu', vmin=-40, vmax=40)
    wl[-5::, inter[0], inter[1]] = np.nan
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
    wllmean0=np.nansum(wll[-5::],axis=0)
    print(scales[0:5])

    bla=day_small.where(day_small>0)
    bla2 = night.where(night>0)
    posi = 50  # 116 ## 118

    f= plt.figure()
    f.add_subplot(2,2,1)
    plt.contourf(lsta['lon'], lsta['lat'], lsta, vmin=-8, vmax=8, nlevels=7, cmap='RdBu_r')
    plt.contour(day_small['lon'], day_small['lat'], day_small, cmap='prism')
    plt.axhline(lsta['lat'][posi], linestyle='--', linewidth=2,
             color='black')
    plt.contour(lsta['lon'], lsta['lat'], wllmean0, cmap='RdBu', levels=[-20, -15, -10, -5, 5, 10, 15,20])
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
    plt.contour(lsta['lon'], lsta['lat'], wllmean0, cmap='RdBu', levels=[-20, -15, -10, -5, 5, 10, 15,20])

    plt.contourf(day_small['lon'], day_small['lat'], bla2, cmap='viridis')

    plt.title(str(pd.to_datetime(lsta['time'].values)) + ' DRY anomaly NIGHTTIME')


    f.add_subplot(2, 2, 3)
    plt.contourf(day_small['lon'], day_small['lat'], day_small)
    plt.contour(lsta['lon'], lsta['lat'], wllmean0, cmap='RdBu', levels=[-20, -15, -10, -5, 5, 10, 15,20])
    plt.colorbar()
    plt.title('Deep conv. frequency, <70C: '+daystring+', Reds: wavelet power')

    f.add_subplot(2, 2, 4)
    plt.contourf(night['lon'], night['lat'], night)
    plt.contour(lsta['lon'], lsta['lat'], wllmean0, cmap='RdBu', levels=[-20, -15, -10, -5, 5, 10, 15,20])
    #plt.gca().invert_yaxis()
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
    mp = ax3.contourf(np.arange(wll.shape[2]) * 3, scales, wll[:, posi, :], cmap='RdBu', levels=[-8,-6,-4,-2,-1,-0.5,-0.25,0,0.25, 0.5 , 1, 2, 4, 6, 8])

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


