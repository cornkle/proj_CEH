import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from wavelet import util
import pdb
import matplotlib
import pandas as pd
import time
import matplotlib.gridspec as gridspec
from scipy.interpolate import griddata
import numpy.ma as ma
from utils import u_plot, u_met, constants, u_arrays
import salem as sm

#nice example 25.06.2006

#24.09.2008 and following day wet/dry place makes no sense
#12.08.2009 is cool
#2.8.2009
#29.07.2007

#27.09.2008 works

# 3.08.2007: ERA-I 10m wind anomaly matches

#28.06.2006 is cool

DATE = {'day' : 28,
        'month' : 6,
        'year' : 2006}

def run_waveletDry():


    file = '/users/global/cornkle/data/OBS/MSG_LSTA/lsta_netcdf/lsta_daily_'+str(DATE['year'])+str(DATE['month']).zfill(2)+str(DATE['day']).zfill(2)+'.nc'

    ds = xr.open_dataset(file)
    ds = ds.sel(lon=slice(-10,10), lat=slice(10,20))

    lsta = ds['LSTA'].squeeze()
    cells = ds['cell'].squeeze()

    lsta_raw = lsta.copy()


    points = np.where(np.isfinite(lsta_raw.values))
    inter1 = np.where(np.isnan(lsta_raw.values))

    # lsta_raw.values[inter1[0], inter1[1]] = 90
    # plt.figure()
    # plt.imshow(lsta_raw, origin='lower')
    # return
    # lsta.values[inter1] = griddata(points, np.ravel((lsta.values)[points]), inter1, method='linear')
    # inter = np.where(np.isnan(lsta))
    # lsta.values[inter] = griddata(points, np.ravel(lsta.values[points]), inter, method='nearest')
    lsta[inter1]=0
    # f = plt.figure()
    # plt.imshow(lsta.values)
    wav = util.waveletLSTA_both(lsta.values,3)

    #wl = wav['dominant']
    power = wav['power_dry']

    power[:, inter1[0], inter1[1]] = np.nan
    #wl[inter1[0], inter1[1]] = np.nan
    wl = 0

    scales = wav['scales']

    print(scales)

    # for id, s in enumerate(scales):
    #     (power[id, :, :])[np.abs(power[id, :, :]) <= s**0.5] = 0

    return wl, scales, power, lsta, inter1, cells, lsta_raw



def wav_checkDry():

    start_time = time.time()

    # smfile = '/users/global/cornkle/data/OBS/AMSRE/day_aqua/amsre_monthly_anomaly.nc'
    # sm = xr.open_dataset(smfile)
    # sm = sm.sel(lon=slice(-10, 10), lat=slice(10, 18))


    file = '/users/global/cornkle/data/OBS/MSG_LSTA/lsta_netcdf/lsta_daily_'+str(DATE['year'])+str(DATE['month']).zfill(2)+str(DATE['day']).zfill(2)+'.nc'
    ftopo = '/users/global/cornkle/data/pythonWorkspace/proj_CEH/topo/gtopo_1min_afr.nc'
    msg_t = '/users/global/cornkle/MCSfiles/blob_map_MCSs_-50_JJAS_points_dominant.nc'
    ds = xr.open_dataset(file)
    top = xr.open_dataarray(ftopo)
    msg = xr.open_dataarray(msg_t)
    msg.name = 'msg'

    pdb.set_trace()
    ds = ds.sel(lon=slice(-10,10), lat=slice(10,20))
    top = top.sel(lon=slice(-10, 10), lat=slice(10, 20))
    msg = msg.sel(lon=slice(-10, 10), lat=slice(10, 20), time=slice('2006-06-28T17:00:00','2006-06-29T15:00:00'))
    msg_mins = msg.min(dim='time')
    lsta2 = ds['LSTA'][0,:,:]

    wll, scales, power, lsta, inter, cells, lsta_raw = run_waveletDry()

    daystring = str(DATE['year'])+'-'+str(DATE['month']).zfill(2)+'-'+str(DATE['day']).zfill(2)

    wllmean0=wll

    ppower_small = np.sum(power[0:4,:,:], axis=0)
    ppower_mid = np.sum(power[4:7, :, :], axis=0)
    ppower_big = np.sum(power[7::, :, :], axis=0)
    print('Small', scales[0:4])
    print('Mid', scales[4:7])
    print('Big', scales[7::])

    map = lsta.salem.get_map()

    xv, yv = np.meshgrid(lsta['lon'], lsta['lat'])
    xl, yl = map.grid.transform(xv, yv)

    cmap = u_plot.discrete_cmap(8, base_cmap='gnuplot2')  #gist_ncar
    cmap_back = u_plot.discrete_cmap(8, base_cmap='Greys')

    # era = xr.open_dataset(constants.ERA_DAILY_PL)
    # era = era.sel(latitude=slice(None, None, -1))
    #
    # eday = era.sel(
    #                time=str(DATE['year']) + str(DATE['month']).zfill(2) + str(DATE['day']).zfill(2), level=850)
    #
    # u = eday['u'].squeeze()
    # v = eday['v'].squeeze()
    # div = eday['d'].squeeze()
    #
    # u_plot.quick_map(div*1000, vmin=-0.01, vmax=0.01, cmap='RdBu')
    #
    #
    # # ws, wd = u_met.u_v_to_ws_wd(u.values,v.values)
    # # eday['ws10'] = (('time', 'latitude', 'longitude'), ws[None,...])
    # # eday['wd10'] = (('time', 'latitude', 'longitude'), wd[None, ...])
    # ws_on_lsta = ds.salem.transform(div)*1000
    # u_on_lsta = ds.salem.transform(u)
    # u_on_lsta = u_on_lsta[4::50, 4::50]
    # v_on_lsta = ds.salem.transform(v)
    # v_on_lsta = v_on_lsta[4::50, 4::50]
    # xe, ye = np.meshgrid(u_on_lsta['lon'], u_on_lsta['lat'])
    # xx, yy = map.grid.transform(xe, ye, crs=u.salem.grid.proj)

    # transform their coordinates to the map reference system and plot the arrows


    xi, yi = np.where(np.isfinite(cells))

    topo_on_lsta = ds.salem.lookup_transform(top)
    grad = np.gradient(topo_on_lsta.values)
    gradsum = abs(grad[0]) + abs(grad[1])
    msg_on_lsta = ds.salem.transform(msg_mins)
    topo_on_lsta.values[gradsum > 30] = np.nan
    topo_on_lsta.values[topo_on_lsta>400] = np.nan



    for rh in [3,6,9,12,15,18,21,24]:
        pos = np.where((cells<=rh) & (cells>rh-3))

        if rh == 24:
            cells[pos] = 0
        else:
            cells[pos] = rh

    cmap2 = plt.get_cmap('RdBu_r')
    cmap2.set_bad('Grey', 1.)

    wavelet = plt.get_cmap('PuRd')
    wavelet.set_bad('Grey', 1.)

    # topobad = plt.get_cmap('topo')
    # topobad.set_bad('Grey', 1.)

    f= plt.figure(figsize=(13,7))

    ax1=f.add_subplot(2,2,1)
    masked = np.ma.masked_invalid(lsta_raw)
    map.set_data(masked, interp='linear')
    map.set_plot_params(cmap=cmap2, levels=[-6,-4,-2,-1,-0.5,0.5,1,2,4,6])
    #map.set_lonlat_contours(interval=0)

    map.visualize(ax=ax1, cbar_title='K', addcbar=True, title='28-06-2006: 0630-1600UTC LSTA & afternoon/nighttime cores')
    cax = ax1.scatter(xl, yl, c=cells, cmap=cmap_back, s=15)
    cax = ax1.scatter(xl, yl, c=cells, cmap=cmap, s=6)
    cbar = plt.colorbar(cax, ticks=[0,3,6,9,12,15,18,21], fraction=0.07, pad=0.5, orientation='horizontal')
    cbar.set_ticklabels([ '22-0','1-3','4-6','7-9','10-12','13-15', '16-18', '19-21'])
    cbar.set_label('Hours of day')
    ax1.get_xaxis().set_visible(False)

    # print('Topo min', np.min(topo_on_lsta))
    ax2 = f.add_subplot(2, 2, 2)
    map.set_plot_params(cmap='topo', vmin=np.nanmin(topo_on_lsta), vmax=1000, extend='both')
    masked = np.ma.masked_invalid(topo_on_lsta)
    map.set_data(masked)
    map.set_contour(topo_on_lsta, levels=np.arange(400,1000,25), cmap='Greys_r')
    map.set_contour(msg_on_lsta, levels=np.arange(-80,-40,5), cmap='jet_r')
    map.visualize(ax=ax2, addcbar=True )
    map.visualize(ax=ax2, cbar_title='m', title='Domain topography, total height difference: 342m')
    #cax = ax2.scatter(xl, yl, c=cells, cmap=cmap_back, s=20)
    #cax = ax2.scatter(xl, yl, c=cells, cmap=cmap, s=10)
    #plt.scatter(xv, yv, c=cells, cmap='jet', s=5)
    map.set_contour()

    ax3 = f.add_subplot(2, 2, 3)
    map.set_plot_params(cmap=wavelet, vmin=0, vmax=5)
    map.set_data(ppower_small)
    map.visualize(ax=ax3, cbar_title='Wavelet power', title='Surface scales 9 - 30km')
    #ax3.scatter(xl, yl, c=cells, cmap=cmap, s=5)
    #plt.scatter(xv, yv, c=cells, cmap='jet', s=5)

    ax4 = f.add_subplot(2, 2, 4)
    map.set_plot_params(cmap=wavelet, vmin=0, vmax=5)
    map.set_data(ppower_big)
    map.visualize(ax=ax4, cbar_title='Wavelet power', title='Surface scales 80 - 250km')
    #ax4.scatter(xl, yl, c=cells, cmap=cmap, s=5)

    plt.tight_layout()
    plt.savefig('/users/global/cornkle/figs/LSTA-bullshit/GEWEX/example_map.png')


    # ax2=f.add_subplot(2,2,2)
    # map.set_data(msg_on_lsta, interp='linear')
    # map.set_plot_params(cmap='viridis', vmin=-80, vmax=-70)
    # map.visualize(ax=ax2, addcbar=True )
    #cax = ax1.scatter(xl, yl, c=cells, cmap=cmap, s=5)

    # qu = ax2.quiver(xx,yy,u_on_lsta,v_on_lsta, scale=60)
    # qk = plt.quiverkey(qu, 0.7, 0.95, 3, '3 m s$^{-1}$',
    #                    labelpos='E', coordinates='figure')













