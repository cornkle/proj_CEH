import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from wavelet import util
import pdb

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

DATE = {'day' : 2,
        'month' : 9,
        'year' : 2008}

def run_waveletDry():


    file = '/users/global/cornkle/data/OBS/MSG_LSTA/lsta_netcdf/lsta_daily_'+str(DATE['year'])+str(DATE['month']).zfill(2)+str(DATE['day']).zfill(2)+'.nc'

    ds = xr.open_dataset(file)
    ds = ds.sel(lon=slice(-10,10), lat=slice(10,18))

    lsta = ds['LSTA'].squeeze()
    cells = ds['cell'].squeeze()

    lsta_raw = lsta.copy()
    #
    # plt.figure()
    # plt.imshow(lsta_raw, origin='lower')


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

    ds = xr.open_dataset(file)

    ds = ds.sel(lon=slice(-10,10), lat=slice(10,18))

    lsta2 = ds['LSTA'][0,:,:]


    wll, scales, power, lsta, inter, cells, lsta_raw = run_waveletDry()

    daystring = str(DATE['year'])+'-'+str(DATE['month']).zfill(2)+'-'+str(DATE['day']).zfill(2)
    #smarr = sm['SM'].sel(time=daystring)

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

    cmap = u_plot.discrete_cmap(8, base_cmap='gist_ncar')  #gist_ncar

    era = xr.open_dataset(constants.ERA_DAILY_SRFC_ANO)
    era = era.sel(latitude=slice(None, None, -1))

    eday = era.sel(longitude=slice(-10, 10), latitude=slice(10, 18),
                   time=str(DATE['year']) + str(DATE['month']).zfill(2) + str(DATE['day']).zfill(2))

    u = eday['u10'].squeeze()
    v = eday['v10'].squeeze()
    div = eday['ws10'].squeeze()


    # ws, wd = u_met.u_v_to_ws_wd(u.values,v.values)
    # eday['ws10'] = (('time', 'latitude', 'longitude'), ws[None,...])
    # eday['wd10'] = (('time', 'latitude', 'longitude'), wd[None, ...])
    ws_on_lsta = ds.salem.transform(div)
    u_on_lsta = ds.salem.transform(u)
    u_on_lsta = u_on_lsta[4::50, 4::50]
    v_on_lsta = ds.salem.transform(v)
    v_on_lsta = v_on_lsta[4::50, 4::50]

    xe, ye = np.meshgrid(u_on_lsta['lon'], u_on_lsta['lat'])

    # transform their coordinates to the map reference system and plot the arrows
    xx, yy = map.grid.transform(xe, ye, crs=u.salem.grid.proj)


    f = plt.figure()
    plt.imshow(div, origin='lower')

    xi, yi = np.where(np.isfinite(cells))

    for r in [3,6,9,12,15,18,21,24]:
        pos = np.where((cells<=r) & (cells>r-3))
        if r == 24:
            cells[pos] = 0
        else:
            cells[pos] = r

    pdb.set_trace()
    data = ppower_small

    f= plt.figure()

    ax1=f.add_subplot(2,2,1)

    map.set_data(lsta_raw)
    map.set_plot_params(cmap='RdBu_r', levels=[-8,-6,-4,-2,-1,-0.5,0.5,1,2,4,6,8])
    #map.set_contour(ppower_small, color='r', levels=[1,3,5])
    map.visualize(ax=ax1, addcbar=False )
    cax = ax1.scatter(xl, yl, c=cells, cmap=cmap, s=5)
    cbar = plt.colorbar(cax, ticks=[0,3,6,9,12,15,18,21])
    cbar.set_ticklabels(['1-3','4-6','7-9','10-12','13-15', '16-18', '19-21', '22-0'])

    plt.title(str(pd.to_datetime(lsta['time'].values))+'Seasonal land surface temperature anomaly')

    # ax2 = f.add_subplot(2, 2, 2)
    #
    # #plt.contourf(day_small['lon'], day_small['lat'], day_small)
    # plt.pcolormesh(xe, ye, div, cmap='RdBu', vmin=-0.005, vmax=0.005)
    # plt.colorbar()
    # # plt.scatter(xe, ye, c=cells, cmap='Greys', s=7)
    # plt.scatter(xv, yv, c=cells, cmap=cmap, s=5)
    # #plt.scatter(xv, yv, c=cells, cmap='jet', s=5)
    # plt.title('Small')

    # ax2 = f.add_subplot(2, 2, 2)
    # map.set_plot_params(cmap='viridis', vmin=0, vmax=5)
    # map.set_data(ppower_small)
    # map.visualize(ax=ax2)
    # ax2.scatter(xl, yl, c=cells, cmap=cmap, s=5)
    # #plt.scatter(xv, yv, c=cells, cmap='jet', s=5)
    # plt.title('Scales 9 - 30km')


    ax3 = f.add_subplot(2, 2, 3)
    map.set_plot_params(cmap='viridis', vmin=0, vmax=5)
    map.set_data(ppower_mid)
    map.visualize(ax=ax3)
    ax3.scatter(xl, yl, c=cells, cmap=cmap, s=5)
    #plt.scatter(xv, yv, c=cells, cmap='jet', s=5)
    plt.title('Scales 31 - 60km')


    ax4 = f.add_subplot(2, 2, 4)
    map.set_plot_params(cmap='viridis', vmin=0, vmax=5)
    map.set_data(ppower_small)
    map.visualize(ax=ax4)
    ax4.scatter(xl, yl, c=cells, cmap=cmap, s=5)

    plt.title('Scales 80 - 250km')


    ax2=f.add_subplot(2,2,2)

    map.set_data(ws_on_lsta)
    map.set_plot_params(cmap='RdBu_r', vmin=-5, vmax=5)
    map.visualize(ax=ax2, addcbar=True )
    cax = ax1.scatter(xl, yl, c=cells, cmap=cmap, s=5)

    qu = ax2.quiver(xx,yy,u_on_lsta,v_on_lsta, scale=60)
    qk = plt.quiverkey(qu, 0.7, 0.95, 3, '3 m s$^{-1}$',
                       labelpos='E', coordinates='figure')





    # ax2 = f.add_subplot(2, 2, 2)
    #
    # #plt.contourf(day_small['lon'], day_small['lat'], day_small)
    # plt.contourf(lsta['lon'], lsta['lat'], ppower_small, cmap='RdBu_r', extend='both')
    # plt.colorbar()
    # plt.scatter(xv, yv, c=cells, cmap=cmap, s=5)
    # #plt.scatter(xv, yv, c=cells, cmap='jet', s=5)
    # plt.title('Small')
    #
    # plt.title(str(pd.to_datetime(lsta['time'].values)) + ' Land surface temperature anomaly, +/-4')
    #
    #
    # f.add_subplot(2, 2, 3)
    # #plt.contourf(day_small['lon'], day_small['lat'], day_small)
    # plt.contourf(lsta['lon'], lsta['lat'], ppower_mid, cmap='RdBu_r',vmin=-1., vmax=1, extend='both')
    # plt.colorbar()
    # plt.scatter(xv, yv, c=cells, cmap=cmap, s=5)
    # #plt.scatter(xv, yv, c=cells, cmap='jet', s=5)
    # plt.title('Mid')
    #
    #
    # f.add_subplot(2, 2, 4)
    # #plt.contourf(night['lon'], night['lat'], night)
    # plt.contourf(lsta['lon'], lsta['lat'], ppower_big, cmap='RdBu_r', vmin=-1, vmax=1, extend='both')
    # plt.colorbar()
    # plt.scatter(xv, yv, c=cells, cmap=cmap, s=5)
    # #plt.scatter(xv, yv, c=cells, cmap='jet', s=5)
    # #plt.gca().invert_yaxis()
    #
    # plt.title('Big')


    # print(np.sum(data[(cells.values>=17) | (cells.values<=21)]>-0.1) / np.sum(((cells.values>=17) | (cells.values<=21)) & (np.isfinite(data))))
    #
    # print('Data dist', np.sum(data>-0.1)/ np.sum(np.isfinite(data)))
    #
    # f = plt.figure()
    # ax = f.add_subplot(131)
    # plt.hist(data[np.isfinite(data)], bins=6)
    # ax = f.add_subplot(132)
    # plt.hist(data[(np.isfinite(data)) & np.isfinite(cells.values)], bins=6)
    #
    # ax = f.add_subplot(133)
    # plt.hist(ppower_small[(np.isfinite(ppower_small))], bins=6)
    #
    #
    # pos = np.where(np.isfinite(cells.values))
    # isgrad = []
    # for y, x in zip(pos[0], pos[1]):
    #     ycirc, xcirc = u_arrays.draw_circle(x,y,1)
    #     vals = data[ycirc,xcirc]
    #     if np.sum(np.isfinite(vals)) < vals.size*0.5:
    #         continue
    #     if (np.min(vals)<0) & (np.max(vals)>0):
    #         isgrad.append(1)
    #     else:
    #         isgrad.append(0)
    #
    # isgrad = np.array(isgrad)
    #
    # pdb.set_trace()








