import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import pdb
from utils import u_arrays as ua, constants
import pandas as pd
import time
from utils import u_plot
import salem

#nice example 25.06.2006

#24.09.2008 and following day wet/dry place makes no sense
#12.08.2009 is cool
#2.8.2009
#29.07.2007
#26.082008 MCS in the morning

DATE = {'day' : 29,
        'month' :8,
        'year' : 2008}

def run_waveletDry():


    file = '/users/global/cornkle/data/OBS/MSG_LSTA/lsta_netcdf_new/lsta_daily_'+str(DATE['year'])+str(DATE['month']).zfill(2)+str(DATE['day']).zfill(2)+'.nc'

    ds = xr.open_dataset(file)
    ds = ds.sel(lon=slice(-10,10), lat=slice(10,20))

    lsta = ds['LSTA'].squeeze()
    cells = ds['cell'].squeeze()

    # lsta.values[inter1] = griddata(points, np.ravel((lsta.values)[points]), inter1, method='linear')
    # inter = np.where(np.isnan(lsta))
    # lsta.values[inter] = griddata(points, np.ravel(lsta.values[points]), inter, method='nearest')

    return lsta, cells



def wav_checkDry():


    lsta, cells = run_waveletDry()


    xv, yv = np.meshgrid(lsta['lon'], lsta['lat'])

    cmap = u_plot.discrete_cmap(24, base_cmap='gist_ncar')

    era = xr.open_dataset(constants.ERA5+'ERA5_2008_12UTCpl.nc')
    eday = era.sel(time=str(DATE['year']) + str(DATE['month']).zfill(2) + str(DATE['day']).zfill(2), level=[950, 700])

    era_on_lsta = lsta.salem.transform(eday)

    shear = era_on_lsta['u'].sel(level=700).values - era_on_lsta['u'].sel(level=950).values
    era_on_lsta = era_on_lsta.sel(level=950)
    f= plt.figure()

    ax1=f.add_subplot(2,2,1)

    map = lsta.salem.get_map()
    xl, yl = map.grid.transform(xv, yv)
    st=20
    xquiv = xl[4::st, 4::st]
    yquiv = yl[4::st, 4::st]

    map.set_data(lsta)

    map.set_plot_params(cmap='RdBu_r', vmin=-8, vmax=8, extend='both')
    map.visualize(ax=ax1, addcbar=False )
    cax = ax1.scatter(xl, yl, c=cells, cmap=cmap, s=5)
    #cbar = plt.colorbar(cax, ticks=np.arange(0,24))
    #cbar.set_ticklabels(np.array(np.arange(0,24), dtype=str))


    plt.title(str(pd.to_datetime(lsta['time'].values))+' LSTA')


    ax2 = f.add_subplot(2, 2, 2)
    # plt.contourf(night['lon'], night['lat'], night)
    tera = shear #era_on_lsta['t'].squeeze()-273.15
    map.set_data( (tera  ))

    map.set_plot_params(cmap='jet',  extend='both', vmin=-25, vmax=-10) #vmin=-0.1, vmax=0.1,
    #map.set_plot_params(cmap='RdBu_r', vmin=-8, vmax=8, extend='both')
    map.visualize(ax=ax2, addcbar=True)

    plt.scatter(xl, yl, c=cells, cmap=cmap, s=5)

    u = era_on_lsta['u'].squeeze()[4::st, 4::st]
    v = era_on_lsta['v'].squeeze()[4::st, 4::st]

    qu = ax2.quiver(xquiv,yquiv, u.values,v.values, scale=120)
    cax = ax2.scatter(xl, yl, c=cells, cmap=cmap, s=5)
    # qk = plt.quiverkey(qu, 0.7, 0.95, 3, '3 m s$^{-1}$',
    #                    labelpos='E', coordinates='figure')

    ax2 = f.add_subplot(2, 2, 3)
    # plt.contourf(night['lon'], night['lat'], night)
    tera =  era_on_lsta['d'].squeeze()*1000
    map.set_data((tera))

    map.set_plot_params(cmap='RdBu', extend='both', vmin=-0.1, vmax=0.1)  # vmin=-0.1, vmax=0.1,
    # map.set_plot_params(cmap='RdBu_r', vmin=-8, vmax=8, extend='both')
    map.visualize(ax=ax2, addcbar=True)

    plt.scatter(xl, yl, c=cells, cmap=cmap, s=5)

    u = era_on_lsta['u'].squeeze()[4::st, 4::st]
    v = era_on_lsta['v'].squeeze()[4::st, 4::st]

    qu = ax2.quiver(xquiv, yquiv, u.values, v.values, scale=120)
    cax = ax2.scatter(xl, yl, c=cells, cmap=cmap, s=5)
    # qk = plt.quiverkey(qu, 0.7, 0.95, 3, '3 m s$^{-1}$',
    #                    labelpos='E', coordinates='figure')





