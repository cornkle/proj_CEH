import numpy as np
import xarray as xr
from wavelet import util
import os
import itertools
import pdb


def run_netcdf():

    for d, m, y in itertools.product(np.arange(1,32), np.arange(6,10), np.arange(2006,2011)):


        DATE = {'day': d,
                'month': m,
                'year': y}


        file = '/users/global/cornkle/data/OBS/MSG_LSTA/lsta_netcdf/lsta_daily_'+str(DATE['year'])+str(DATE['month']).zfill(2)+str(DATE['day']).zfill(2)+'.nc'
        outfile = '/users/global/cornkle/data/OBS/MSG_LSTA/power_maps/lsta_daily_powerMaxscale_'+str(DATE['year'])+str(DATE['month']).zfill(2)+str(DATE['day']).zfill(2)+'.nc'

        ds = xr.open_dataset(file)
        #ds = ds.sel(lon=slice(-10, 10), lat=slice(10, 18))

        lsta = ds['LSTA'].squeeze()
        container = ds.copy()

        inter1 = np.where(np.isnan(lsta.values))
        lsta[inter1] = 0

        wav = util.LSTA_maxPowerScale(lsta.values, 3)

        wl = wav['dominant']

        power = wav['power']

        power[inter1[0], inter1[1]] = np.nan
        wl[inter1[0], inter1[1]] = np.nan

        scales = wav['scales']

        container['LSTA'].values = wl[None, ...]
        container['power'] = (('time', 'lat', 'lon', ), power[None, ...])
        container.drop('cell')

        container.attrs['scales'] = scales

        try:
            os.remove(outfile)
        except OSError:
            pass
        container.to_netcdf(path=outfile, mode='w')


        print('Saved ' + outfile)


