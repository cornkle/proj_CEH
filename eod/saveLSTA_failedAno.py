import xarray as xr
import pdb
import pandas as pd
import numpy as np

def saveAnomaly():
    mf = xr.open_mfdataset('/users/global/cornkle/data/OBS/MSG_LSTA/lst_new/*.nc', concat_dim='time')

    mf['ymonth'] = ('time', [str(y)+'-'+str(m) for (y,m) in zip(mf['time.year'].values,mf['time.month'].values)])

    grouped='ymonth'

    valid_days = mf.groupby(grouped).count(dim='time') # number of valid days per month

    minus =  mf.groupby(grouped).mean(dim='time')
    arr = minus['LSTA'].values

    arr[valid_days['LSTA'].values<10] = np.nan
    minus['LSTA'].values = arr

    dso = mf.groupby(grouped) - minus

    for d in dso['time']:
        try:
            arr = dso.sel(time=d.values).drop('ymonth')
        except ValueError:
            arr = dso.sel(time=d.values)
        day = arr['time.day'].values
        month = arr['time.month'].values
        year = arr['time.year'].values

        date = [pd.datetime(year, month, day, 0, 0)]
        da = xr.DataArray(arr['LSTA'].values[None, ...],
                          coords={'time': date, 'lat': arr.lat, 'lon': arr.lon},
                          dims=['time', 'lat', 'lon'])  # [np.newaxis, :]
        ds = xr.Dataset({'LSTA': da})

        date = str(arr['time.year'].values)+str(arr['time.month'].values).zfill(2)+str(arr['time.day'].values).zfill(2)
        ds.to_netcdf('/users/global/cornkle/data/OBS/MSG_LSTA/lsta_new/lsta_daily_'+date+'.nc')