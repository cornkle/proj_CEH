import numpy as np
import xarray as xr
from wavelet import util

#input = '/media/ck/Elements/global/annual_trend_data/'
input = '/prj/global_water/LST_trend_global/'

#all_data = xr.open_dataset(input+'lst_monthly_trends.nc')
all_data = xr.open_dataset(input+'lst_monthly_trends.nc')
#box = [-18,20,3,12]
all_data = all_data#.where((all_data.Longitude>box[0])&(all_data.Longitude<box[1])&(all_data.Latitude>box[2])&(all_data.Latitude<box[3]), drop=True)
months = all_data['LST_trend']

ylist = []
xlist = []
dlist = []
cnt = 0
for ts in months:
    print(cnt)
    mask = np.isfinite(ts.values)
    ts.values[np.isnan(ts.values)] = 0.01
    dic = util.waveletT1D(ts.values, dataset='LSTATREND5K_GLOBAL', mask=mask, sign='positive')
    # del ts

    powery = dic['powery']
    powerx = dic['powerx']
    powery[powery < 0.1] = 0
    powerx[powerx < 0.1] = 0

    ylist = (xr.DataArray(powery, coords={'scales': dic['scales'], 'lon': all_data.Longitude[:, 0],
                                          'lat': all_data.Latitude[:, 0]},
                          dims=['scales', 'lon', 'lat']))  # [np.newaxis, :])
    xlist = (xr.DataArray(powerx, coords={'scales': dic['scales'], 'lon': all_data.Longitude[:, 0],
                                          'lat': all_data.Latitude[:, 0]},
                          dims=['scales', 'lon', 'lat']))

    dlist = (xr.DataArray(ts.values, coords={'lon': all_data.Longitude[:, 0], 'lat': all_data.Latitude[:, 0]},
                          dims=['lon', 'lat']))
    cnt = cnt + 1
    del dic
    ds = xr.Dataset()
    #     ds['power_y'] = xr.concat(ylist, 'time')
    #     ds['power_x'] = xr.concat(xlist, 'time')
    #     ds['trend'] = xr.concat(dlist,'time')
    # ipdb.set_trace()
    ds['power_y'] = ylist
    ds['power_x'] = xlist
    ds['trend'] = dlist
    comp = dict(zlib=True, complevel=5)
    enc = {var: comp for var in ds.data_vars}

    ds.to_netcdf(
        path=input+'LSTtrend_scales_perMonth_POSITIVE_' + str(
            cnt).zfill(2) + '.nc', mode='w', encoding=enc, format='NETCDF4')
    del ds