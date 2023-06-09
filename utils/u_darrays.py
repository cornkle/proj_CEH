import numpy as np
import xarray as xr
from utils import u_mann_kendall as mk, u_arrays as ua
from scipy import stats
import matplotlib.pyplot as plt
import bottleneck
import scipy
import ipdb
import salem



def coarsen(data, factor):
    grid = data.salem.grid.regrid(factor=factor)
    coarse = grid.lookup_transform(data)
    grid = grid.to_dataset()
    da = xr.DataArray(coarse, coords=[data['time'], grid['y'], grid['x']], dims=['time', 'lat', 'lon'])
    return da


def shift_lons(ds, lon_dim='lon', save=False):
    """ Shift longitudes from [0, 360] to [-180, 180] """
    lons = ds[lon_dim].values
    new_lons = np.empty_like(lons)
    mask = lons > 180
    new_lons[mask] = -(360 - lons[mask])
    new_lons[~mask] = lons[~mask]
    new_lons = np.sort(new_lons)
    #ds[lon_dim].values = new_lons
    #ipdb.set_trace()
    #ipdb.set_trace()
    ds = ds.assign_coords({lon_dim:new_lons})
   # ipdb.set_trace()s
    if save:
        ds.to_netcdf(save)
    return ds


def roll_lons(ds, lon_dim='lon'):
    lon_name = lon_dim  # whatever name is in the data

    # Adjust lon values to make sure they are within (-180, 180)
    ds['_longitude_adjusted'] = xr.where(
        ds[lon_name] > 180,
        ds[lon_name] - 360,
        ds[lon_name])

    # reassign the new coords to as the main lon coords
    # and sort DataArray using new coordinate values
    ds = (
        ds
        .swap_dims({lon_name: '_longitude_adjusted'})
        .sel(**{'_longitude_adjusted': sorted(ds._longitude_adjusted)})
        .drop(lon_name))

    ds = ds.rename({'_longitude_adjusted': lon_name})

    return ds




def shift_lons_data(ds, lon_dim='lon', save=False):
    """ Shift longitudes from [0, 360] to [-180, 180] """
    lons = ds[lon_dim].values

    mask = np.where(lons > 180)
    mask2 = np.where((lons < 180)  & (lons != 0))
    #ipdb.set_trace()
    argpos = np.argmin(np.abs(ds[lon_dim].values-180))
    new_lons = lons - ds[lon_dim].values[argpos]#179.0625
    #ipdb.set_trace()
    for dv in ds.data_vars:
        new_data = np.empty_like(ds[dv].values)
        ndims = new_data.ndim
        if ndims == 2:
            try:
                new_data[:,mask2[0]] = ds[dv].values[:,mask[0]]
                new_data[:,mask[0]] = ds[dv].values[:,mask2[0]]
            except:
                print('2dim prob')
                ipdb.set_trace()
        if ndims == 3:
            try:
                new_data[:,:,mask2[0]] = ds[dv].values[:,:,mask[0]]
                new_data[:,:,mask[0]] = ds[dv].values[:,:,mask2[0]]
            except:
                print('3dim prob')
                ipdb.set_trace()
        if ndims == 4:
            try:
                new_data[:,:, :, mask2[0]] = ds[dv].values[:, :, :, mask[0]]
                new_data[:, : , :, mask[0]] = ds[dv].values[:,:, :, mask2[0]]
            except:
                print('4dim prob')
                ipdb.set_trace()

        ds[dv] = ds[dv].assign_coords({lon_dim: new_lons})
        ds[dv].values = new_data
    #ipdb.set_trace()
    if save:
        ds.to_netcdf(save)
    return ds

def linear_trend_mk(x, eps=0.001, alpha=0.01, nb_valid=None):


    xin = x.where(np.isfinite(x), drop=True).squeeze()

    min_val = 3
    if nb_valid is not None:
        min_val = nb_valid

    if np.sum(np.isfinite(xin)) < min_val:

        pf, slope, int, p, ind = mk.test(np.arange(len(x)), x.squeeze().values, eps=eps, alpha=alpha, Ha='upordown')

        ds = xr.Dataset()
        ds['slope'] = xr.DataArray(slope*np.nan, )
        ds['pval'] = xr.DataArray(p*np.nan, )
        ds['ind'] = xr.DataArray(ind*np.nan, )

        # try:
        #     ds['allpoints']
        # except KeyError:
        #     ipdb.set_trace()

        return ds

    pf, slope, int, p, ind = mk.test(np.arange(len(xin)),xin.values, eps=eps, alpha=alpha, Ha='upordown')

    # we need to return a dataarray or else xarray's groupby won't be happy

    ds = xr.Dataset()
    ds['slope'] = xr.DataArray(slope,)
    ds['pval'] = xr.DataArray(p, )
    ds['ind'] = xr.DataArray(ind,)

    # try:
    #     ds['allpoints']
    # except KeyError:
    #     ipdb.set_trace()

    return ds

def linear_trend_lingress(x, nb_valid=None):

    xin = x.where(np.isfinite(x), drop=True).squeeze()

    min_val = 3
    if nb_valid is not None:
        min_val = nb_valid

    if np.sum(np.isfinite(xin)) < min_val:
        ds = xr.Dataset()
        ds['slope'] = xr.DataArray(np.nan, )
        ds['pval'] = xr.DataArray(np.nan, )
        ds['r'] = xr.DataArray(np.nan, )

        # try:
        #     ds['allpoints']
        # except KeyError:
        #     ipdb.set_trace()

        return ds

    slope, intercept, r, p, std_err = stats.linregress(np.arange(len(xin)), xin)

    # we need to return a dataarray or else xarray's groupby won't be happy

    ds = xr.Dataset()
    ds['slope'] = xr.DataArray(slope,)
    ds['pval'] = xr.DataArray(p, )
    ds['r'] = xr.DataArray(r,)

    # try:
    #     ds['allpoints']
    # except KeyError:
    #     ipdb.set_trace()

    return ds



def linear_trend_polyfit(x, nb_valid=None):

    xin = x.where(np.isfinite(x), drop=True).squeeze()

    min_val = 3
    if nb_valid is not None:
        min_val = nb_valid

    if np.sum(np.isfinite(xin)) < min_val:
        ds = xr.Dataset()
        ds['slope'] = xr.DataArray(np.nan, )
        ds['pval'] = xr.DataArray(np.nan, )
        ds['r'] = xr.DataArray(np.nan, )

        # try:
        #     ds['allpoints']
        # except KeyError:
        #     ipdb.set_trace()

        return ds

    slope, inter = np.polyfit(np.arange(len(xin)), xin,1)

    # we need to return a dataarray or else xarray's groupby won't be happy

    ds = xr.Dataset()
    ds['slope'] = xr.DataArray(slope,)

    # try:
    #     ds['allpoints']
    # except KeyError:
    #     ipdb.set_trace()

    return ds


def flip_lat_write(filepath):

    ds = xr.open_dataset(filepath)
    for var in ds.data_vars:
        ds[var].values = ds[var].values[:,::-1,:]
    comp = dict(zlib=True, complevel=5)
    enc = {var: comp for var in ds.data_vars}
    savefile = filepath.replace('.nc', '_yflip.nc')
    ds.to_netcdf(path=savefile, mode='w', encoding=enc, format='NETCDF4')
    print('Saved ' + savefile)



def flip_lat(ds):
    """
    Flips latitude of dataset. Only works with correct latitude name...
    :param ds:
    :return:
    """
    try:
        ds = ds.sel(latitude=slice(None, None, -1))
    except ValueError:
        try:
            ds = ds.sel(lat=slice(None, None, -1))
        except ValueError:
            try:
                ds = ds.sel(Latitude=slice(None, None, -1))
            except ValueError:
                try:
                    ds = ds.sel(allpoints_level_0=slice(None, None, -1))
                except ValueError:
                    print('Cant flip lat, coord name not found')

    return ds


def to_newarray(da):

    das = xr.DataArray(da.values,
                      coords={'time': da.time, 'latitude': da.latitude.values,
                              'longitude': da.longitude.values},
                      dims=['time', 'latitude', 'longitude'])  # [np.newaxis, :]
    das.attrs = da.attrs
    return das



################### correlation computation with quick parallel dask usage
def covariance_gufunc(x, y):
    return ((x - x.mean(axis=-1, keepdims=True))
            * (y - y.mean(axis=-1, keepdims=True))).mean(axis=-1)

def pearson_correlation_gufunc(x, y):
    return covariance_gufunc(x, y) / (x.std(axis=-1) * y.std(axis=-1))

def spearman_correlation_gufunc(x, y):
    x_ranks = bottleneck.rankdata(x, axis=-1)
    y_ranks = bottleneck.rankdata(y, axis=-1)
    return pearson_correlation_gufunc(x_ranks, y_ranks)

def spearman_correlation(x, y, dim):
    return xr.apply_ufunc(
        spearman_correlation_gufunc, x, y,
        input_core_dims=[[dim], [dim]],
        dask='parallelized',
        output_dtypes=[float])

def pearson_correlation(x, y, dim):
    return xr.apply_ufunc(
        pearson_correlation_gufunc, x, y,
        input_core_dims=[[dim], [dim]],
        dask='parallelized',
        output_dtypes=[float])


def cut_box(arr, xpos=None, ypos=None, dist=None):
    """

    :param xpos: x coordinate in domain for kernel centre point
    :param ypos: y coordinate in domain for kernel centre point
    :param arr: numpy array (2d)
    :param dist: distance from kernel centre point to kernel edge (total width = 2*dist+1)
    :return: the kernel of dimensions (2*dist+1, 2*dist+1)
    """

    if dist == None:
        'Distance missing. Please provide distance from kernel centre to edge (number of pixels).'
        return

    if arr.ndim == 2:
        kernel = ua.cut_kernel(arr.values,xpos, ypos,dist)
        if kernel.shape != (dist * 2 + 1, dist * 2 + 1):
            print("Please check kernel dimensions, there is something wrong")
            ipdb.set_trace()
    elif arr.ndim == 3:
        kernel = ua.cut_kernel_3d(arr.values,xpos, ypos,dist)

        if kernel.shape != (arr.shape[0], dist * 2 + 1, dist * 2 + 1):
            print("Please check kernel dimensions, there is something wrong")
            ipdb.set_trace()
    else:
        print('Dimension of array not supported, please check')
        ipdb.set_trace()

    if arr.ndim == 3:
        try:
            levels = arr.level.values
        except AttributeError:
            levels = arr.pressure.values

        return xr.DataArray(kernel, dims=['level','y','x'],
                            coords={'level' : arr.level.values})
    else:
        return xr.DataArray(kernel, dims=['y','x'])



def box_correlation(a, b, bsingle=None, c_box=None):
    ds = xr.Dataset()
    ds['pval'] = a.copy(deep=True).sum('year') * np.nan
    ds['r'] = a.copy(deep=True).sum('year') * np.nan
    ds['slope'] = a.copy(deep=True).sum('year') * np.nan

    corr_box = c_box

    if bsingle:
        bb = b
    else:
        bb = b.sel(latitude=slice(corr_box[2], corr_box[3]), longitude=slice(corr_box[0], corr_box[1])).mean(dim=['latitude', 'longitude'])

    for lat in a.latitude.values:
        for lon in a.longitude.values:
            aa = a.sel(latitude=lat, longitude=lon)
            if bsingle:
                r, p = stats.pearsonr(aa.values, bb)

                pf = np.polyfit(aa.values, bb, 1)
            else:
                r, p = stats.pearsonr(aa.values, bb.values)
                pf = np.polyfit(aa.values, bb.values, 1)


            slope = pf[0]

            if (np.nansum(aa.values == 0) >= 10):
                p = np.nan
                r = np.nan

            ds['r'].loc[{'latitude': lat, 'longitude': lon}] = r
            ds['pval'].loc[{'latitude': lat, 'longitude': lon}] = p
            ds['slope'].loc[{'latitude': lat, 'longitude': lon}] = slope

    return ds



def leap_year(year, calendar='standard'):
    """Determine if year is a leap year"""
    leap = False
    if ((calendar in ['standard', 'gregorian',
        'proleptic_gregorian', 'julian']) and
        (year % 4 == 0)):
        leap = True
        if ((calendar == 'proleptic_gregorian') and
            (year % 100 == 0) and
            (year % 400 != 0)):
            leap = False
        elif ((calendar in ['standard', 'gregorian']) and
                 (year % 100 == 0) and (year % 400 != 0) and
                 (year < 1583)):
            leap = False
    return leap


dpm = {'noleap': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
           '365_day': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
           'standard': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
           'gregorian': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
           'proleptic_gregorian': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
           'all_leap': [0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
           '366_day': [0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
           '360_day': [0, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30]}
def get_dpm(time, calendar='standard'):
        """
        return a array of days per month corresponding to the months provided in `months`
        """
        month_length = np.zeros(len(time), dtype=np.int)

        cal_days = dpm[calendar]

        for i, (month, year) in enumerate(zip(time.month, time.year)):
            month_length[i] = cal_days[month]
            if leap_year(year, calendar=calendar):
                month_length[i] += 1
        return month_length




def season_mean(ds, calendar='standard', mstring='QS-DEC', mpick=12):
    # Make a DataArray of season/year groups
    #     year_season = xr.DataArray(ds.time.to_index().to_period(freq='Q-NOV').to_timestamp(how='E'),
    #                                coords=[ds.time], name='year_season')

    month_length = xr.DataArray(get_dpm(ds.time.to_index(), calendar=calendar), coords=[ds.time],
                                name='month_length')
    result = ((ds * month_length).resample(time=mstring).sum() / month_length.where(ds.notnull()).resample(
        time=mstring).sum())  # QS-DEC

    #ipdb.set_trace()

    # Calculate the weighted average
    # ipdb.set_trace()
    return result.isel(time=result['time.month'] == 12)

