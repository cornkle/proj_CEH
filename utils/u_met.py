import numpy as np
from scipy import stats
import xarray as xr
import scipy.ndimage.interpolation as inter
from utils import constants
import pdb


def u_v_to_ws_wd(u,v):
    """
    U and V wind coords to wind speed and direction
    :param u:
    :param v:
    :return: tuple (wind speed, wind direction (convention "wind coming from")
    """

    ws = np.array(np.sqrt(u*u + v*v))
    wd = np.array(180. + np.arctan2(u, v) * 180./ np.pi)  # dir where wind is coming from
    pos = np.where(ws == 0)
    try:
        wd[pos] = np.nan
    except IndexError:
        pass

    return ws,wd


def era_Tlapse_height(month, temp, lon, lat):
    """
    Estimates the height (in metres) of a given atmospheric temperature from ERA-I monthly temperatures and
    geopotential heights on pressure levels

    :param month:
    :param temp: temperature of point
    :param lon: longitude
    :param lat: latitude
    :return: approximate height above surface
    """
    file =constants.ERA_MONTHLY_TUVWZ_AFRICA
    da = xr.open_dataset(file)
    da = da.sel(longitude=lon, latitude=lat, method='nearest')
    da = da.isel(month=month-1)
    print('Read ERA data')

    g = 9.80665
    t = da['t']-273.15
    z = da['z']
    zm = z / g  ## geopotential / gravity constant
    ismin = np.argmin(t.values)

    gradient, intercept, r_value, p_value, std_err = stats.linregress(zm[ismin:ismin+2], t[ismin:ismin+2])
    t[0:ismin + 1] = gradient * zm[0:ismin + 1] + intercept ## linear tropopause correction (after tmin, temp rises!)

    X = np.abs(t-temp)
    idx = np.argmin(X.values)
    height = zm.values[idx]

    #plt.plot(zm, t)
    # plt.plot(zm[0:ismin+1], gradient*zm[0:ismin+1]+intercept)
    return height

def era_wind_rotate(array, ptime, lat, lon, level=None, ref_angle=None):
    """

    :param array: 2d array
    :param ptime: pandas time step to choose ERA-time (daily at 12UTC)
    :param lat: latitude
    :param lon: longitude
    :param level:  ERA pressure level in hPa, note: level 0 means surface. Available: 0, 925, 850, 700, 600
    :param ref_angle: the reference angle to rotate to (direction where all wind should come from)
    :return: rotated array
    """

    if array.ndim != 2:
        raise IndexError('Cut kernel only allows 2D arrays.')

    if ref_angle==None:
        ref_angle=0

    if level==None:
        level= 0
        ustr = 'u10'
        vstr = 'v10'
    if level == 0:
        era = xr.open_dataset(constants.ERA_DAILY_SURFACE)
        point = era.sel(latitude=lat, longitude=lon, method='nearest', time=ptime)
    else:
        era = xr.open_dataset(constants.ERA_DAILY_PL)
        point = era.sel(latitude=lat, longitude=lon, method='nearest', time=ptime, level=level)
        ustr = 'u'
        vstr = 'v'

    u = point[ustr].values
    v = point[vstr].values
    ws, wd = u_v_to_ws_wd(u, v)
    rot_array = inter.rotate(array, ref_angle - wd, reshape=False, cval=np.nan, prefilter=False)

    return rot_array

def era_wind_rotate3d(array, ptime, lat, lon, level=None, ref_angle=None):
    """

    :param array: 2d array
    :param ptime: pandas time step to choose ERA-time (daily at 12UTC)
    :param lat: latitude
    :param lon: longitude
    :param level:  ERA pressure level in hPa, note: level 0 means surface. Available: 0, 925, 850, 700, 600
    :param ref_angle: the reference angle to rotate to (direction where all wind should come from)
    :return: rotated array
    """
    if array.ndim != 3:
        raise IndexError('Cut kernel only allows 3D arrays.')

    if ref_angle==None:
        ref_angle=0

    if level==None:
        level= 0
        ustr = 'u10'
        vstr = 'v10'
    if level == 0:
        era = xr.open_dataset(constants.ERA_DAILY_SURFACE)
        point = era.sel(latitude=lat, longitude=lon, method='nearest', time=ptime)
    else:
        era = xr.open_dataset(constants.ERA_DAILY_PL)
        point = era.sel(latitude=lat, longitude=lon, method='nearest', time=ptime, level=level)
        ustr = 'u'
        vstr = 'v'

    u = point[ustr].values
    v = point[vstr].values
    ws, wd = u_v_to_ws_wd(u, v)
    rot_array = np.zeros_like(array)

    for s in range(array.shape[0]):
        torot = array[s,:,:]
        rot = inter.rotate(torot, ref_angle - wd, reshape=False, cval=np.nan, prefilter=False)
        rot_array[s,:,:] = rot

    return rot_array
