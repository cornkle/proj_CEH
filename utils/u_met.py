import numpy as np
from scipy import stats
import xarray as xr
import scipy.ndimage.interpolation as inter
from utils import constants
import metpy
from metpy import calc
from metpy.units import units
import ipdb

Lv = 2.501e6 # heat of vapourisation
Cp = 1005 # heat capacity of dry air at static pressure
g = 9.80665
Rv =462 #[J kg-1 K-1]


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


def theta_factor(pz):

    return (1000 /pz) ** 0.286


def theta(pz, t):

    if np.max(t) >100:
        print('t is very big, please check, should be given in C')

    ist = t + 273.15
    try:
        ist =  ist * ((1000 / pz) ** 0.286)  # kappa = 0.286
    except ValueError:
        ist = (ist.T * ((1000 / pz) ** 0.286)).T

    return  ist-273.15


def olr_to_bt_uncorrected(olr):
    sigma = 5.670373e-8
    return ((olr/sigma)**0.25)-273.15


def olr_to_bt(olr):
    # after Yang and Slingo (2001)
    a = 1.228
    b = -1.106e-3
    sigma = 5.67e-8
    tf = (olr/sigma)**0.25
    tb = (-a + np.sqrt(a**2 + 4*b*tf))/(2*b)
    return tb



def moist_static_energy(T, q, geop=None, z=None):

    if geop is None:
        try:
            geop = z * g
        except:
            'Please provide z (height above ground) for calculation of geopotential'

    mse = (T+273.15) * Cp  + geop + q * Lv
    return mse  #J/kg

def dry_static_energy(T, geop=None, z=None):

    if geop is None:
        try:
            geop = z * g
        except:
            'Please provide z (height above ground) for calculation of geopotential'

    mse = (T+273.15) * Cp  + geop
    return mse


def theta_e(pz, t, q):

    if np.max(q) >10:
        print('q is very big, please check, should be given in kg/kg')

    if np.max(t) >100:
        print('t is very big, please check, should be given in C')

    ist = t + 273.15
    try:
        ist =  (ist + Lv/Cp*q) * ((1000 / pz) ** 0.286)  # kappa = 0.286
    except ValueError:
        ist = ((ist.T + Lv/Cp*q.T)* ((1000 / pz) ** 0.286)).T  # (T + Lv/cpd * mixRat)(p0/p)**kappa, Stull 1988

    return  ist-273.15


def qdeficit(q1,t2, p2):

    pup = units.Quantity(p2, 'hPa')
    tup = units.Quantity(t2, 'K')

    thetae_up = np.array(calc.saturation_equivalent_potential_temperature(pup, tup))
    theta_up = theta(p2, t2 - 273.15) + 273.15

    qsat_inK = thetae_up - theta_up

    try:
        qsat =  ((qsat_inK * Cp / Lv)/ ((1000 / p2) ** 0.286))  # kappa = 0.286
    except ValueError:
        qsat = ((qsat_inK.T * Cp / Lv)/ ((1000 / p2) ** 0.286)).T  # (T + Lv/cpd * mixRat)(p0/p)**kappa, Stull 1988

    qdef = qsat - q1

    return qdef


