import numpy as np
from scipy import stats
import xarray as xr
import scipy.ndimage.interpolation as inter
from utils import constants
from metpy import calc as metcalc


Lv = 2.501e6 # heat of vapourisation
Cp = 1005 # heat capacity of dry air at static pressure
g = 9.80665


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


def olr_to_bt(olr):
    sigma = 5.670373e-8
    return ((olr/sigma)**0.25)-273.15


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




def saturation_equivalent_potential_temperature(pressure, temperature):
    r"""Calculate saturation equivalent potential temperature.

    This calculation must be given an air parcel's pressure and temperature.
    The implementation uses the formula outlined in [Bolton1980]_ for the
    equivalent potential temperature, and assumes a saturated process.

    First, because we assume a saturated process, the temperature at the LCL is
    equivalent to the current temperature. Therefore the following equation

    .. math:: T_{L}=\frac{1}{\frac{1}{T_{D}-56}+\frac{ln(T_{K}/T_{D})}{800}}+56

    reduces to

    .. math:: T_{L} = T_{K}

    Then the potential temperature at the temperature/LCL is calculated:

    .. math:: \theta_{DL}=T_{K}\left(\frac{1000}{p-e}\right)^k
              \left(\frac{T_{K}}{T_{L}}\right)^{.28r}

    However, because

    .. math:: T_{L} = T_{K}

    it follows that

    .. math:: \theta_{DL}=T_{K}\left(\frac{1000}{p-e}\right)^k

    Both of these are used to calculate the final equivalent potential temperature:

    .. math:: \theta_{E}=\theta_{DL}\exp\left[\left(\frac{3036.}{T_{K}}
                                              -1.78\right)*r(1+.448r)\right]

    Parameters
    ----------
    pressure: `pint.Quantity`
        Total atmospheric pressure
    temperature: `pint.Quantity`
        Temperature of parcel

    Returns
    -------
    `pint.Quantity`
        The saturation equivalent potential temperature of the parcel

    Notes
    -----
    [Bolton1980]_ formula for Theta-e is used (for saturated case), since according to
    [DaviesJones2009]_ it is the most accurate non-iterative formulation
    available.

    """
    t = temperature.to('kelvin').magnitude
    p = pressure.to('hPa').magnitude
    e = saturation_vapor_pressure(temperature).to('hPa').magnitude
    r = saturation_mixing_ratio(pressure, temperature).magnitude

    th_l = t * (1000 / (p - e)) ** mpconsts.kappa
    th_es = th_l * np.exp((3036. / t - 1.78) * r * (1 + 0.448 * r))

    return th_es * units.kelvin


