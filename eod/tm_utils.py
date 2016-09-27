import numpy as np
import math

HOD=range(24)   # hours of day
CFAC=-781648343
LFAC=-781648343
COFF=1856
LOFF=1856


# msg grid: 1856 x 1856 pixels, index start south eastern corner, x rises from right to left!
# ==============================================================================
# To MSG indices
# ==============================================================================
def ll_toMSG(lons, lats, cfac=CFAC, lfac=LFAC, coff=COFF, loff=LOFF):
    lats = np.array(lats)
    lons = np.array(lons)

    if not lats.shape == lons.shape:
        print('Lats lons must have same dimensions')
        return
    if not lats.size == lons.size:
        print('Lats lons must have same size')
        return

    if (np.min(lats) < -90.) | (np.max(lats) > 90.):
        print('Lats are out of range')
        return
    if (np.min(lons) < -180.) | (np.max(lons) > 180.):
        print('Lons are out of range')
        return

    pi = 3.14159265359  # Define as double precision.
    sat_height = 42164.0  # Height of satellite.
    r_eq = 6378.169  # Radius of Earth at equator.
    r_pol = 6356.5838  # Radius of Earth at pole.
    sub_lon = 0.0  # Longitude of sub-satellite point.

    # Convert lats and longs to radians.
    lats_r = lats * pi / 180.
    lons_r = lons * pi / 180.

    # Calculate geocentric latitude from the geographic one.
    c_lat = np.arctan(0.993243 * np.sin(lats_r) / np.cos(lats_r))

    # Use c_lat to calculate the length from the Earth centre
    # to the surface of the Earth ellipsoid equations.
    re = r_pol / np.sqrt(1. - 0.00675701 * np.cos(c_lat) * np.cos(c_lat))

    # Calculate the forward projection.
    r1 = sat_height - re * np.cos(c_lat) * np.cos(lons_r - sub_lon)
    r2 = -re * np.cos(c_lat) * np.sin(lons_r - sub_lon)
    r3 = re * np.sin(c_lat)
    rn = np.sqrt(r1 * r1 + r2 * r2 + r3 * r3)

    # Create output arrays.
    cols = np.empty_like(lats)
    rows = np.empty_like(lats)

    # Check for visibility, whether the point is visible from the satellite.
    dotprod = np.array([r1 * re * np.cos(c_lat) * np.cos(lons_r - sub_lon) - r2 * r2 - r3 * r3 * (r_eq / r_pol) ** 2.])

    cols = np.arctan(-r2 / r1)
    rows = np.arcsin(-r3 / rn)

    cols = np.round(coff + cols * cfac / 2. ** 16)  # This seems to be incorrect in the example program.
    rows = np.round(loff + rows * lfac / 2. ** 16)

    dic = {'x': cols, 'y': rows}

    return dic


def ll_toMSG_rev(lon, lat, rflon=0):
    #    """
    #    input
    #    vlon: longitude
    #    vlat: latitude
    #    rflon: satellite longitude, default=0degrees
    #    output xr: x position
    #    yr: y position
    #    """
    # geo2file_geos(lon, lat)

    # Setup constants
    re = 6378.160  # equatorial radius
    h = 42164.0 - re  # Reference altitude
    rp = 6356.5838
    lpsi2 = 1  # spin direction
    resol = 3712.0  # pixel/line number for MSG
    deltax = 17.83 / resol  # E-W scanning step   scan amplitude is 17.83 deg. for MSG
    deltay = 17.83 / resol  #

    dtor = math.radians(1.0)
    xlat = dtor * lat
    xlon = dtor * lon
    cosxlat = math.cos(xlat)
    sinxlat = math.sin(xlat)
    cosxlon = math.cos(xlon)

    rom = (re * rp) / math.sqrt(rp * rp * cosxlat * cosxlat + re * re * sinxlat * sinxlat)
    y = math.sqrt(h * h + rom * rom - 2.0 * h * rom * cosxlat * cosxlon)
    r1 = y * y + rom * rom
    r2 = h * h

    if (r1 > r2):
        dic = {'x': -1, 'y': -1}
        return dic

    rs = re + h
    reph = re
    rpph = rp
    coslo = math.cos(rflon * dtor)
    sinlo = math.sin(rflon * dtor)
    teta = math.atan((rpph / reph) * math.tan(xlat))
    xt = reph * math.cos(teta) * cosxlon
    yt = reph * math.cos(teta) * math.sin(xlon)
    zt = rpph * math.sin(teta)
    px = math.atan((coslo * (yt - rs * sinlo) - (xt - rs * coslo) * sinlo) / (
    sinlo * (yt - rs * sinlo) + (xt - rs * coslo) * coslo))
    py = math.atan(zt * ((math.tan(px) * sinlo - coslo) / (xt - rs * coslo)) * math.cos(px))
    px = px / dtor
    py = py / dtor
    xr = px / (deltax * lpsi2)
    yr = py / (deltay * lpsi2)
    if (xr >= 0.):
        xr = int(px / (deltax * lpsi2)) + 0.5
    else:
        xr = int(px / (deltax * lpsi2)) - 0.5
    if (yr >= 0.):
        yr = int(py / (deltax * lpsi2)) + 0.5
    else:
        yr = int(py / (deltax * lpsi2)) - 0.5
    xr = xr + 0.5 * resol + 0.5
    yr = yr + 0.5 * resol + 0.5

    dic = {'x': xr, 'y': yr}

    return dic


def getTRMMconv(flags):
    bla = flags.astype(int)
    npfalse = []

    for b, i in zip(np.nditer(bla), range(bla.size)):
        bb = '{0:016b}'.format(int(b))
        npfalse.append(int(bb[-6]))

    return npfalse


def getTRMMstrat(flags):
    bla = flags.astype(int)
    npfalse = []

    for b, i in zip(np.nditer(bla), range(bla.size)):
        bb = '{0:016b}'.format(int(b))
        npfalse.append(int(bb[-5]))

    return npfalse

"""
Returns the smallest minute difference of arbitrary minutes with respect
to the usual MSG interval of 15 or 30 minutes
"""
def minute_delta(tmin, msg_interval):

    if msg_interval != 15 and msg_interval != 30:
        print('Not allowed msg interval')
        exit()

    if msg_interval == 15:
        arr = np.array([15,30,45,60,0])
    if msg_interval == 30:
        arr = np.array([30,60,0])

    dm = arr - tmin
    ind = (np.abs(dm)).argmin()

    dt0 = dm[ind]

    return dt0

"""create one unique integer from two positive integers
 Cantor pairing function"""
def unique_of_pair(x,y):

    uni = (x + y) * (x + y + 1) / 2 + y
    return uni

"""search for non zero value in +-1 kernel"""
def kernel_no_zero(arr, xx, yy):

    shape = arr.shape
    a = 1
    b = 2
    c = 1
    d =2

    if arr[yy,xx]:
        return arr[yy, xx]

    if (yy==0): a = 0
    if (yy==shape[0]): b=0
    if (xx==0): c=0
    if (xx == shape[1]): d = 0

    if arr[yy - a:yy + b, xx - c:xx + d].any():
        sub = arr[yy - a:yy + b, xx - c:xx + d]
        nb = sub[sub > 0][0]
    else:
        nb = False

    return nb

"""Cuts out a kernel box surrounding the index xx, yy with pixel distance r.
    r is at least 1
"""
def cut_kernel(arr, xx, yy, r):

    return arr[yy - r:yy + r+1, xx - r:xx + r+1]