import numpy as np

HOD=range(24)   # hours of day
CFAC=-781648343
LFAC=-781648343
COFF=1856
LOFF=1856


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
