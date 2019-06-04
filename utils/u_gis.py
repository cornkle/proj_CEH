
import numpy as np
from salem import get_demo_file, open_xr_dataset, GeoTiff, wgs84
from utils import u_met
from utils import u_gis
import pdb
import math
from math import radians, cos, sin, asin, sqrt

def read_tiff(file):

    g = GeoTiff(file)
    ls = g.get_vardata()

    return ls

def grid_tiff(file, grid):

    g = GeoTiff(file)
    # Spare memory
    ex = grid.extent_in_crs(crs=wgs84)  # l, r, b, t
    g.set_subset(corners=((ex[0], ex[2]), (ex[1], ex[3])),
                 crs=wgs84, margin=10)
    ls = g.get_vardata()
    ls_on_grid = grid.lookup_transform(ls, grid=grid)

    return ls_on_grid

def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a))
    r = 6371 # Radius of earth in kilometers. Use 3956 for miles

    ## gives back distance in km, degrees
    return c * r, math.degrees(c)


def parallax_corr_msg_old(slon, slat, plon, plat, height):

    er = 6378.137 # earth radius km
    mh = 35786 #msg satellite height
    tot = er + mh
    lat_diff = np.abs(slat-plat)
    lon_diff = np.abs(slon-plon)
    lax_lat = (height * tot * math.sin(math.radians(lat_diff))) / (er*(tot*math.cos(math.radians(lat_diff))-(er+height)))
    lax_lon = (height * tot * math.sin(math.radians(lon_diff))) / (er * (tot * math.cos(math.radians(lon_diff)) - (er + height)))

    lax_y = er * lax_lat
    lax_x = er * lax_lon

    #print('Degree_lat_lon' , lax_lat, lax_lon)
    #print('Km_y_x', lax_y, lax_x)

    if plon < slon:
        lax_x = lax_x * -1

    if plat < slat:
        lax_y = lax_y * -1

    return (lax_x, lax_y), (math.degrees(lax_lon), math.degrees(lax_lat))

def parallax_corr_msg(slon, slat, plon, plat, height):
    """(
    :param slon: Satellite longitude
    :param slat: Satellite latitude
    :param plon: Cloud point longitude
    :param plat: Cloud point latitude
    :param height: Cloud top height
    :return: Tuple, (parallax in x direction (km), parallax in y direction (km)), (parallax in lon direction (deg), parallax in lat direction (deg))
             Note: Parallax is the absolute distance between the point the satellite "assumes" to see and the actual location
             of the cloud. For coordinate correction, the parallax of cloud points to the West and South of the satellite are
             defined to be negative here, such that the correction follows as: cloud point - parallax = corrected location.
    """
    er = 6378.077 # earth radius equator km
    er_po = 6356.577        # earth radius pole km
    mh = 35786 #msg satellite height
    tot = er + mh
    r_ratio = er / er_po
    # geodetic latitude of satellite
    slat_g = math.degrees(math.atan(math.tan(math.radians(slat)) * r_ratio ** 2))
    # geodetic latitude of point
    plat_g = math.degrees(math.atan(math.tan(math.radians(plat)) * r_ratio ** 2))
    r_local = er / math.sqrt(math.cos(math.radians(plat_g)) ** 2 + r_ratio ** 2 * math.sin(
        math.radians(plat_g)) ** 2)  # local radius at cloud point

    lat_diff = np.abs(slat_g-plat_g)
    lon_diff = np.abs(slon-plon)
    bothkm, both_diff = haversine(slon, slat, plon, plat)

    ### parallax correction on the sphere
    lax_lat = (height * tot * math.sin(math.radians(lat_diff))) / (er*(tot*math.cos(math.radians(lat_diff))-(r_local+height)))
    lax_lon_single = (height * tot * math.sin(math.radians(lon_diff))) / (er * (tot * math.cos(math.radians(lon_diff)) - (r_local + height)))
    lax_both = (height * tot * math.sin(math.radians(both_diff))) / (er * (tot * math.cos(math.radians(both_diff)) - (r_local + height)))

    if lax_both > lax_lat:
        lax_lon = np.sqrt(lax_both**2 - lax_lat**2)  ## assume trigonometric dependance cause distances are small
    else:
        lax_lon = lax_lon_single

    if plat < slat:
        lax_lat = lax_lat - 1

    ### deg to km
    lax_y = r_local * lax_lat
    lax_x = r_local * lax_lon * math.cos(math.radians(lat_diff))

    if plon < slon:
        lax_x = lax_x * -1
        lax_lon = lax_lon * -1
    if plat < slat:
        lax_y = lax_y * -1

    return (lax_x, lax_y), (math.degrees(lax_lon), math.degrees(lax_lat))


def call_parallax_era(month, t_cloud, lon_cloud, lat_cloud, lon_sat, lat_sat):

    height = u_met.era_Tlapse_height(month, t_cloud, lon_cloud, lat_cloud)  # height in meters
    km, coords = u_gis.parallax_corr_msg(lon_sat, lat_sat, lon_cloud, lat_cloud, height / 1000)

    return km, coords


def parallax_correction(slon, slat, plon, plat, height, sheight):
    """(
    :param slon: Satellite longitude
    :param slat: Satellite latitude
    :param plon: Cloud point longitude
    :param plat: Cloud point latitude
    :param height: Cloud top height
    :param sheight: Satellite height
    :return: Tuple, (parallax in x direction (km), parallax in y direction (km)), (parallax in lon direction (deg), parallax in lat direction (deg))
             Note: Parallax is the absolute distance between the point the satellite "assumes" to see and the actual location
             of the cloud. For coordinate correction, the parallax of cloud points to the West and South of the satellite are
             defined to be negative here, such that the correction follows as: cloud point - parallax = corrected location.
    """
    er = 6378.077 # earth radius equator km
    er_po = 6356.577        # earth radius pole km
    mh = sheight # satellite height
    tot = er + mh
    r_ratio = er / er_po
    # geodetic latitude of satellite
    slat_g = math.degrees(math.atan(math.tan(math.radians(slat)) * r_ratio ** 2))
    # geodetic latitude of point
    plat_g = math.degrees(math.atan(math.tan(math.radians(plat)) * r_ratio ** 2))
    r_local = er / math.sqrt(math.cos(math.radians(plat_g)) ** 2 + r_ratio ** 2 * math.sin(
        math.radians(plat_g)) ** 2)  # local radius at cloud point

    lat_diff = np.abs(slat_g-plat_g)
    lon_diff = np.abs(slon-plon)
    bothkm, both_diff = haversine(slon, slat, plon, plat)

    ### parallax correction on the sphere
    lax_lat = (height * tot * math.sin(math.radians(lat_diff))) / (er*(tot*math.cos(math.radians(lat_diff))-(r_local+height)))
    lax_lon_single = (height * tot * math.sin(math.radians(lon_diff))) / (er * (tot * math.cos(math.radians(lon_diff)) - (r_local + height)))
    lax_both = (height * tot * math.sin(math.radians(both_diff))) / (er * (tot * math.cos(math.radians(both_diff)) - (r_local + height)))

    if lax_both > lax_lat:
        lax_lon = np.sqrt(lax_both**2 - lax_lat**2)  ## assume trigonometric dependance cause distances are small
    else:
        lax_lon = lax_lon_single

    if plat < slat:
        lax_lat = lax_lat - 1

    ### deg to km
    lax_y = r_local * lax_lat
    lax_x = r_local * lax_lon * math.cos(math.radians(lat_diff))

    if plon < slon:
        lax_x = lax_x * -1
        lax_lon = lax_lon * -1
    if plat < slat:
        lax_y = lax_y * -1

    return (lax_x, lax_y), (math.degrees(lax_lon), math.degrees(lax_lat))





