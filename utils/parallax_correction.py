import math
import numpy as np
import pdb
import matplotlib.pyplot as plt
from utils import u_gis as ug
def parallax_corr_msg(slon, slat, plon, plat, height):
    """

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

    er = 6378.137 # earth radius equator km
    er_po = 6356.8        # earth radius pole km
    mh = 35786 #msg satellite height
    tot = er + mh
    lat_diff = np.abs(slat-plat)
    lon_diff = np.abs(slon-plon)
    ### parallax correction on the sphere
    lax_lat = (height * tot * math.sin(math.radians(lat_diff))) / (er*(tot*math.cos(math.radians(lat_diff))-(er+height)))
    lax_lon = (height * tot * math.sin(math.radians(lon_diff))) / (er * (tot * math.cos(math.radians(lon_diff)) - (er + height)))

    ### deg to km
    lax_y = er * lax_lat
    lax_x = er * lax_lon

    if plon < slon:
        lax_x = lax_x * -1
        lax_lon = lax_lon * -1

    if plat < slat:
        lax_y = lax_y * -1
        lax_lat = lax_lat -1

    return (lax_x, lax_y), (lax_lon, lax_lat)

def own_test(slon, slat, plon, plat, height):
    er = 6378.077  # earth radius equator km
    er_po = 6356.577  # earth radius pole km
    mh = 35786 #msg satellite height
    tot = er + mh
    totc = er + height
    alpha = np.abs(slon-plon)
    alpha2 = np.abs(slat-plat)

    splon = math.sqrt(tot**2+er**2-2*tot*er*math.cos(math.radians(alpha)))
    splat = math.sqrt(tot ** 2 + er ** 2 - 2 * tot * er * math.cos(math.radians(alpha2)))

    sat_angle_lon = math.degrees(math.asin(er*math.sin(math.radians(alpha))/splon))
    sat_angle_lat = math.degrees(math.asin(er * math.sin(math.radians(alpha2)) / splat))

    beta_lon = math.degrees(math.asin(tot*math.sin(math.radians(sat_angle_lon))/totc))
    beta_lon = 180-beta_lon

    beta_lat = math.degrees(math.asin(tot*math.sin(math.radians(sat_angle_lat))/totc))
    beta_lat = 180-beta_lat

    point_lon = 180 - beta_lon - sat_angle_lon
    offset_lon = alpha - point_lon

    point_lat = 180 - beta_lat - sat_angle_lat
    offset_lat = alpha2 - point_lat

    km, deg = ug.haversine(slon, slat, plon, plat)
    sboth = math.sqrt(tot ** 2 + er ** 2 - 2 * tot * er * math.cos(math.radians(deg)))
    sat_angle_both = math.degrees(math.asin(er * math.sin(math.radians(deg)) / sboth))
    beta_both = math.degrees(math.asin(tot * math.sin(math.radians(sat_angle_both)) / totc))
    beta_both = 180 - beta_both
    point_both = 180 - beta_both - sat_angle_both
    offset_both = deg - point_both

    offkm_both =  offset_both * (2 * math.pi / 360) * er


    ### to km

    offkm_lon = offset_lon  * math.cos(math.radians(plat)) * (2*math.pi/360) * er
    offkm_lat = offset_lat * (2 * math.pi / 360) * er

    return (offkm_lon, offkm_lat), (offset_lon, offset_lat)



def parallax_corr_msg_impr(slon, slat, plon, plat, height):
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
    #lon_diff = np.abs(slon-plon)
    bothkm, both_diff = ug.haversine(slon, slat, plon, plat)

    ### parallax correction on the sphere
    lax_lat = (height * tot * math.sin(math.radians(lat_diff))) / (er*(tot*math.cos(math.radians(lat_diff))-(r_local+height)))
    ##lax_lon = (height * tot * math.sin(math.radians(lon_diff))) / (er * (tot * math.cos(math.radians(lon_diff)) - (r_local + height)))
    lax_both = (height * tot * math.sin(math.radians(both_diff))) / (er * (tot * math.cos(math.radians(both_diff)) - (r_local + height)))

    lax_lat = lax_lat
    ##lax_lon = math.degrees(lax_lon)
    lax_both = lax_both
    lax_lon = np.sqrt(lax_both**2 - lax_lat**2)  ## assume trigonometric dependance cause distances are small

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



def para_check():


    x= np.arange(10,21)

    ylon = []
    ylat = []

    er = 6378.137
    circ = (2*math.pi/360)

    valslon = [6.9, 7.6, 8.3, 9,9.7, 10.4, 11.2, 11.9, 12.6, 13.3, 14]

    valslat = [14.9, 16.4, 17.9, 19.4, 20.9, 22.4, 23.9, 25.4, 26.9, 28.4, 29.9 ]

    for xx in x:

        km, degs = parallax_corr_msg_impr(-3.4, 0, 15.8, 48.5,xx)
        ylon.append(km[0])
        ylat.append(km[1])



    f = plt.figure()
    plt.plot(x, valslon)
    plt.plot(x,np.array(ylon), color='red')
    plt.title('lons')

    f = plt.figure()
    plt.plot(x, np.array(valslon)/(er*circ))
    plt.plot(x,np.array(ylon)/(er*circ), color='red')
    plt.title('lons')


    f = plt.figure()
    plt.plot(x, valslat)
    plt.plot(x, ylat, color='red')
    plt.title('lats')

