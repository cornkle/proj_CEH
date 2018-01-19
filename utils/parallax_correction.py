import math
import numpy as np
import pdb

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

def parallax_corr_msg_impr(slon, slat, plon, plat, height):
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
    ### parallax correction on the sphere
    lax_lat = (height * tot * math.sin(math.radians(lat_diff))) / (er*(tot*math.cos(math.radians(lat_diff))-(r_local+height)))
    lax_lon = (height * tot * math.sin(math.radians(lon_diff))) / (er * (tot * math.cos(math.radians(lon_diff)) - (r_local + height)))

    ### deg to km
    lax_y = r_local * lax_lat
    lax_x = r_local * lax_lon

    lax_lat = math.degrees(lax_lat)
    lax_lon = math.degrees(lax_lon)

    if plon < slon:
        lax_x = lax_x * -1
        lax_lon = lax_lon * -1

    if plat < slat:
        lax_y = lax_y * -1
        lax_lat = lax_lat -1

    return (lax_x, lax_y), (lax_lon, lax_lat)


def parallax_corr_msgEU(slon, slat, plon, plat, height):
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
    r_mean = 0.5*(er+er_po) # mean radius
    r_ratio = er / er_po
    # geodetic latitude of satellite
    slat_g = math.degrees(math.atan(math.tan(math.radians(slat)) * r_ratio**2))
    # geodetic latitude of point
    plat_g = math.degrees(math.atan(math.tan(math.radians(plat)) * r_ratio**2))
    r_local = er / math.sqrt(math.cos(math.radians(plat_g))**2 + r_ratio**2*math.sin(math.radians(plat_g))**2) #local radius at cloud point
    r_ratio_local = ((er + height) / (er_po + height))**2

    mh = 35786 + er #msg satellite height

    # satellite position cartesian coords
    xs = mh * math.cos(math.radians(slat_g))*math.sin(math.radians(slon))
    ys = mh * math.sin(math.radians(slon))
    zs = mh * math.cos(math.radians(slat_g))*math.cos(math.radians(slon))

    # point position cartesian coords
    xp = r_local*math.cos(math.radians(plat_g))*math.sin(math.radians(plon))
    yp = r_local*math.sin(math.radians(plon))
    zp = r_local*math.cos(math.radians(plat_g))*math.cos(math.radians(plon))

    # difference vector between satellite and cloud
    xdiff = xs - xp
    ydiff = ys - yp
    zdiff = zs - zp

    # Correction for the line of sight at cloud top height
    e1 = xdiff**2 + r_ratio_local * ydiff**2 + zdiff**2
    e2 = 2*(xp*xdiff+r_ratio_local*yp*ydiff+zp*zdiff)
    e3 = xp**2+r_ratio_local*yp**2+zp**2-(er+height)**2

    c = (math.sqrt(e2**2-4*e1*e3)-e2) / 2*e1

    xoff = c*xdiff
    yoff = c*ydiff
    zoff = c*zdiff

    latoff = math.degrees(math.atan(math.tan(math.tan(math.radians(yoff/math.sqrt(xoff**2+zoff**2))))/ r_ratio ** 2))
    lonoff = math.degrees(math.atan2(xoff,yoff))

    # apply correction to cartesian coordinates of cloud position
    xcorr = xp + c*xdiff
    ycorr = yp + c*ydiff
    zcorr = zp + c*zdiff

    latcorr = math.degrees(math.atan(math.tan(math.tan(math.radians(ycorr/math.sqrt(xcorr**2+zcorr**2))))/ r_ratio ** 2))
    loncorr = math.degrees(math.atan2(xcorr,ycorr))


    pdb.set_trace()


    return
