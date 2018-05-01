import numpy as np
import os
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs


#========================================================================================
# Rewrites msg lat lon to something nicer
#  file: lat lon grads file
#  ny : pixel in y direction
#  nx : pixel in x direction
#========================================================================================
def rewriteMsgLonLat(file, nx, ny):
    llFile = file

    llShape = (ny,nx)
    llMDI = np.float32(13.5)
    ll = np.fromfile(llFile,dtype=llMDI.dtype)
    lon = ll[0:ny*nx]
    lat = ll[ny*nx:]
    lat.shape = llShape
    lon.shape = llShape

    llsavefile = file.replace('.gra', '')
    np.savez(llsavefile,lon=lon,lat=lat)

# ========================================================================================
# Reads the METEOSAT raw files
# bfile: binary meteosat file
# nx: pixel in x direction
# ny: number in y direction
#
# ========================================================================================
def readMSGraw(bfile, nx, ny, llfile):
    if not os.path.isfile(bfile):
        return np.array([False])

    rrShape = (ny, nx)  # msg shape
    rrMDI = np.uint8(255)
    rr = np.fromfile(bfile, dtype=rrMDI.dtype)
    rr.shape = rrShape
    rr = rr.astype(np.int32) - 173
    msg_latlon = np.load(llfile)

    mlon = msg_latlon['lon']
    mlat = msg_latlon['lat']

    msg_obj = {'t': rr, 'lons': mlon, 'lats': mlat}  # lats lons numpy arrays!

    return msg_obj


def draw_map(t, lat, lon):
    f=plt.figure()
    ax = f.add_subplot(111, projection=ccrs.PlateCarree())
    plt.contourf(lon, lat, t, transform=ccrs.PlateCarree())
    ax.coastlines()
    # Gridlines
    xl = ax.gridlines(draw_labels=True);
    xl.xlabels_top = False
    xl.ylabels_right = False
    # Countries
    ax.add_feature(cartopy.feature.BORDERS, linestyle='--');
    plt.colorbar()
    plt.show()
