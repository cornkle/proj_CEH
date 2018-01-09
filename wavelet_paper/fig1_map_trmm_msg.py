import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import salem
import pdb
from utils import u_grid
from scipy.interpolate import griddata
import xarray as xr
from numpy import ma
from scipy.ndimage.measurements import label


def create_map_data():
    # read only trmm files that I need and give out proper lons lats etc
    files = "/users/global/cornkle/data/OBS/TRMM/trmm_swaths_WA/2011/06/2A25.20110612.77322.7.gra"

     # /2011/06/2A25.20110612.77322.7.gra"  good to show
    trr = np.fromfile(files,dtype=np.int16)
    x = 49
    nb = trr.size
    single = int(nb/4) # variables lon lat rainrate flag

    lons = trr[0:single]
    lats = trr[single:2*single]
    rainrs = trr[2*single:3*single]
    flags = trr[3*single:4*single]

    y = int(lons.size/x)
    lons = np.resize(lons, (y,x))
    lats = np.resize(lats, (y,x))
    rainrs = np.resize(rainrs, (y,x))
    flags = np.resize(flags, (y,x))
    lon=lons/100.
    lat=lats/100.
    rainr=rainrs/10.
    lonmin, lonmax=np.amin(lon),np.amax(lon)
    latmin, latmax=np.amin(lat),np.amax(lat)
    lonx=lon[0,:]
    laty=lat[:,0]
    rainrs.shape

    path = "/users/global/cornkle/data/OBS/meteosat_WA30/cell_blob_files/2011/06/"  # 201106122130 good to show
    filename = "201106122130.gra"
    files = path + filename
    rrShape = (580,1640)
    rrMDI = np.uint16()
    rr = np.fromfile(files,dtype=rrMDI.dtype)
    rr.shape = rrShape

    path = "/users/global/cornkle/data/OBS/meteosat_WA30/msg_raw_binary/2011/06/"  # 201106122130 good to show
    filename = "201106122130.gra"
    files = path + filename
    rrMDI = np.uint8(255)
    rr2 = np.fromfile(files, dtype=rrMDI.dtype)
    rr2.shape = rrShape
    rr2 = rr2.astype(np.int32) - 173

    msg_latlon=np.load('/users/global/cornkle/data/OBS/meteosat_WA30/MSG_1640_580_lat_lon.npz')
    mlon = msg_latlon['lon']
    mlat = msg_latlon['lat']

    # make salem grid
    grid = u_grid.make(mlon, mlat, 5000)
    xi, yi = grid.ij_coordinates
    glon, glat = grid.ll_coordinates

    # Transform lons, lats to grid
    xt, yt = grid.transform(lon.flatten(), lat.flatten(), crs=salem.wgs84)

    # Convert for griddata input
    tpoints = np.array((yt, xt)).T
    inter = np.array((np.ravel(yi), np.ravel(xi))).T

    # Interpolate using delaunay triangularization
    dummyt = griddata(tpoints, rainrs.flatten(), inter, method='linear')
    outt = dummyt.reshape((grid.ny, grid.nx))

    for nb in range(5):
        boole = np.isnan(outt)
        outt[boole] = -1000
        grad = np.gradient(outt)
        outt[boole] = np.nan
        outt[abs(grad[1]) > 300] = np.nan
        outt[abs(grad[0]) > 300] = np.nan

    xm, ym = grid.transform(mlon.flatten(), mlat.flatten(), crs=salem.wgs84)
    mpoints = np.array((ym, xm)).T
    out = griddata(mpoints, rr2.flatten(), inter, method='linear')
    outp = out.reshape((grid.ny, grid.nx))

    out = griddata(mpoints, rr.flatten(), inter, method='nearest')
    outb = out.reshape((grid.ny, grid.nx))

    data = xr.Dataset({'trmm': (['lat', 'lon'], outt),
                       'tir': (['lat', 'lon'], outp),
                       'tblob' : (['lat', 'lon'], outb)},
             coords={ 'lat': glat[:,0], 'lon':glon[0,:]}) #[np.newaxis, :]

    data.to_netcdf('/users/global/cornkle/C_paper/wavelet/saves/maps/trmm_msg_map.nc')

def plot_data():

    fpath = '/users/global/cornkle/C_paper/wavelet/figs/paper/'

    data = xr.open_dataset('/users/global/cornkle/C_paper/wavelet/saves/maps/trmm_msg_map.nc')
    data=data.sel(lon=slice(-17,20), lat=slice(4.,20))
    map = data.salem.get_map(cmap='inferno')

    outb = data['tir']


    map.set_data(data['tir'])
    map.set_shapefile(countries=True, linewidth=0.1 )

    f= plt.figure(figsize=(9, 6.8), dpi=300)
    ax1 = f.add_subplot(211)
    ax2 = f.add_subplot(212)
    map.set_lonlat_contours(add_xtick_labels=False)
    map.set_plot_params(cmap='inferno',  vmin=-80, vmax=20, extend=None)
    map.visualize(ax=ax1, cbar_title='Brightness temperature ($^{\circ}$C)',)

    outb[outb > -40] = 0
    labels, numL = label(outb.values)


    for i in np.unique(labels):
        pos = np.where(labels == i)
        if len(pos[0]) < 65:
            labels[pos]=0

    labels, numL = label(labels)
    labels = np.array(labels, dtype=float)
    labels[labels == 0] = np.nan
    pos = np.where(labels>0)

    labels[pos] = outb.values[pos]
    #
    # f = plt.figure()
    # plt.imshow(labels)

    outb = np.array(data['tblob'], dtype=float)
    outb[outb == 0] = np.nan

    outboth = data['trmm'].copy().values
    outboth[outboth<=0] = 10

   # outboth[np.isfinite(labels)] = labels[np.isfinite(labels)]

    map.set_plot_params(cmap='Greys',  vmin =-150, vmax = 20,  extend=None)
    map.set_data(outboth)
    # map.set_contour(data['trmm'].values, cmap='winter')
    map.visualize(ax=ax2, addcbar=False)

    map.set_plot_params(cmap='Greys', vmin =-150, vmax = 20, extend=None)
    map.set_lonlat_contours(add_xtick_labels=True)
    map.set_data(labels)
    map.visualize(ax=ax2, addcbar=False)

    map.set_contour(data['trmm'].values, cmap='winter')
    map.visualize(ax=ax2, addcbar=False)

    plt.annotate('a)', xy=(0.01, 0.96), xytext=(0, 4), size=15, xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points')
    plt.annotate('b)', xy=(0.01, 0.5), xytext=(0, 4), size=15, xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points')

    plt.tight_layout()

    pos2 = ax2.get_position()
    pos2 = [pos2.x0-0.025, pos2.y0, pos2.width, pos2.height]
    ax2.set_position(pos2)

    plt.savefig(fpath + 'map_example.png')
    plt.close('all')

if __name__ == "__main__":
    plot_data()