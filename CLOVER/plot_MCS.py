import xarray as xr
import glob
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
from utils import u_met
import numpy as np


def draw_map(t, p, lat, lon):
    f=plt.figure()
    ax = f.add_subplot(111, projection=ccrs.PlateCarree())
    plt.contourf(lon, lat, t, transform=ccrs.PlateCarree(), levels=np.arange(-75,-39,5))
    plt.colorbar()
    plt.contour(lon, lat, p, transform=ccrs.PlateCarree(), levels=np.arange(10,51,10), cmap='jet')
    #ax.coastlines()
    # Gridlines
    xl = ax.gridlines(draw_labels=True);
    xl.xlabels_top = False
    xl.ylabels_right = False
    # Countries
    #ax.add_feature(cartopy.feature.BORDERS, linestyle='--');

    plt.show()


def get_file(f, var):

    file = xr.open_dataset(f)

    data= file[var]

    if var == 'lw_out_PBLtop':

        data.values = u_met.OLR_to_Tb(data)-273.15
        p = file['lsRain']*3600

        draw_map(data.values, p.values, data.latitude, data.longitude)
    else:
        p = file['p']
        draw_map(data.values, p.values, data.lat, data.lon)



def plot():

    var = 'lw_out_PBLtop'
    #var = 'tc_lag0'
    CP4files = glob.glob('/users/global/cornkle/data/CP4/CLOVER/MCS/*.nc')
    #MSGfiles = glob.glob('/users/global/cornkle/MCSfiles/WA5000_4-8N_13W-13E_-40_18UTC/*.nc')

    for a in range(32,50):    #  # 300,325
        get_file(CP4files[a],var)