import xarray as xr
import glob
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
from utils import u_met
import numpy as np
from utils import constants as cnst


def draw_map(t, p, lat, lon):
    f=plt.figure()
    ax = f.add_subplot(111, projection=ccrs.PlateCarree())
    plt.contourf(lon, lat, t, transform=ccrs.PlateCarree(), levels=np.arange(-65,-50))  #levels=np.arange(-75,-39,5)
    plt.colorbar()
    plt.contour(lon, lat, p, transform=ccrs.PlateCarree(), levels=np.arange(2,15,2), cmap='jet') #np.arange(10,51,10)
    #ax.coastlines()
    # Gridlines
    xl = ax.gridlines(draw_labels=True);
    xl.xlabels_top = False
    xl.ylabels_right = False
    # Countries
    #ax.add_feature(cartopy.feature.BORDERS, linestyle='--');
    plt.show()


def get_file(f, var):
    #f = CP4files[6]
    file = xr.open_dataset(f)
    #file['lw_out_PBLtop'].plot.contourf()
    data= file[var]
    if var == 'lw_out_PBLtop':

        #data.values = u_met.OLR_to_Tb(data.values)-273.15

        if np.sum(np.isfinite(data.values)) <=0:
            return

        p = file['totRain']*3600

        draw_map(data.values, p.values, data.latitude, data.longitude)
    else:
        p = file['p']
        draw_map(data.values, p.values, data.lat, data.lon)



def plot():

    var = 'lw_out_PBLtop'
    #var = 'tc_lag0'
    CP4files = glob.glob('/users/global/cornkle/data/CP4/CLOVER/MCS/*.nc')
    CP4files = glob.glob(cnst.network_data + 'data/CP4/CLOVER/MCS25_-50_5000km2/*.nc')
    #MSGfiles = glob.glob('/users/global/cornkle/MCSfiles/WA5000_4-8N_13W-13E_-40_18UTC/*.nc')

    CP4files

    for a in range(230,240):    #  # 300,325
        get_file(CP4files[a],var)
