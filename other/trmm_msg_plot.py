import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

def plot_data():

    fpath = '/users/global/cornkle/C_paper/wavelet/figs/paper/'

    data = xr.open_dataset('/users/global/cornkle/C_paper/wavelet/saves/maps/trmm_msg_map.nc')
    map = data.salem.get_map(cmap='inferno')

    outb = np.array(data['tblob'], dtype=float)
    outb[outb==0] = np.nan

    outboth = data['trmm'].copy().values
    outboth[outb>0] = outb[outb>0]

    map.set_data(data['tir'])
    map.set_shapefile(countries=True, linewidth=0.1 )

    f= plt.figure()
    ax1 = f.add_subplot(211)
    ax2 = f.add_subplot(212)
    map.visualize(ax=ax1, cbar_title='Cloud top temperature ($^{\circ}$C)')


    map.set_plot_params(cmap='hot', vmax=91)
    map.set_data(outboth)
    map.set_contour(data['trmm'], cmap='winter')
    map.visualize(ax=ax2, cbar_title='Cloud ID | < -40$^{\circ}$C')

    plt.tight_layout()
    plt.savefig(fpath + 'map_example.png')