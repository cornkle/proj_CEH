# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 17:27:16 2016

@author: cornkle
"""

import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy
import numpy as np
import pdb


def quick_map(xar, save = None):

    f = plt.figure(figsize=(15, 7), dpi=400)
    ax = plt.axes(projection=ccrs.PlateCarree())
    xar.plot.contourf('lon', 'lat', projection=ccrs.PlateCarree(), vmax=50)
    ax.coastlines()
    # Gridlines
    xl = ax.gridlines(draw_labels=True);
    xl.xlabels_top = False
    xl.ylabels_right = False
    # Countries
    ax.add_feature(cartopy.feature.BORDERS, linestyle='--');
    if save:
        plt.savefig(save)
    else:
        plt.show()

def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:

    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)

    #return base.from_list(cmap_name, color_list, N)
    return plt.cm.get_cmap(base_cmap, N)