# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 17:27:16 2016

@author: cornkle
"""

import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy


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