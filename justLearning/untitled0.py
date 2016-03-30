# -*- coding: utf-8 -*-
"""
Created on Sat Mar  5 22:02:01 2016

@author: cklein
"""
# Define the tools we are going to need today
import matplotlib.pyplot as plt  # plotting library
import numpy as np  # numerical library
import xray  # NetCDF library
import cartopy  # Plotting libary
import cartopy.crs as ccrs  # Projections


def reformERA():

    netcdf = xray.open_dataset('/media/cklein/Elements/sens_stud/files/ERA/ERA_prcp_monthly_79-2013.nc')  
    t2_avg = netcdf.tp.mean(dim='time')*24.*31.
    t2_max = netcdf.tp.max(dim='time')*24.*31.
    box=t2_avg.sel()  #latitude=slice(18, 4), longitude=slice(-10,10)
    # Define the map projection.
    ax = plt.axes(projection=ccrs.Robinson())
    #box=t2_max
    # ax is an empty plot. We now plot the variable t2_avg onto ax
    box.plot(ax=ax, transform=ccrs.PlateCarree(), cmap=plt.get_cmap('summer'), levels=10) 
    # The keywords "origin" and "transform" are projection details 
    # The keyword aspect just ensures that the plot aspect ratio is preserved
    ax.coastlines()  # Add coastlines to the plot
    plt.show(ax)

reformERA()
