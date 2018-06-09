#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 31 11:04:44 2018

@author: adam
"""

import iris
import datetime
import cf_units as unit
import iris.coord_categorisation
from iris.coord_categorisation import add_categorised_coord
import iris.coords as coord
import iris.quickplot as qplt
import iris.plot as iplt
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt
import os # operating system interface
import pdb


def iris_original():
    cubes= iris.load('/users/global/cornkle/w2018_bamba/mini_forest/pb20140407_forest.nc')
    print (cubes)

    atemp = cubes[2]

    print (atemp)

    veg_lat=iris.Constraint(grid_latitude_uv=lambda cell: 5 <= cell <= 6.0) # latitude constraint 4-25N
    veg_lon=iris.Constraint(grid_longitude_uv=lambda cell: -7.5 <= cell <= -6.9) # longitude constraint 18W - 25E
    temp_box = atemp.extract(veg_lat & veg_lon) # apply constraints
    ts = temp_box.collapsed(['grid_latitude_uv', 'grid_longitude_uv'], iris.analysis.MEAN)

    #print (ts.coord('time'))

    U = unit.Unit('hours since 2014-04-07 01:00:00', calendar='gregorian')  # define a new Iris unit instance whose timestep is days
    ts.coord('time').convert_units(U)

    print (ts)
    #print (ts.units.date2num(datetime.datetime(2014, 4, 7)))

    #pdb.set_trace()

    nveg_lat=iris.Constraint(grid_latitude_uv=lambda cell: 5.2 <= cell <= 6.2) # latitude constraint 4-25N
    nveg_lon=iris.Constraint(grid_longitude_uv=lambda cell: -6.5 <= cell <= -6.0) # longitude constraint 18W - 25E
    ntemp_box = atemp.extract(nveg_lat & nveg_lon) # apply constraints

    nts = ntemp_box.collapsed(['grid_latitude_uv', 'grid_longitude_uv'], iris.analysis.MEAN)

    sU = unit.Unit('hours since 2014-04-07 01:00:00', calendar='gregorian')  # define a new Iris unit instance whose timestep is days
    nts.coord('time').convert_units(sU)

    diff = ts - nts

    pdb.set_trace()
    #temp = ts[0,:,:]
    #print (ts)

    #pdb.set_trace()
    fig = plt.figure(figsize = (12,8))
    ax = fig.add_subplot(3,1,1)
    ax.contourf(ts, nts.coord('time').points, nts.coord(), cmap=cm.jet) #coord=['grid_longitude', 'grid_latitude'])
    plt.colorbar()
    plt.title('relative humidity over vegetated area 2014/04/07')

    plt.subplot(3,1,2)
    iplt.contourf(nts, 15, cmap=cm.jet) #coord=['grid_longitude', 'grid_latitude'])
    plt.colorbar()
    plt.title('relative humidity over non-vegetated area 2014/04/07')

    plt.subplot(3,1,3)
    iplt.contourf(diff, 15, cmap=cm.seismic) #coord=['grid_longitude', 'grid_latitude'])
    plt.colorbar()
    plt.title('rh diff [Veg-NonVeg] 2014/04/07')

    plt.subplots_adjust(hspace=0.3)
    plt.savefig('rh_veg_nveg07_diff.png')
    plt.show()




    #print pr
    #print pr.shape
    #print data.dimensions.keys()
    #plt.imshow(pr[0], origin='lower') #data for the first time-step (January 1948)
    #plt.show()



#for v in data.variables:
#    print v

#data.variables['pseudo_level'][:]
#print data.shape
