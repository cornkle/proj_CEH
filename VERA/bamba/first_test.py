#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 29 16:08:33 2018

@author: adam
"""

import iris
import iris.coord_categorisation
import iris.coords as coord
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import numpy as np
from iris.time import PartialDateTime
import datetime
import ipdb


filenames = [('/media/adam/CHEICK/w2018/xmhkga.pc20140405_00'),('/media/adam/CHEICK/w2018/xmhkga.pc20140405_12')]
cubes = iris.load(filenames)
#print (cubes)
precipitation_flux = cubes[11]

print (precipitation_flux.coord('time'))

#reg_cubes = precipitation_flux.intersection(grid_longitude=(-15.0, 25.0), grid_latitude=(4.0, 31.0))
#reg_cubes.coord('latitude').guess_bounds()
#reg_cubes.coord('longitude').guess_bounds()

lat_WA = iris.Constraint(grid_latitude=lambda cell: 0.0 <= cell <= 16.0)
lon_WA = iris.Constraint(grid_longitude=lambda cell: 340.0 <= cell <= 380.0)
reg_cubes = precipitation_flux.extract(lat_WA & lon_WA)
print (reg_cubes)

iris.coord_categorisation.add_hour(reg_cubes, 'time', name='minute')

cubes_time_cumul = reg_cubes.aggregated_by('minute', iris.analysis.SUM)

print (cubes_time_cumul)

#reg_cubes.convert_units('kg s-1')
#reg_cubes.rename('Precipitation_rate')
#reg_cubes.units='mm s-1'




#
#
#print (precipitation_flux)

#cumul = precipitation_flux.collapsed('time', iris.analysis.SUM)

