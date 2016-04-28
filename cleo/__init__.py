"""
Cleo
====

Provides some graphical tools.

Copyright: Fabien Maussion, 2014-2015

License: GPLv3+
"""
from __future__ import division

from os import path

# Path to the file directory
file_dir = path.join(path.abspath(path.dirname(__file__)), 'files')
files = dict()
files['world_borders'] = path.join(file_dir, 'shapes', 'world_borders',
                                   'world_borders.shp')
files['oceans'] = path.join(file_dir, 'shapes', 'oceans',
                                    'ne_50m_ocean.shp')
files['rivers'] = path.join(file_dir, 'shapes', 'rivers',
                                    'ne_50m_rivers_lake_centerlines.shp')

# API
from cleo.colors import get_cm
from cleo.graphics import DataLevels
from cleo.graphics import Map