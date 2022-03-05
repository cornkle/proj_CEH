# -*- coding: utf-8 -*-


import multiprocessing
import glob
from eod import rewrite_data
import ipdb
import os
import xarray as xr
import numpy as np
import pandas as pd
from utils import u_darrays as uda
from scipy.interpolate import griddata


VARS = {  'rlut' : ['toa_outgoing_longwave_flux', 60, '2d'],  # fname, varname, interval (minutes)
          'pr' : ['precipitation_flux', 15, '2d'],
          'prw' : ['atmosphere_water_vapor_content', 15, '2d'],
          'clivi' : ['atmosphere_mass_content_of_cloud_ice', 15, '2d'],
          'hfls' : ['surface_upward_latent_heat_flux', 15, '2d'],
          'hfss' : ['surface_upward_sensible_heat_flux', 15, '2d'],
          'tas' : ['air_temperature', 15, '2d'], # 2m temp
          'ua' : ['eastward_wind', 390, '3d'],
          'va' : ['northward_wind', 390, '3d'],
          'uas': ['eastward_wind', 390, '2d'], # at 10m
          'vas': ['northward_wind', 15, '2d'], # at 10m
          'wap' : ['lagrangian_tendency_of_air_pressure', 15, '3d'], # omega on plev/uv grid
          'ps' : ['surface_air_pressure', 15, '2d'],
          'hur' : ['relative_humidity', 15, '3d'],
          'huss': ['specific_humidity', 15, '2d'], # at 1.5m
          'mrsol' : ['moisture_content_of_soil_layer', 15, '3d'], # 4 levels
          'rlds'  : ['surface_downwelling_longwave_flux_in_air', 60, '2d'],
          'rlus'  : ['surface_upwelling_longwave_flux_in_air', 60, '2d']

          }


