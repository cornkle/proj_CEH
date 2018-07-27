import salem
import pyproj
import numpy as np
from scipy.interpolate import griddata
from scipy.ndimage.measurements import label
import datetime as dt
from eod import msg, trmm, tm_utils, trmm_clover
import xarray as xr
import os
import ipdb
import matplotlib.pyplot as plt
from utils import constants
from utils import u_grid
import pandas as pd
import pdb


def test_scatter():


path = constants.CLOVER_HIST
rain = constants.CP4_RAIN
olr = constants.CP4_OLR

file_save(rain, olr)