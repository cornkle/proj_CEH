import xarray as xr
import pdb
import numpy as np
import glob
import os
import salem



file = '/users/global/cornkle/figs/VERA/bamba/w2018_bamba/intest/xmhkja.pc20140405_12.nc'

data = salem.open_metum_dataset(file)

data = data['STASH_m01s03i234_2']

