import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from utils import u_arrays as ua
from collections import OrderedDict
import pandas as pd
import multiprocessing
import pickle as pkl
from scipy.ndimage.measurements import label
import pdb


dic = pkl.load( open ('/users/global/cornkle/C_paper/wavelet/saves/bulk_40big_zR.p', 'rb'))

p30 = np.array(dic['po30'])
lat = np.array(dic['clat'])
mcs_count = np.sum(p30[(lat >= 4) & (lat<=7.5)])



files = ua.locate(".nc", '/users/global/cornkle/TRMMfiles')
cnt = 0
for f in files:
    print('Doing ', f)
    xa = xr.open_dataset(f)

    lat = xa.lat.values
    lon = xa.lon.values
    arr = xa['p'].values
    arr = arr[(lat>=4) & (lat <= 7.8) & (lon >=-17) & (lon <=20)]
    nb = np.sum(arr >= 30)

    cnt += nb


print('MCS frac', mcs_count/cnt)