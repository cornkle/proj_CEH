import seaborn as sns
import pickle as pkl
import glob
import xarray as xr
import numpy as np
import pdb
files = glob.glob('/users/global/cornkle/data/CP4/CLOVER/MCS_-60_1000km2/*.nc' )




for f in files:

     data = xr.open_dataset(f)

     if int(data['time.month'].values) != 10:
         continue
     tis = data['lw_out_PBLtop'].values
     if (np.sum(np.isfinite(tis))*(4.4**2))>=500000:
         print (f)
