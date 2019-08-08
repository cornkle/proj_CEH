import numpy
import xarray as xr
from wavelet import util


data = xr.open_dataarray('/home/ck/DIR/cornkle/data/VERA/lsta_trend_chris/mean_lst_trend.nc')

data.values[data.values<0] = 0

