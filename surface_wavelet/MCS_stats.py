import numpy as np
import xarray as xr
from utils import constants
import matplotlib.pyplot as plt
import pdb


MCSdom = xr.open_dataarray(constants.MCS_POINTS_DOM)
MCSdom = MCSdom.sel(lat=slice(10,18), lon=slice(-10,10))

# f = plt.figure()
# plt.hist(MCSdom.values[np.isfinite(MCSdom.values)])
plt.title('MCS dominant scale distribution')

count = MCSdom.groupby(MCSdom['time.hour']).count()

f = plt.figure()
plt.plot(count.hour, count)

u35 = MCSdom.where((MCSdom.values<35) & (MCSdom.values>7)).groupby(MCSdom['time.hour']).count()
o60 = MCSdom.where((MCSdom.values>60)).groupby(MCSdom['time.hour']).count()

u35_frac = u35/count
o60_frac = o60/count

f = plt.figure()
plt.plot(count.hour, u35, 'ro')
plt.title('Under 35km')

f = plt.figure()
plt.plot(count.hour, o60, 'ko')
plt.title('Over 60km')

f = plt.figure()
plt.plot(count.hour, u35_frac, 'ro')
plt.plot(count.hour, o60_frac, 'ko')

