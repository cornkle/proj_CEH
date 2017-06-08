import salem
from salem.utils import get_demo_file
import xarray as xr
import matplotlib.pyplot as plt
import pdb
import numpy as np
from functools import partial
from salem import get_demo_file, open_xr_dataset, GeoTiff, wgs84



dat = xr.open_dataarray('/users/global/cornkle/data/pythonWorkspace/proj_CEH/topo/gtopo_1min.nc')
#dat=dat.sel(lon=slice(-18,120), lat=slice(-30,60))

#dat=dat.sel(lon=slice(-100,-40), lat=slice(-30,60))

grid = dat.salem.grid

grid50 = grid.regrid(factor=0.03)
lakes = salem.read_shapefile(salem.get_demo_file('ne_50m_rivers_lake_centerlines.shp'), cached=True)


top_on_grid50 = grid50.lookup_transform(dat, method=np.std)

sm = dat.salem.get_map(cmap='topo')
lakes = salem.read_shapefile(salem.get_demo_file('ne_50m_rivers_lake_centerlines.shp'), cached=True)
sm.set_shapefile(lakes, edgecolor='k', facecolor='none', linewidth=2,)
mask_lakes = grid.region_of_interest(shape=lakes)
sm.set_data(top_on_grid50, grid50)
sm.set_plot_params(vmin=20, vmax=500)
#sm.set_data(dat, grid)
sm.visualize()

f = plt.figure()

sm.set_data(dat, grid)
sm.set_plot_params(vmax=1000)
sm.visualize()
