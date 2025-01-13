import rioxarray
import xarray as xr
import numpy as np
import glob
from pyproj import CRS
from utils import u_interpolate_small as inter

# MTG VIS file location
files = glob.glob('/prj/nflics/MTG_vis_testbed/FCI_tif_NRT/*.tif')

#read reprojection information
with open("/prj/nflics/MTG_vis_testbed/crs_proj_geostationary_MTG.txt", "r") as f:
    crs = CRS.from_proj4(f.read())

# read tif file function
def read_geotiff(file):
    data = rioxarray.open_rasterio(file)
    data_rio = data.rio.write_crs(crs)
    data_reprojected = data_rio.rio.reproject("EPSG:4326")

    data_reprojected = data_reprojected.sel(x=slice(13,36), y=slice(-9, -35))
    data_reprojected = data_reprojected.where(data_reprojected.values<60000, other=0)

    #normalise data between 0 and 1 for visualisation
    data_reprojected = (data_reprojected - np.min(data_reprojected)) / (np.max(data_reprojected) - np.min(data_reprojected))

    return data_reprojected.squeeze()

# read in specified file (SHOULD BE CHANGED TO DYNAMIC USER INPUT)
da_box = read_geotiff(files[120])
#############

xmax, xmin = np.max(da_box.x), np.min(da_box.x)
ymax, ymin = np.max(da_box.y), np.min(da_box.y)

# create new coordinates going from ~500m to ~2km by dividing by 4
new_x = np.linspace(xmin.values, xmax.values, int(len(da_box.x)/4))
new_y = np.linspace(ymin.values, ymax.values, int(len(da_box.y)/4))


# usual interpolation routine, input is not irregular
# THIS SHOULD BE IMPLEMENTED IN WORKFLOW TO BE RUN ONLY ONCE AND "ind, weights, shape" OUTPUT SAVED FOR REUSE
inds, weights, shape = inter.interpolation_weights(da_box.x.values, da_box.y.values, new_x, new_y, irregular_1d=False)
##################

# read reprojected data as numpy array
data_2km = inter.interpolate_data(da_box.values, inds, weights, shape)

# convert to data array for subsetting
da_2km = xr.DataArray(data_2km,
    coords={"y": new_y, "x": new_x},
    dims=["y", "x"],
)

# subsetting to remove undefined edges for portal plotting
da_2km = da_2km.sel(y=slice(-33,-10), x=slice(15,35))



######### da_2km should be saved in format needed for portal hosting.



