import rioxarray
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
from pyproj import CRS


# MTG vis file
file = '/scratch/cmt/FCIL1HRFI_20250326T074107Z_20250326T074726Z_epct_ae3ffc0b_FC.tif'
# Directory path and name for new figure
outfile = '/scratch/cmt/vis_view.jpg'

# read reprojection information - obtained from dummy netcdf and written to txt for general use.
with open("/prj/nflics/MTG_vis_testbed/crs_proj_geostationary_MTG.txt", "r") as f:
    crs = CRS.from_proj4(f.read())
    
# read tif file function   
def read_geotiff(file):
    data = rioxarray.open_rasterio(file)
    data_rio = data.rio.write_crs(crs)
    data_reprojected = data_rio.rio.reproject("EPSG:4326")

    # to cut down domain size or mask out missing values if needed:
    #data_reprojected = data_reprojected.sel(x=slice(13,36), y=slice(-9, -35))
    #data_reprojected = data_reprojected.where(data_reprojected.values<60000, other=0)
    
    return data_reprojected.squeeze()

# read in file, da_box data array can be used for analyses.
da_box = read_geotiff(file)

# plot it!
f = plt.figure(figsize=(13,9), dpi=200)
ax = f.add_subplot(111, projection=ccrs.PlateCarree())
plt.contourf(da_box.x.values, da_box.y.values, da_box, transform=ccrs.PlateCarree(), cmap='viridis', extend='both')
ax.coastlines()

xl = ax.gridlines(draw_labels=True);

ax.add_feature(cartopy.feature.BORDERS, linestyle='--', color='grey');

cbar = plt.colorbar()

f.savefig(outfile)

print('Written figure ', outfile)