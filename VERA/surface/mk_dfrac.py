import numpy as np
import numpy.ma as ma
#import matplotlib
#matplotlib.use("Agg")
from scipy import misc
from scipy.ndimage import gaussian_filter, uniform_filter
import iris
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import iris.plot as iplt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from netCDF4 import Dataset, num2date, MFDataset
import cf_units
import sys

iris.FUTURE.netcdf_no_unlimited = True
file_hurtt = "/prj/vera/ssf/hurtt_dataset/LUH2_v1.0/states.nc"
#file_hurtt_cmip5 = "/scratch/ssf/LUHv1t1_fractions.nc"

names = ["primf", "primn", "secdf", "secdn", "urban", "c3ann", "c4ann", "c3per", "c4per", "c3nfx","pastr", "range", "secmb", "secma"]
names_sub= ["pastr", "c3ann", "c4ann", "c3per", "c4per"]
names_for = ["primf", "primn"]
#names_for = ["primn"]
ncid = Dataset(file_hurtt, "r")

lat = ncid.variables["lat"][:]
lon = ncid.variables["lon"][:]
#print lat

nctime = ncid.variables["time"]

mlon, mlat = np.meshgrid(lon, lat)

croptot = np.zeros_like((mlon))
fortot = np.zeros_like((mlon))
croptot50 = np.zeros_like((mlon))
fortot50 = np.zeros_like((mlon))
land = np.empty_like((mlon))
t1 = 1100
t2 = 1150
print mlon.shape

t1 = 1100 #1950
t2 = 1150 #2000
t2 = 1160 #2010
for i, name in enumerate(names_for):
    data = ncid.variables[name][t2,:,:]
    croptot = croptot + np.ma.masked_invalid(data)
    data50 = ncid.variables[name][t1,:,:]
    croptot50 = croptot50 + np.ma.masked_invalid(data50)
    
print nctime[t1], nctime[t2]
#for i, name in enumerate(names_for):
#    data = ncid.variables[name][1160,:,:]
#    fortot = fortot + data

### create a land sea mask by compbining all datasets
for i, name in enumerate(names):
    data = ncid.variables[name][t2,:,:]
    land = land + np.ma.masked_invalid(data)
    

plt.pcolormesh(mlon, mlat, land.mask)
plt.show()


#sys.exit()



    
## plot comparing the land cover fractions for managed grassland and forest cover in 2010 for CCI-agg (model resolution) and Hurtt (0.5 degress)
dcrop = np.ma.masked_array(croptot - croptot50, mask = land.mask)
dcrop[land.mask] = 0.
land.data[:,:] = 1
land[land.mask] = 0

#dcropfrac = (croptot50 - croptot)/croptot50
dcropfrac = np.ma.masked_array((croptot - croptot50)/croptot50, mask = land.mask)
dcropfrac[land.mask] = 0.
fig = plt.figure(figsize=(10,5))
fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
plt.rcParams.update({"font.size":8})
levels = np.arange(0,1.1,0.1)
cmap = plt.get_cmap('Greens')
norm = mcolors.BoundaryNorm(levels, ncolors=cmap.N, clip=True)

fig.suptitle("Change in primary land")
ax1 = fig.add_subplot(131,title =  "Primary land fraction {}".format(nctime[t1]), projection = ccrs.PlateCarree())
ax1.set_extent([-20, 20, 0, 25])
ax1.pcolormesh(mlon.T, mlat.T, croptot50.T ,vmin= 0, vmax = 1, cmap = cmap, norm = norm)
ax1.coastlines()
ax1.add_feature(cfeature.BORDERS)
ax1.gridlines()

ax2 = fig.add_subplot(132,title =  "Primary land fraction {}".format(nctime[t2]), projection = ccrs.PlateCarree())
ax2.set_extent([-20, 20, 0, 25])
pmap = ax2.pcolormesh(mlon.T, mlat.T, croptot.T ,vmin= 0, vmax = 1, cmap = cmap, norm = norm)
ax2.coastlines()
ax2.add_feature(cfeature.BORDERS)
ax2.gridlines(draw_labels = True)


cax = fig.add_axes([0.2,0.2,0.32,0.02])
fig.colorbar(pmap, cax = cax, orientation = "horizontal")

ax3 = fig.add_subplot(133,title =  "LULCC {}-{}".format(nctime[t1],nctime[t2]), projection = ccrs.PlateCarree())
#levels = np.arange(-0.6,0.7,0.1)
levels = [-0.8,-0.6,-0.4,-0.2,-0.01,0.01,0.2,0.4,0.6,0.8]
cmap = plt.get_cmap('BrBG')
norm = mcolors.BoundaryNorm(levels, ncolors=cmap.N, clip=True)
ax3.set_extent([-20, 20, 0, 25])
pmap = ax3.pcolormesh(mlon.T, mlat.T, dcrop.T ,cmap = cmap, norm = norm)
ax3.coastlines()
ax1.gridlines()
ax3.add_feature(cfeature.BORDERS)

cax = fig.add_axes([0.66,0.2,0.26,0.02])
fig.colorbar(pmap, cax = cax, orientation = "horizontal",spacing = "proportional" )
#plt.tight_layout()
plt.savefig("Hurtt_LUH2_change_in_primary_land_2010_1950.png")
plt.show()

fout = "dfrac_primary_land_2010_1950.nc"

fig = plt.figure(figsize=(10,5))
fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
plt.rcParams.update({"font.size":8})
levels = [-2.,-1.,-0.9,-0.8,-0.6,-0.4,-0.2,-0.01,0.01,0.2,0.4,0.6,0.8,0.9,1.,2.]
cmap = plt.get_cmap('RdBu')
norm = mcolors.BoundaryNorm(levels, ncolors=cmap.N, clip=True)

fig.suptitle("Fractional change in primary land since 1950")
ax1 = fig.add_subplot(111,title =  "Forest fraction {}".format(nctime[t1]), projection = ccrs.PlateCarree())
ax1.set_extent([-20, 20, 0, 25])
pmap = ax1.pcolormesh(mlon.T, mlat.T, dcropfrac.T , cmap = cmap, norm = norm)
ax1.coastlines()
ax1.add_feature(cfeature.BORDERS)
plt.colorbar(pmap)
plt.savefig("dfrac_land.png")
plt.show()

### craete netcdf file of changes
print dcrop.shape
latitude = iris.coords.DimCoord(np.flipud(lat), standard_name = "latitude", units = "degrees")
longitude = iris.coords.DimCoord(lon, standard_name = "longitude", units = "degrees")

c2 = iris.cube.Cube(np.flipud(dcrop), long_name = "absolute_dfrac", 
                        dim_coords_and_dims = [(latitude,0),(longitude,1)])
                        
c3 = iris.cube.Cube(np.flipud(dcropfrac), long_name = "fractional_dfrac", 
                        dim_coords_and_dims = [(latitude,0),(longitude,1)])

c4 = iris.cube.Cube(np.flipud(land.data), long_name = "lsm", 
                        dim_coords_and_dims = [(latitude,0),(longitude,1)])

cube_list = iris.cube.CubeList([c2, c3,c4])
iris.save(cube_list, fout)


## test smoothing the dfrac

dcropmk = np.ma.masked_array(dcrop, mask = data ==  1.e20)

local_mean = uniform_filter(dcropmk, size=3)
sm1 = gaussian_filter(dcrop, sigma=0.5)
sm1mk = np.ma.masked_array(dcrop, mask = data ==  1.e20)
fig = plt.figure(figsize=(10,5))


levels = [-0.8,-0.6,-0.4,-0.2,-0.1,-0.01,0.01,0.1,0.2,0.4,0.6,0.8]
cmap = plt.get_cmap('BrBG')
norm = mcolors.BoundaryNorm(levels, ncolors=cmap.N, clip=True)


fig.suptitle("Change in primary land since 1950-smoothed")

ax1 = fig.add_subplot(121,title =  "Primary land change ", projection = ccrs.PlateCarree())
ax1.set_extent([-20, 20, 0, 25])
ax1.pcolormesh(mlon.T, mlat.T, dcrop.T ,cmap = cmap, norm = norm)
ax1.coastlines()
ax1.add_feature(cfeature.BORDERS)
ax1.gridlines()


ax1 = fig.add_subplot(122,title =  "Primary land change", projection = ccrs.PlateCarree())
ax1.set_extent([-20, 20, 0, 25])
pmap = ax1.pcolormesh(mlon.T, mlat.T, sm1.T , cmap = cmap, norm = norm)
ax1.coastlines()
ax1.add_feature(cfeature.BORDERS)
plt.colorbar(pmap)
plt.savefig("dfrac_smooth.png")
plt.show()

sm1frac = gaussian_filter(dcropfrac, sigma=1)
sm1mkfrac = np.ma.masked_array(sm1frac, mask = data ==  1.e20)



levels = [-2.,-1.,-0.9,-0.8,-0.6,-0.4,-0.2,-0.01,0.01,0.2,0.4,0.6,0.8,0.9,1.,2.]
cmap = plt.get_cmap('RdBu')
norm = mcolors.BoundaryNorm(levels, ncolors=cmap.N, clip=True)
fig.suptitle("Change in primary land since 1950-smoothed")

ax1 = fig.add_subplot(121,title =  "Primary land change ", projection = ccrs.PlateCarree())
ax1.set_extent([-20, 20, 0, 25])
ax1.pcolormesh(mlon.T, mlat.T, dcropfrac.T ,cmap = cmap, norm = norm)
ax1.coastlines()
ax1.add_feature(cfeature.BORDERS)
ax1.gridlines()


ax1 = fig.add_subplot(122,title =  "Primary land change", projection = ccrs.PlateCarree())
ax1.set_extent([-20, 20, 0, 25])
pmap = ax1.pcolormesh(mlon.T, mlat.T, sm1mkfrac.T , cmap = cmap, norm = norm)
ax1.coastlines()
ax1.add_feature(cfeature.BORDERS)
plt.colorbar(pmap)
plt.savefig("dfrac_fractional_smooth.png")
plt.show()



## save smmothed data set to file 
fout = "dfrac_primary_land_2010_1950_guassian.nc"
latitude = iris.coords.DimCoord(np.flipud(lat), standard_name = "latitude", units = "degrees")
longitude = iris.coords.DimCoord(lon, standard_name = "longitude", units = "degrees")

c4 = iris.cube.Cube(np.flipud(sm1), long_name = "absolute_dfrac", 
                        dim_coords_and_dims = [(latitude,0),(longitude,1)])
                        
c5 = iris.cube.Cube(np.flipud(sm1mkfrac), long_name = "fractional_dfrac", 
                        dim_coords_and_dims = [(latitude,0),(longitude,1)])

cube_list = iris.cube.CubeList([c4, c5])
iris.save(cube_list, fout)


## plot histgram over region


