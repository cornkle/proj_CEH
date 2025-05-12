import geopandas, pickle, time
import numpy as np
import numpy.ma as ma
from datetime import datetime
from netCDF4 import Dataset
from sklearn.neighbors import BallTree
from tqdm import tqdm


def get_nearest(src_points, candidates, k_neighbors=1):
    """Find nearest neighbors for all source points from a set of candidate points"""
    # Create tree from the candidate points
    tree = BallTree(candidates, leaf_size=15, metric='haversine')
    # Find closest points and distances
    distances, indices = tree.query(src_points, k=k_neighbors)
    # Transpose to get distances and indices into arrays
    distances = distances.transpose()
    indices = indices.transpose()
    # Get closest indices and distances (i.e. array at index 0)
    # note: for the second closest points, you would take index 1, etc.
    closest = indices[0]
    closest_dist = distances[0]
    # Return indices and distances
    return (closest, closest_dist)


def nearest_neighbour(left_gdf, right_gdf, return_dist=False):
    """
    For each point in left_gdf, find closest point in right GeoDataFrame and return them.
    NOTICE: Assumes that the input Points are in WGS84 projection (lat/lon).
    """
    left_geom_col = left_gdf.geometry.name
    right_geom_col = right_gdf.geometry.name
    # Ensure that index in right gdf is formed of sequential numbers
    right = right_gdf.copy().reset_index(drop=True)
    # Parse coordinates from points and insert them into a numpy array as RADIANS
    left_radians = np.array(left_gdf[left_geom_col].apply(lambda geom: (geom.y * np.pi / 180, geom.x * np.pi / 180)).to_list())
    right_radians = np.array(right[right_geom_col].apply(lambda geom: (geom.y * np.pi / 180, geom.x * np.pi / 180)).to_list())
    # Find the nearest points
    # -----------------------
    # closest ==> index in right_gdf that corresponds to the closest point
    # dist ==> distance between the nearest neighbors (in meters)
    closest, dist = get_nearest(src_points=left_radians, candidates=right_radians)
    # Return points from right GeoDataFrame that are closest to points in left GeoDataFrame
    closest_points = right.loc[closest]
    # Ensure that the index corresponds the one in left_gdf
    closest_points = closest_points.reset_index(drop=True)
    # Add distance if requested
    if return_dist:
        # Convert to meters from radians
        earth_radius = 6371  # km
        closest_points['distance'] = dist * earth_radius
    return closest_points


def neighbour_tree():
    """
    Use a Ball Tree to find the nearest neigbour ocean pixel for each deforested pixel
    and compute the distance between them.

    Assumes that aggregate_forest_loss.py has been run to generate a map of forest loss
    at 0.05 degree resolution.
    Contains code to read MODIS Percent_land_in_grid layer from MOD11C3 to use as land/ocean mask.
    Alternative masks at 0.05 degree resolution may be substituted. 
    Distance data is saved using pickle for later use.
    Parameters:
    None
    Returns:
    None
    """
    # Create coordinates for global map at 0.05 degree resolution
    lats = np.arange(-90, 90, 0.05) + 0.5 * 0.05
    lons = np.arange(-180, 180, 0.05) + 0.5 * 0.05
    # Construct ocean mask from MODIS land percentage (ocean if land covers 0% of 0.05 deg pixel)
    with Dataset("../data/modis_cmg_land_percentage.nc", 'r') as data:
        ocean_map = (data.variables['Land_percent'][:] == 0.).T
    # Read forest loss data at 0.05 degree resolution 
    with Dataset("../data/hansen_0pt05/global_percent_forest_loss_0pt05deg.nc", 'r') as data:
        defn_pc = data.variables['forest_cover_loss'][:]
    # Find centre coordinates of 0.05 degree grid boxes with >0% forest loss
    lats_deforestation = ma.where(deforestation_pc > 0.)[0]
    lons_deforestation = ma.where(deforestation_pc > 0.)[1]
    # Construct geopandas data frame of coordinates
    deforestation_points = geopandas.points_from_xy(lons[lons_deforestation], 
                                                    lats[lats_deforestation], crs='EPSG:4326')
    deforestation_geom = {'geometry': deforestation_points}
    gdf_deforestation = geopandas.GeoDataFrame(deforestation_geom, geometry=deforestation_points)
    print('Created data frame of deforested grid box centres.')
    # Find centre coordinates of 0.05 degree grid boxes with no land
    lats_ocean = np.where(ocean_map)[0]
    lons_ocean = np.where(ocean_map)[1]
    # Construct geopandas data frame of coordinates
    ocean_points = geopandas.points_from_xy(lons[lons_ocean], lats[lats_ocean], crs='EPSG:4326')
    ocean_geom = {'geometry': ocean_points}
    gdf_ocean = geopandas.GeoDataFrame(ocean_geom, geometry=ocean_points)
    print('Created data frame of ocean grid box centres.')
    # For each deforested grid box, compute nearest ocean grid box and distance
    closest_ocean = nearest_neighbour(gdf_deforestation, gdf_ocean, return_dist=True)
    closest_ocean = closest_ocean.rename(columns={'geometry': 'closest_ocean_geom'})
    # Match up neighbour/distance data with the corresponding deforested grid box coordinates
    deforestation_distance = gdf_deforestation.join(closest_ocean)
    # Save data frame of deforested grid box centres, nearest neighbour ocean grid box centres, distances
    pickle.dump(deforestation_distance, open('../data/deforestation_distance.p', 'wb'))
    print('Saved nearest neighbour locations and distances.')
    

def write_distance_map():
    deforestation_distance = pickle.load(open('../data/deforestation_distance.p','rb'))
    global_lat = np.arange(-90, 90, 0.05) + 0.5 * 0.05
    global_lon = np.arange(-180, 180, 0.05) + 0.5 * 0.05
    distance_map = np.zeros((global_lat.size, global_lon.size)) * np.nan
    n = len(deforestation_distance)
    for px in tqdm(range(n), desc='rasterising points'):
        deforestation_point = deforestation_distance['geometry'][px]
        px_lon = deforestation_point.x
        px_lat = deforestation_point.y
        distance = deforestation_distance['distance'][px]
        lat_idx = np.argmin(np.abs(global_lat - px_lat))
        lon_idx = np.argmin(np.abs(global_lon - px_lon))
        distance_map[lat_idx, lon_idx] = distance
    f = Dataset('../data/distance_from_ocean_0pt05deg.nc', 'w', format='NETCDF4')
    f.description = f'Distance from ocean on 0.05deg CMG. Pixels with >0% forest cover loss only (Hansen 2013).'
    today = datetime.today()
    f.history = "Created " + today.strftime("%d/%m/%y")
    f.createDimension('lat', global_lat.size)
    f.createDimension('lon', global_lon.size)
    latitude = f.createVariable('lat', 'f4', 'lat')
    latitude[:] = global_lat
    latitude.standard_name = 'lat'
    latitude.long_name = 'latitude'
    latitude.units = 'degrees_north'
    longitude = f.createVariable('lon', 'f4', 'lon')
    longitude[:] = global_lon
    longitude.standard_name = 'lon'
    longitude.long_name = 'longitude'
    longitude.units = 'degrees_east'
    distance_km = f.createVariable(f'distance_from_ocean', 'f4', ('lat', 'lon'))
    distance_km.standard_name = 'distance_from_ocean'
    distance_km.long_name = 'distance_from_ocean'
    distance_km.units = 'km'
    distance_km[:] = distance_map
    f.close()


if __name__ == '__main__':
    neighbour_tree()
    write_distance_map()
