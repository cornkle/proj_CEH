import xesmf as xe
import xarray as xr

# regridding for unstructured to structured lat/lon. Attention, does not handle reprojections.
# Only from lat/lon regridding from one projection into same projection

def interpolation_weights(x, y, new_x, new_y, regridder_path=None, method='conservative'):
    """
    x : numpy array for longitudes (can be 1d or 2d)
    y : numpy array for latitudes (can be 1d or 2d)
    new_x : new longitudes (can be 1d or 2d)
    new_y : new latitudes (can be 1d or 2d)
    regridder_path : option filepath for regridder file to be saved. If chosen, an intermediate Netcdf file will be saved from which later regridding
                     from/to the exact same grids can be repeated without recomputing the regridder.
                     For example: "regridder_path = /my/path/to/file/mtg_onWGS84_2km_xesmf.nc"

    method: string. interpolation method can be "bilinear" or "conservative".
            When regridding from fine to much coarser, it is better to use conservative.
    """
    # rename coordinates to lat/lon so xesmf can detect coordinates
    ds_in = xr.Dataset({"lat": (["lat"], y),"lon": (["lon"], x)})
    ds_out = xr.Dataset({"lat": (["lat"], new_y),"lon": (["lon"], new_x),})

    if regridder_path is None:
        regridder = xe.Regridder(ds_in, ds_out, method=method)
    else:
        regridder = xe.Regridder(ds_in, ds_out, method=method, filename=regridder_path
    )

    return regridder


def interpolate_data(data, regridder, x_name='lon', y_name='lat'):
    """
    Function to use if regridder is loaded in the python session already as object "regridder".

    data : xarray data array containing the variable to be regridded
    regridder: regridding object from xesmf, either from active session or loaded from intermediate file.
    x_name : longitude coordinate name in "data" data array.
    y_name: latitude coordinate name in "data" data array.
    """

    ds_in = data.rename({y_name: 'lat', x_name: 'lon'})

    ds_regridded = regridder(ds_in)  # apply regridding

    return ds_regridded



def interpolate_data_fromFile(data, regridder_path, x_name='lon', y_name='lat'):
    """
    Function to use instead of "interpolate data" if regridder is loaded from file as provided in "regridder_path".

    data : xarray data array containing the variable to be regridded
    regridder_path: path to regridding object, loaded from intermediate file.
                    For example: "regridder_path = /my/path/to/file/mtg_onWGS84_2km_xesmf.nc"
    x_name : longitude coordinate name in "data" data array.
    y_name: latitude coordinate name in "data" data array.
    """

    # rename coordinates in case is needed
    data_for_xesmf = data.rename({y_name: 'lat', x_name: 'lon'})

    #Load the existing regridder weights
    regridder = xe.Regridder(
        None, None,
        method="bilinear",  # must match method used when weights were created
        filename=regridder_path,
        reuse_weights=True
    )

    # apply the regridder
    data_out = regridder(data_for_xesmf)
    return data_out


