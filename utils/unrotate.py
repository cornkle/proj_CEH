import xarray as xr

LATLON_ATTRIBUTES = {
    'lat': {
        'standard_name': 'latitude',
        'long_name': 'latitude',
        'units': 'degrees_north'},
    'lon': {
        'standard_name': 'longitude',
        'long_name': 'longitude',
        'units': 'degrees_east'}}

TO_RENAME = {
    'rlat': 'lat',
    'rlon': 'lon'}

TIME_ENCODING = {
    'units':  'seconds since 1970-01-01 00:00:00',
    'calendar': 'standard',
    'dtype': 'float64'}


def unrotate_ds(ds, shift_lon=False):
    """Convert Dataset from roated-pole (lon only) to regular lat/lon"""
    try:
        rp = ds['rotated_pole']
    except ValueError:
        raise ValueError('Dataset contains not variable `rotated_pole`.')
    assert rp.attrs['grid_mapping_name'] == 'rotated_latitude_longitude'
    assert rp.attrs['grid_north_pole_latitude'] == 90
    lon_offset = rp.attrs['grid_north_pole_longitude']

    # rotate longitudes
    dsnew = ds.rename(TO_RENAME).drop('rotated_pole')
    if shift_lon:
        dsnew.coords['lon'] = ds['rlon'] + lon_offset
        dsnew = dsnew.swap_dims({'rlon': 'lon'}).drop('rlon')

    # modify attributes
    for da in dsnew.data_vars.values():
        da.attrs.pop('grid_mapping')
    for varn in ['lat', 'lon']:
        dsnew[varn].attrs = LATLON_ATTRIBUTES[varn]
    dsnew.attrs = {'Conventions': 'CF-1.7'}
    return dsnew


def unrotate_netcdf(infile, outfile, shift_lon=False):
    """Convert netCDF file from roated-pole (lon only) to regular lat/lon"""
    ds = xr.open_dataset(infile)
    dsnew = unrotate_ds(ds, shift_lon=shift_lon)
    encoding = {}
    if 'time' in dsnew:
        encoding.update(time=TIME_ENCODING)
    dsnew.to_netcdf(outfile, encoding=encoding)