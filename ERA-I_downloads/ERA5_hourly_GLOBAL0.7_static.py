import cdsapi

c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-single-levels',
    {
        'product_type': 'reanalysis',
        'variable': [
            'land_sea_mask',
            'geopotential',
            'standard_deviation_of_filtered_subgrid_orography',
            'slope_of_subgridscale_orography'
        ],
        'year': '2020',
        'month': '01',
        'day': '01',
        'time': '00:00',
        'format': 'netcdf',
        'grid': [0.25, 0.25],
    },
    '/prj/global_water/ERA5_global_0.7/static/era5_invariant_025deg.nc'
)