import cdsapi

def download(year):
    c = cdsapi.Client()

    c.retrieve(
        'reanalysis-era5-single-levels',
        {
            'product_type': 'reanalysis',
            'format': 'netcdf',
            'variable': [
                'land_sea_mask', 'orography', 'slope_of_sub_gridscale_orography',
                'standard_deviation_of_orography',
            ],
            'area': '1/-82/-18/-57',
            'year': '2019',
            'month': '08',
            'day': '28',
            'time': '12:00',
        },   '/media/ck/Elements/SouthAmerica/ERA5/hourly/ERA5_static_hourly_0.25deg.nc')

for y in range(2019,2020):
    print('Doing year', y)
    download(y)