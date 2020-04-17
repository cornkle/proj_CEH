import cdsapi

def download(year):
    c = cdsapi.Client()

    c.retrieve(
        'reanalysis-era5-single-levels-monthly-means',
        {
            'format': 'netcdf',
            'variable': [
                'land_sea_mask', 'orography', 'slope_of_sub_gridscale_orography',
                'standard_deviation_of_orography'
            ],
            'year': '2019',
            'month': '08',
            'time': '12:00',
            'area': '6/-83/-40/-32',#'25/-18.5/3.5/17',  # pick domain upper/left/lower/right
            'grid': '0.7/0.7',

        },   '/media/ck/Elements/SouthAmerica/ERA5/monthly/pressure_levels/synop/ERA5_static_synop_0.7deg.nc')

for y in range(2019,2020):
    print('Doing year', y)
    download(y)