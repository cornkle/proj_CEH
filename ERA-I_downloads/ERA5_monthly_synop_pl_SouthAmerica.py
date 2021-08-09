import cdsapi

def download(year):
    c = cdsapi.Client()

    c.retrieve(
        'reanalysis-era5-pressure-levels-monthly-means',
        {
            'format': 'netcdf',
            'product_type':'monthly_averaged_reanalysis_by_hour_of_day',
            'variable': [
                'divergence', 'geopotential',
                'relative_humidity', 'specific_humidity',
                'temperature', 'u_component_of_wind', 'v_component_of_wind',
                'vertical_velocity'
            ],
            'pressure_level': [
                '200',
                '300', '400',
                '500', '550', '600',
                '650', '700', '750', '800',
                '825', '850', '875',
                '900', '925', '950',
                '975'
            ],
            'year': [str(year)],
            'month': [
                '01', '02', '03',
                '04', '05', '06',
                '07', '08', '09',
                '10', '11', '12'
            ],
            'area': '6/-83/-40/-32',#'25/-18.5/3.5/17',  # pick domain upper/left/lower/right
            'grid': '0.7/0.7',
            'time': [
                '00:00', '03:00',
                '06:00',
                '09:00',
                '12:00', '15:00',
                '18:00',
                '21:00'
            ]
        },   '/media/ck/Elements/SouthAmerica/ERA5/monthly/pressure_levels/synop/ERA5_monthly_pl_'+str(y)+'_synop.nc')

for y in range(2020,2022):
    print('Doing year', y)
    download(y)