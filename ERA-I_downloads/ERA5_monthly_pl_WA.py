import cdsapi
import os

def download(year, month, file):
    c = cdsapi.Client()

    c.retrieve(
        'reanalysis-era5-pressure-levels-monthly-means',
        {
            'format': 'netcdf',
            'product_type': 'monthly_averaged_reanalysis',
            'variable': [
                'divergence', 'geopotential',
                'relative_humidity', 'specific_humidity',
                'temperature', 'u_component_of_wind', 'v_component_of_wind',
                'vertical_velocity'
            ],
            'pressure_level': [
                '200', '250', '300',
                '350', '400', '450',
                '500', '550', '600',
                '650', '700', '750',
                '825', '850', '875',
                '900', '925', '950',
                '975'
            ],
            'area' : '25/-18.5/3.5/17',   # pick domain upper/left/lower/right
            'year': [str(year)],
            'month': [str(month).zfill(2)],
            'time': '00:00'
        },  file)

for y in range(1979,2020):
    for m in range(1, 13):

        out_dir = '/home/ck/DIR/mymachine/ERA5/monthly/pressure_levels/'
        path_file =  out_dir + 'ERA5_' + str(y) + '_' + str(m).zfill(2) + '_pl.nc'
        print('Doing ' + path_file)

        if os.path.isfile(path_file):
            print('File exists, continue')
            continue

        download(y, m, path_file)