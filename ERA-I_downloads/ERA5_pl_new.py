import cdsapi
import os

def download(year, month, file):
    c = cdsapi.Client()

    c.retrieve(
        'reanalysis-era5-pressure-levels',
        {
            'product_type': 'reanalysis',
            'format': 'netcdf',
             'variable': [
                 'divergence', 'geopotential', 'potential_vorticity',
                 'relative_humidity', 'specific_humidity', 'temperature',
                 'u_component_of_wind', 'v_component_of_wind', 'vertical_velocity'
             ],
            #'variable': [
            #        'divergence',
            #        'specific_humidity', 'temperature'
            #    ],
             'pressure_level': [
                 '250',
                 '350', '450',
                 '500', '550', '600',
                 '650', '700', '750', '800',
                 '825', '850', '875',
                 '900', '925', '950',
                 '975'
             ],
            #'pressure_level': ['925'],
            'area' : '25/-18.5/3.5/17',   # pick domain upper/left/lower/right
            'year': [str(year)],
            'month': [str(month).zfill(2)],
            'day': [
                '01', '02', '03',
                '04', '05', '06',
                '07', '08', '09',
                '10', '11', '12',
                '13', '14', '15',
                '16', '17', '18',
                '19', '20', '21',
                '22', '23', '24',
                '25', '26', '27',
                '28', '29', '30',
                '31'
            ],
            'time':  ['00:00','03:00','06:00','09:00','12:00','15:00','18:00','21:00']
        },  file)

for y in range(2012,2014): # (1979,2020)
    for m in range(1, 13):

        out_dir = '/prj/AMMA2050/ERA5/hourly/pressure_levels/'
        path_file =  out_dir + 'ERA5_' + str(y) + '_' + str(m).zfill(2) + '_pl.nc'
        print('Doing ' + path_file)

        if os.path.isfile(path_file):
            try:
                os.remove(path_file)
            except OSError:
                pass
            # print('File exists, continue')
            # continue

        download(y, m, path_file)
