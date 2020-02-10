import cdsapi
import os
import multiprocessing

def download(year, month, day, file):
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
            'pressure_level': [
                '250',
                '350', '450',
                '500', '550', '600',
                '650', '700', '750', '800',
                '825', '850', '875',
                '900', '925', '950',
                '975'
            ],
            'area' : '1/-82/-18/-57', #'10/-60/-50/90', #'25/-18.5/3.5/17',   # pick domain upper/left/lower/right
            'year': [str(year)],
            'month': [str(month).zfill(2)],
            'day': [str(day).zfill(2)],
            'time': ['00:00','03:00','06:00','09:00','12:00','15:00','18:00','21:00']
            # ['00:00','01:00', '02:00', '03:00',
            #          '04:00', '05:00', '06:00', '07:00',
            #          '08:00', '09:00', '10:00', '11:00',
            #          '12:00','13:00', '14:00', '15:00',
            #          '16:00', '17:00', '18:00', '19:00',
            #          '20:00', '21:00', '22:00', '23:00'] #'time':  ['00:00','03:00','06:00','09:00','12:00','15:00','18:00','21:00']
        },  file)

#dates = []
for y in range(2007,2018): # (1979,2020)
    for m in range(1, 13):
        for d in range(1,32):

            out_dir = '/prj/nflics/ERA5_scratch/hourly/pressure_levels/'
            path_file =  out_dir + 'ERA5_' + str(y) + '_' + str(m).zfill(2) + '_' + str(d).zfill(2) + '_pl.nc'
            print('Doing ' + path_file)

            if os.path.isfile(path_file):

                # try:
                #     os.remove(path_file)
                # except OSError:
                #     pass
                print('File exists, continue')
                continue

 #           dates.append((y,m,d,path_file))

            # pool = multiprocessing.Pool(processes=5)
            #
            # res = pool.map(func, data)
            # pool.close()
            try:
                download(y, m, d, path_file)
            except:
                continue