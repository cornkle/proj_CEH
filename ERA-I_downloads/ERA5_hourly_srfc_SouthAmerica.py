import cdsapi
import os

def download(year, month, day, file):
    c = cdsapi.Client()

    c.retrieve(
        'reanalysis-era5-single-levels',
        {
            'product_type':'reanalysis',
            'format':'netcdf',
            'variable':[
                '10m_u_component_of_wind',
                '10m_v_component_of_wind','2m_dewpoint_temperature','2m_temperature',
                'convective_available_potential_energy','convective_inhibition','mean_sea_level_pressure',
                'mean_surface_latent_heat_flux',
                'mean_surface_net_long_wave_radiation_flux','mean_surface_net_short_wave_radiation_flux','mean_surface_sensible_heat_flux',
                'mean_total_precipitation_rate',
                'surface_pressure', 'total_cloud_cover', 'total_column_cloud_ice_water',
                'total_column_water_vapour', 'vertical_integral_of_divergence_of_moisture_flux',

            ],
            'area' : '1/-82/-18/-57', #'25/-18.5/3.5/17',   # pick domain upper/left/lower/right
            'year': [str(year)],
            'month': [str(month).zfill(2)],
            'day': [str(day).zfill(2)],

            'time': ['00:00','03:00','06:00','09:00','12:00','15:00','18:00','21:00']
            #['00:00','01:00', '02:00', '03:00',
                     # '04:00', '05:00', '06:00', '07:00',
                     # '08:00', '09:00', '10:00', '11:00',
                     # '12:00','13:00', '14:00', '15:00',
                     # '16:00', '17:00', '18:00', '19:00',
                     # '20:00', '21:00', '22:00', '23:00'] #['00:00','03:00','06:00','09:00','12:00','15:00','18:00','21:00']
        },
        file)


mdic = {1: range(1, 32), 2: range(1, 30), 3: range(1, 32),
       4: range(1, 31), 5: range(1, 32), 6: range(1, 31),
       7: range(1, 32), 8: range(1, 32), 9: range(1, 31),
       10: range(1, 32), 11: range(1, 31), 12: range(1, 32)
       }

for y in range(2020,2022): # (1979,2020)
    for m in range(1, 13):
        for d in mdic[m]:

            out_dir = '/media/ck/Elements/SouthAmerica/ERA5/hourly/surface/'
            path_file =  out_dir + 'ERA5_' + str(y) + '_' + str(m).zfill(2) + '_' + str(d).zfill(2) + '_srfc.nc'
            print('Doing ' + path_file)

            if os.path.isfile(path_file):
                # try:
                #     os.remove(path_file)
                # except OSError:
                #     pass
                print('File exists, continue')
                continue
            try:
                download(y, m, d, path_file)
            except:
                continue

