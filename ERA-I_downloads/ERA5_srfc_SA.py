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
                'mean_surface_downward_long_wave_radiation_flux','mean_surface_downward_short_wave_radiation_flux','mean_surface_latent_heat_flux',
                'mean_surface_net_long_wave_radiation_flux','mean_surface_net_short_wave_radiation_flux','mean_surface_sensible_heat_flux',
                'mean_total_precipitation_rate','sea_surface_temperature',
                'surface_pressure', 'total_cloud_cover', 'total_column_cloud_ice_water',
                'total_column_water_vapour', 'vertical_integral_of_divergence_of_moisture_flux',
                'soil_temperature_level_1', 'volumetric_soil_water_layer_1'
            ],
            'area' : '0/-20/-40/60', #'25/-18.5/3.5/17',   # pick domain upper/left/lower/right
            'year': [str(year)],
            'month': [str(month).zfill(2)],
            'day': [str(day).zfill(2)],

            'time': ['00:00','01:00', '02:00', '03:00',
                     '04:00', '05:00', '06:00', '07:00',
                     '08:00', '09:00', '10:00', '11:00',
                     '12:00','13:00', '14:00', '15:00',
                     '16:00', '17:00', '18:00', '19:00',
                     '20:00', '21:00', '22:00', '23:00'] #['00:00','03:00','06:00','09:00','12:00','15:00','18:00','21:00']
        },
        file)

for y in range(1995,2000): # (1979,2020)
    for m in range(1, 13):
        for d in range(1,32):

            out_dir = '/prj/nflics/ERA5_scratch/hourly/surface/'
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

