import cdsapi
import os

def download(year, month, file):
    c = cdsapi.Client()

    c.retrieve(
        'reanalysis-era5-single-levels',
        {
            'product_type':'reanalysis',
            'format':'netcdf',
            # 'variable':[
            #     '100m_u_component_of_wind','100m_v_component_of_wind','10m_u_component_of_wind',
            #     '10m_v_component_of_wind','2m_dewpoint_temperature','2m_temperature',
            #     'convective_available_potential_energy','convective_inhibition','mean_sea_level_pressure',
            #     'mean_surface_downward_long_wave_radiation_flux','mean_surface_downward_short_wave_radiation_flux','mean_surface_latent_heat_flux',
            #     'mean_surface_net_long_wave_radiation_flux','mean_surface_net_short_wave_radiation_flux','mean_surface_sensible_heat_flux',
            #     'mean_total_precipitation_rate','mean_vertically_integrated_moisture_divergence','sea_surface_temperature',
            #     'skin_temperature','surface_pressure','total_column_cloud_ice_water',
            #     'total_column_water','total_column_water_vapour','vertical_integral_of_divergence_of_mass_flux',
            #     'vertical_integral_of_divergence_of_moisture_flux','vertical_integral_of_divergence_of_total_energy_flux','vertical_integral_of_thermal_energy',
            #     'vertically_integrated_moisture_divergence'
            # ],
            'variable' : ['volumetric_soil_water_layer_1'],
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
            'time': ['00:00', '09:00', '10:00', '11:00','12:00', '18:00']#['00:00','03:00','06:00','09:00','12:00','15:00','18:00','21:00']
        },
        file)

for y in range(2006,2011): # (1979,2020)
    for m in range(1, 13):

        out_dir = '/prj/AMMA2050/ERA5/hourly/surface/'
        path_file =  out_dir + 'ERA5_' + str(y) + '_' + str(m).zfill(2) + '_srfc_SM.nc'
        print('Doing ' + path_file)

        if os.path.isfile(path_file):
            try:
                os.remove(path_file)
            except OSError:
                pass
            # print('File exists, continue')
            # continue

        download(y, m, path_file)
