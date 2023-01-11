import cdsapi
import os

c = cdsapi.Client()
def download(y):
    c.retrieve(
        #'reanalysis-era5-pressure-levels',
        'reanalysis-era5-single-levels',
        {
            'product_type': 'reanalysis',
            'format': 'netcdf',
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
           # 'pressure_level': ['200'],
            'year': [
                str(y)
            ],
            'month': [
                '01', '02', '03',
                '04', '05', '06',
                '07', '08', '09',
                '10', '11', '12',
            ],
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
                '31',
            ],
            'time': ['15:00'],  #'06:00','09:00','12:00',
            # 'area': [
            #     6, -83, -40,
            #     -32,
            # ],
            'area': [
                6, -82, -25, # [upper,left,lower,right]
                -58,   # -58
            ],
        },
        '/media/ck/Elements/SouthAmerica/ERA5/hourly/qr_15UTC/rain_15UTC_'+str(y)+'_peru2.nc')


strct = {

    'u200', [']
}


for y in range(1985,2019):  # 89 missing

    for

    if os.path.isfile('/media/ck/Elements/SouthAmerica/ERA5/hourly/qr_15UTC/rain_15UTC_'+str(y)+'_peru2.nc'):
        print('File exists, continue!')
        continue

    download(y)
