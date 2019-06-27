import cdsapi

def download(year):
    c = cdsapi.Client()

    c.retrieve(
        'reanalysis-era5-single-levels-monthly-means',
        {
            'format':'netcdf',
            'product_type':'monthly_averaged_reanalysis_by_hour_of_day',
            'variable':[
                '100m_u_component_of_wind','100m_v_component_of_wind','10m_u_component_of_wind',
                '10m_v_component_of_wind','2m_dewpoint_temperature',
                '2m_temperature','boundary_layer_height','convective_available_potential_energy',
                'convective_inhibition','high_cloud_cover',
                'mean_sea_level_pressure',
                'surface_latent_heat_flux','surface_pressure','surface_sensible_heat_flux',
                'total_column_water_vapour','vertical_integral_of_divergence_of_moisture_flux',
                'vertically_integrated_moisture_divergence'
            ],
            'year':[str(year)],
            'month':[
                '01','02','03',
                '04','05','06',
                '07','08','09',
                '10','11','12'
            ],
            'area': '25/-18.5/3.5/17',  # pick domain upper/left/lower/right
            'grid': '0.7/0.7',
            'time': [
                '00:00',  '03:00',
                 '06:00',
                '09:00',
                 '12:00', '15:00',
                '18:00',
                '21:00'
            ]
        },
        '/media/ck/Elements/ERA5/monthly/synoptic/surface/ERA5_monthly_srfc_'+str(y)+'_synop.nc')

for y in range(1979,2020):
    print('Doing year', y)
    download(y)