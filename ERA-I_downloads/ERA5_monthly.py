import cdsapi

def download(year):
    c = cdsapi.Client()

    c.retrieve(
        'reanalysis-era5-single-levels-monthly-means',
        {
            'format':'netcdf',
            'product_type':'reanalysis-monthly-means-of-daily-means',
            'variable':[
                '10m_u_component_of_wind',
                '10m_v_component_of_wind','2m_dewpoint_temperature',
                '2m_temperature','convective_available_potential_energy',
                'convective_inhibition','high_cloud_cover',
                'mean_sea_level_pressure','sea_surface_temperature',
                'surface_latent_heat_flux','surface_pressure','surface_sensible_heat_flux',
                'total_column_water','total_precipitation',
                'vertically_integrated_moisture_divergence'
            ],
            'year':[str(year) ],
            'month':[
                '03',
                '04','05','06',
                '07','08','09',
                '10','11'
            ],
            'grid' : '0.7/0.7' ,
            'time' : ['00:00'],
        },
        '/prj/AMMA2050/ERA5/synop/srfc/ERA5_monthly_'+str(year)+'.nc')

for y in range(1979,2020):
    download(y)
