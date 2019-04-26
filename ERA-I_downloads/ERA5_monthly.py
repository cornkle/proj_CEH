import cdsapi

c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-single-levels-monthly-means',
    {
        'format':'netcdf',
        'product_type':'reanalysis-monthly-means-of-daily-means',
        'variable':[
            '100m_u_component_of_wind','100m_v_component_of_wind','10m_u_component_of_wind',
            '10m_v_component_of_neutral_wind','10m_v_component_of_wind','2m_dewpoint_temperature',
            '2m_temperature','boundary_layer_height','convective_available_potential_energy',
            'convective_inhibition','high_cloud_cover','instantaneous_10m_wind_gust',
            'mean_sea_level_pressure','orography','sea_surface_temperature',
            'surface_latent_heat_flux','surface_pressure','surface_sensible_heat_flux',
            'total_column_water','total_precipitation','vertical_integral_of_divergence_of_moisture_flux',
            'vertical_integral_of_thermal_energy','vertically_integrated_moisture_divergence'
        ],
        'year':[
            '1979','1980','1981',
            '1982','1983','1984',
            '1985','1986','1987',
            '1988','1989','1990',
            '1991','1992','1993',
            '1994','1995','1996',
            '1997','1998','1999',
            '2000','2001','2002',
            '2003','2004','2005',
            '2006','2007','2008',
            '2009','2010','2011',
            '2012','2013','2014',
            '2015','2016','2017',
            '2018','2019'
        ],
        'month':[
            '01','02','03',
            '04','05','06',
            '07','08','09',
            '10','11','12'
        ]
    },
    'ERA5_monthly_1997-2019.nc')