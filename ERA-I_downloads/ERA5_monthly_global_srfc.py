import cdsapi
import os
from utils import constants as cnst
def download(year,file):
    c = cdsapi.Client()

    c.retrieve(
        'reanalysis-era5-single-levels-monthly-means',
        {
            'product_type': 'monthly_averaged_reanalysis',
            'variable': [
                '100m_u_component_of_wind', '100m_v_component_of_wind', '10m_u_component_of_wind',
                '10m_v_component_of_wind', '2m_dewpoint_temperature', '2m_temperature',
                'boundary_layer_height', 'convective_available_potential_energy', 'mean_sea_level_pressure',
                'sea_surface_temperature', 'surface_latent_heat_flux', 'surface_pressure',
                'surface_sensible_heat_flux', 'total_column_water_vapour', 'total_precipitation',
                'vertically_integrated_moisture_divergence', 'volumetric_soil_water_layer_1', 'volumetric_soil_water_layer_4',
            ],
            'year': [
                '2017', '2018', '2019',
                '2020', '2021',
            ],
            'month': [
                '01', '02', '03',
                '04', '05', '06',
                '07', '08', '09',
                '10', '11', '12',
            ],
            'time': '00:00',
            'format': 'netcdf',
        },
            file)



    for y in range(2017,2022):

                filename = "ERA5_" + str(y) + "_pl.nc"

                path = "/media/ck/Elements/global/ERA5/monthly_new/"

                if os.path.isfile(path+filename):
                    print('File exists, continue!')
                    continue

                download(y, path+filename)
