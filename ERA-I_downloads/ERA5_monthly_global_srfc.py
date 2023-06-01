import cdsapi
import os
from utils import constants as cnst

# @ct.application(title='Download data')
# @ct.output.download()

def download(year):
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
                   str(year),
                ],
                'month': [
                    '01', '02', '03',
                    '04', '05', '06',
                    '07', '08', '09',
                    '10', '11', '12',
                ],
                'format': 'netcdf',
                'grid': '0.7/0.7',
                'area': '90/-180/-90/180',
                'time': ['00:00'],
        },
        cnst.lmcs_drive+'/ERA5_global_0.7/monthly/surface/ERA5_monthly_0.7deg_'+str(year)+'.nc')

for y in range(2000,2023):
    if os.path.isfile(cnst.lmcs_drive+'/ERA5_global_0.7/monthly/surface/ERA5_monthly_0.7deg_'+str(y)+'.nc'):
        continue
    download(y)
