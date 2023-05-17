import cdsapi
import os
from utils import constants as cnst
def download(year, month, day, box, file):
    c = cdsapi.Client()

    c.retrieve(
        'reanalysis-era5-pressure-levels',
        {
            'format': 'netcdf',
            'product_type': 'reanalysis',
            'pressure_level': [
                '925', '850', '650'
            ],
            # 'time': [
            #     str(h-3).zfill(2)+':00', str(h).zfill(2)+':00',
            #     str(h+3).zfill(2)+':00', str(h+6).zfill(2)+':00', str(h+9).zfill(2)+':00'
            # ],
            'time': [
                '10:00', '12:00', '18:00', ],
            'variable': [vv],
            'year': [str(year)],
            'day': [str(day)
                    ],
            'month': [
                str(month)
            ],  # [upper,left,lower,right]
            'grid': '0.7/0.7',
            'area': '90/-180/-90/180',
        },
        file)

    # for y in range(1979,2020):
    #     download(y)

mdays = {1: 31, 2: 28, 3: 31, 4: 30, 5: 31, 6: 30, 7: 31, 8: 31, 9: 30, 10: 31, 11: 30, 12: 31}
variable_list = [
            '10m_u_component_of_wind',
            '10m_v_component_of_wind','2m_dewpoint_temperature','2m_temperature',
            'convective_available_potential_energy','convective_inhibition','mean_sea_level_pressure',
            'mean_surface_latent_heat_flux', 'volumetric_soil_water_layer_1',
            'mean_surface_net_short_wave_radiation_flux','mean_surface_sensible_heat_flux',
            'mean_total_precipitation_rate',
            'surface_pressure', 'total_cloud_cover', 'total_column_cloud_ice_water',
            'total_column_water_vapour', 'vertical_integral_of_divergence_of_moisture_flux' ]

for vv in variable_list:
    for y in range(2020, 2021):
        for m in range(1, 13):
            for d in range(1, mdays[m] + 1):

                filename = vv + "_" + str(y) + "_" + str(m).zfill(2) + "_" + str(d).zfill(2) + "_pl.nc"

                path = cnst.lmcs_drive + "ERA5_global_0.7/hourly/pressure_levels/"

                if not os.path.isdir(path):
                    os.mkdir(path)

                if os.path.isfile(path + filename):
                    print('File exists, continue!')
                    continue

                download(y, m, d, vv, path + filename)