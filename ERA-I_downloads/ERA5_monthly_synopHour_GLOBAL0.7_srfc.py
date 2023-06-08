import cdsapi
import os
from utils import constants as cnst
def download(year, month, var,file):
    c = cdsapi.Client()

    c.retrieve(
        'reanalysis-era5-single-levels-monthly-means',
        {
            'format':'netcdf',
            'product_type':'monthly_averaged_reanalysis_by_hour_of_day',
            # 'time': [
            #     str(h-3).zfill(2)+':00', str(h).zfill(2)+':00',
            #     str(h+3).zfill(2)+':00', str(h+6).zfill(2)+':00', str(h+9).zfill(2)+':00'
            # ],
            'time': ['00:00', '01:00', '02:00', '03:00',
                     '04:00', '05:00', '06:00', '07:00',
                     '08:00', '09:00', '10:00', '11:00',
                     '12:00', '13:00', '14:00', '15:00',
                     '16:00', '17:00', '18:00', '19:00',
                     '20:00', '21:00', '22:00', '23:00'],
            'variable': [var

            ],
            'year':[str(year) ],

            'month':[
                str(month)
            ], # [upper,left,lower,right]
            'grid': '0.7/0.7',
            'area': "70/-160/-60/160"
        },
        file)

# for y in range(1979,2020):
#     download(y)

var = [
                '10m_u_component_of_wind', '100m_u_component_of_wind', '100m_v_component_of_wind',
                '10m_v_component_of_wind','2m_dewpoint_temperature','2m_temperature',
                'convective_available_potential_energy','convective_inhibition',
                'mean_surface_latent_heat_flux', 'volumetric_soil_water_layer_1',
                'mean_surface_net_short_wave_radiation_flux','mean_surface_sensible_heat_flux',
                'mean_total_precipitation_rate',
                'surface_pressure', 'total_column_cloud_ice_water',
                'total_column_water_vapour', 'vertical_integral_of_divergence_of_moisture_flux',

            ]



for vv in var:
    for y in range(2000,2020):
        for m in range(1, 13):

                filename = vv+"_ERA5_" + str(y) + "_" + str(m).zfill(2) + "_srfc.nc"
                path = cnst.lmcs_drive + "ERA5_global_0.7/monthly_synopticHour/surface/"+vv +"/"

                if not os.path.isdir(path):
                    os.mkdir(path)

                if os.path.isfile(path+filename):
                    print('File exists, continue!')
                    continue


                download(y, m, vv, path+filename)
