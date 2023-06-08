import cdsapi
import os
from utils import constants as cnst
def download(year, month, var,file):
    c = cdsapi.Client()

    c.retrieve(
        'reanalysis-era5-pressure-levels-monthly-means',
        {
            'format':'netcdf',
            'product_type':'monthly_averaged_reanalysis_by_hour_of_day',
            'pressure_level': [
                '200', '250', '300',
                '350', '400', '450',
                '500', '550', '600',
                '650', '700', '750',
                '825', '850', '875',
                '900', '925', '950',
                '975'
            ],
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
    'divergence', 'geopotential', 'potential_vorticity',
     'specific_humidity', 'temperature', 'relative_humidity',
    'u_component_of_wind', 'v_component_of_wind', 'vertical_velocity'
]


for vv in var:
    for y in range(2000, 2020):
        for m in range(1, 13):

                filename = vv + "_ERA5_" + str(y) + "_" + str(m).zfill(2) + "_pl.nc"
                path = cnst.lmcs_drive + "ERA5_global_0.7/monthly_synopticHour/pressure_levels/" + vv +"/"

                if not os.path.isdir(path):
                    os.mkdir(path)

                if os.path.isfile(path + filename):
                    print('File exists, continue!')
                    continue

                download(y, m, vv, path + filename)
