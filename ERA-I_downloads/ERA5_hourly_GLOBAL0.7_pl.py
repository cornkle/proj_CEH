import cdsapi
import os
from utils import constants as cnst
def download(year, month, day, var,file):
    c = cdsapi.Client()

    c.retrieve(
        'reanalysis-era5-pressure-levels',
        {
            'format':'netcdf',
            'product_type':'reanalysis',
            'pressure_level': [
                '925' , '850', '750','650', '500'
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
            'day': [str(day)
            ],
            'month':[
                str(month)
            ], # [upper,left,lower,right]
            'grid': '0.7/0.7',
            'area': "70/-139/-60/160"
        },
        file)

# for y in range(1979,2020):
#     download(y)

var = [
    'divergence', 'geopotential', 'potential_vorticity',
     'specific_humidity', 'temperature',
    'u_component_of_wind', 'v_component_of_wind', 'vertical_velocity'
]

mdays = {1: 31, 2: 28, 3: 31, 4: 30, 5: 31, 6: 30, 7: 31, 8: 31, 9: 30, 10: 31, 11: 30, 12: 31}

for vv in var:
    for y in range(2000, 2020):
        for m in range(1, 13):
            for d in range(1, mdays[m] + 1):

                filename = vv + "_ERA5_" + str(y) + "_" + str(m).zfill(2) + "_" + str(d).zfill(2) + "_pl.nc"
                path = cnst.lmcs_drive + "ERA5_global_0.7/hourly/pressure_levels/" + vv +"/"

                if not os.path.isdir(path):
                    os.mkdir(path)

                if os.path.isfile(path + filename):
                    print('File exists, continue!')
                    continue

                download(y, m, d, vv, path + filename)
