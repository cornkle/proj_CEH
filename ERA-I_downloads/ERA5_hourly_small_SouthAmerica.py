import cdsapi
import os

c = cdsapi.Client()
def download(y):
    c.retrieve(
        'reanalysis-era5-pressure-levels',
        {
            'product_type': 'reanalysis',
            'format': 'netcdf',
            'variable': ['u', 'v'],
            'pressure_level': ['200'],
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
            'time': '15:00',
            # 'area': [
            #     6, -83, -40,
            #     -32,
            # ],
            'area': [
                6, -82, -25, # [upper,left,lower,right]
                -45,   # -58
            ],
        },
        '/media/ck/Elements/SouthAmerica/ERA5/hourly/qr_15UTC/uv_15UTC_'+str(y)+'_peru2.nc')



for y in range(2000,2021):  # 89 missing

    if os.path.isfile('/media/ck/Elements/SouthAmerica/ERA5/hourly/qr_15UTC/uv_15UTC_'+str(y)+'_peru2.nc'):
        print('File exists, continue!')
        continue

    download(y)
