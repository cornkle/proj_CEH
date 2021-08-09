import cdsapi
import os


# @ct.application(title='Download data')
# @ct.output.download()

def download(year):
    c = cdsapi.Client()

    c.retrieve(
            'reanalysis-era5-pressure-levels-monthly-means',
            {
                'product_type': 'monthly_averaged_reanalysis',
                'variable': [
                    'u', 'v', 'z',
                    't', 'q',
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
                'pressure_level': [
                    '200',
                    '500', '850'
                ],
                'format': 'netcdf',
                'grid': '0.7/0.7',
                'time': ['00:00'],
        },
        '/media/ck/Elements/global/ERA5/monthly/ERA5_monthly_0.7deg_uv_'+str(year)+'.nc')

for y in range(1979,2021):
    if os.path.isfile('/media/ck/Elements/global/ERA5/monthly/ERA5_monthly_0.7deg_uv_'+str(y)+'.nc'):
        continue
    download(y)
