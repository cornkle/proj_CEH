import cdsapi
import os
from utils import constants as cnst

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
                    '500', '650','850', '925'
                ],
                'format': 'netcdf',
                'grid': '0.7/0.7',
                'area': '90/-180/-90/180',
                'time': ['00:00'],
        },
        cnst.lmcs_drive+'/ERA5_global_0.7/monthly/pressure_levels/ERA5_monthly_0.7deg_uv_'+str(year)+'.nc')

for y in range(2000,2023):
    if os.path.isfile(cnst.lmcs_drive+'/ERA5_global_0.7/monthly/pressure_levels/ERA5_monthly_0.7deg_'+str(y)+'.nc'):
        continue
    download(y)
