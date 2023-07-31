import cdsapi
from utils import constants as cnst

# @ct.application(title='Download data')
# @ct.output.download()

c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-single-levels',
    {
        'product_type': 'reanalysis',
        'format': 'netcdf',
        'variable': 'land_sea_mask',
        'day': [
            '01'
        ],
        'time': [
            '00:00'
        ],
        'year': '2011',
        'month': '06',
        'grid': '0.7/0.7',
        'area': '90/-180/-90/180',
    },
    cnst.lmcs_drive+'/ERA5_global_0.7/monthly/ERA5_monthly_0.7deg_static.nc')