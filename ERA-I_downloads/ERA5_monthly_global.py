import cdsapi
import os


# @ct.application(title='Download data')
# @ct.output.download()

def download(year):
    c = cdsapi.Client()

    c.retrieve(
            'reanalysis-era5-single-levels-monthly-means',
            {
                'product_type': 'monthly_averaged_reanalysis',
                'variable': [
                    '10m_u_component_of_wind', '10m_v_component_of_wind', 'orography',
                    'sea_surface_temperature', 'surface_pressure',
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
                'time': ['00:00'],
        },
        '/media/ck/Elements/global/ERA5/monthly/ERA5_monthly_0.7deg_'+str(year)+'.nc')

for y in range(2020,2021):
    if os.path.isfile('/media/ck/Elements/global/ERA5/monthly/ERA5_monthly_0.7deg_'+str(y)+'.nc'):
        continue
    download(y)
