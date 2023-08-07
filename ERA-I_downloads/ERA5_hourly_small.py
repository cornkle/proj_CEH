import cdsapi

def download(year, month, file):
    c = cdsapi.Client()

    c.retrieve(
        'reanalysis-era5-pressure-levels',
        {
            'format':'netcdf',
            'product_type':'reanalysis',
            'pressure_level': [
                '650',
                '925'
            ],
            'time': [
                '09:00', '12:00',
                '15:00', '18:00', '21:00'
            ],
            'variable': [
                'specific_humidity', 'u_component_of_wind', 'v_component_of_wind'
            ],
            'year':[str(year) ],
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
                '31'
            ],
            'month':[
                str(month)
            ],
            'area': '25/-18.5/3.5/17',
        },
        file)

# for y in range(1979,2020):
#     download(y)


for y in range(2000,2015):
    for m in range(1, 13):
        file = "/prj/AMMA2050/ERA5/pressure_levels_small/ERA5_" + str(y) + "_" + str(m).zfill(2) + "_pl_small.nc"

        download(y, m,file)