import cdsapi
import os

c = cdsapi.Client()
def download(input, y):
    c.retrieve(

        #'reanalysis-era5-pressure-levels',
        #'reanalysis-era5-single-levels',
        input[2],
        {
            'product_type': 'reanalysis',
            'format': 'netcdf',
            'variable': [
                input[0]
            ],
            'pressure_level': input[1],
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
            'time': ['15:00'],  #'06:00','09:00','12:00',
            # 'area': [
            #     6, -83, -40,
            #     -32,
            # ],
            'area': [
                6, -82, -25, # [upper,left,lower,right]
                -58,   # -58
            ],
        },
        '/media/ck/LStorage/SouthAmerica/ERA5/hourly/15UTC/'+st+'_15UTC_'+str(y)+'_peru.nc')


  #
  # 'divergence', 'geopotential', 'potential_vorticity',
  #               'relative_humidity', 'specific_humidity', 'temperature',
  #               'u_component_of_wind', 'v_component_of_wind', 'vertical_velocity'
strct = {
    'tcwv': ['total_column_water_vapour', [850], 'reanalysis-era5-single-levels'],
    'u200': ['u_component_of_wind', [200],'reanalysis-era5-pressure-levels'],
    #'u550', ['u_component_of_wind', 550],
    #'u850', ['u_component_of_wind', 850],
    'v200':['v_component_of_wind', [200],'reanalysis-era5-pressure-levels'],
    #'v550', ['v_component_of_wind', 550],
    #'v850', ['v_component_of_wind', 850],
    #'w650', ['w_component_of_wind', 650],
    #'v550', ['v_component_of_wind', 550],
    'v850': ['v_component_of_wind', [850],'reanalysis-era5-pressure-levels'],
    'u850': ['u_component_of_wind', [850],'reanalysis-era5-pressure-levels'],
    'q850': ['specific_humidity', [850], 'reanalysis-era5-pressure-levels']
    #'rh850', ['relative_humidity', 850],


}


#for y in range(1985,2019):  # 89 missing

for st in strct.keys():
    input = strct[st]
    for y in range(1985,2022):
        if os.path.isfile('/media/ck/LStorage/SouthAmerica/ERA5/hourly/15UTC/'+st+'_15UTC_'+str(y)+'_peru.nc'):
            print('File exists, continue!')
            continue

        download(input, y)
