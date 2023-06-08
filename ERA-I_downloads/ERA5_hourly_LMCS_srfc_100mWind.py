import cdsapi
import os
from utils import constants as cnst
def download(year, month, day, box,file):
    c = cdsapi.Client()

    c.retrieve(
        'reanalysis-era5-single-levels',
        {
            'format':'netcdf',
            'product_type':'reanalysis',
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
            'variable': [
                '10m_u_component_of_wind',
                '10m_v_component_of_wind','2m_dewpoint_temperature','2m_temperature',
                'convective_available_potential_energy','convective_inhibition','mean_sea_level_pressure',
                'mean_surface_latent_heat_flux', 'volumetric_soil_water_layer_1',
                'mean_surface_net_short_wave_radiation_flux','mean_surface_sensible_heat_flux',
                'mean_total_precipitation_rate',
                'surface_pressure', 'total_cloud_cover', 'total_column_cloud_ice_water',
                'total_column_water_vapour', 'vertical_integral_of_divergence_of_moisture_flux',

            ],
            'year':[str(year) ],
            'day': [str(day)
            ],
            'month':[
                str(month)
            ], # [upper,left,lower,right]
            'area': str(box[3]+2)+'/'+str(box[0]-2)+'/'+str(box[2]-2)+'/'+str(box[1]+2)
        },
        file)

# for y in range(1979,2020):
#     download(y)



mregions = {'WAf' : [[-18,25,4,25], 'spac', 0], # last is hourly offset to UCT # 12
'SAf' : [[20,35, -35,-15], 'spac', 2], # 10
'india' : [[70,90, 5,30], 'asia', 5], # 7
'china' : [[105,115,25,40], 'asia', 8 ], # 4
'australia' : [[120,140,-23, -11], 'asia', 9], # 3
#'sub_SA' : [[-68,-47, -40, -20.5], 'spac', -4] , # 16
# 'trop_SA' : [[-75, -50, -20, -5], 'spac', -5], # 17
#'GPlains' : [[-100,-90,32,47], 'nam', -6] # # 18

}


mdays = {1:31, 2:28, 3:31, 4:30, 5:31,6:30, 7:31, 8:31, 9:30, 10:31, 11:30,12:31}

for mm in mregions.keys():
    mreg = mm
    box = mregions[mm][0]
    for y in range(2000,2020):
        for m in range(1, 13):
            for d in range(1,mdays[m]+1):

                filename = "ERA5_" + str(y) + "_" + str(m).zfill(2) +"_" + str(d).zfill(2) + "_"+ mreg + "_srfc.nc"

                path = cnst.lmcs_drive+"ERA5/hourly/surface/"+mreg+"/"

                if not os.path.isdir(path):
                    os.mkdir(path)

                if os.path.isfile(path+filename):
                    print('File exists, continue!')
                    continue


                download(y, m, d, box, path+filename)
