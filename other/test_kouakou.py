
import numpy as np
import xarray as xr
import array as arr
import ipdb

############## set directories ##################################################

dir_in_T_xojua ='/media/ck/Elements/kouadio/T/'

#dir_coord = '/media/kouadio/Seagate/VERADATA/VERA_4km/'

###################################################################################
####### set files to open #########################################################

file_T_xojua = dir_in_T_xojua+'xojua_30204_T_hi_*'

# cord_file_lat = dir_coord+'latitude_uv.nc'
# cord_file_lon = dir_coord+'longitude_uv.nc'

####################################################################################
########### open files  ############################################################
T_xojua = xr.open_mfdataset(file_T_xojua)

# lon_uv = xr.open_dataarray(cord_file_lon)
# lat_uv = xr.open_dataarray(cord_file_lat)


####################################################################################
###### replace coordinates name by time lat lon lev ################################

T_xojua  = T_xojua.rename({'grid_latitude_uv':'lat', 'grid_longitude_uv':'lon','TH1':'time','P_ECMWF':'lev'}).set_coords(['lon', 'lat','time','lev'])

lev = arr.array('d', [1000.,975.,950.,925.,850.,800.,750.,700.,650.,600.])
################################################################################################################
######## select time range 00H and 12H ########################################################3
T_xojua_00H = T_xojua['STASH_m01s30i204'][23::24,0:10,:,:]
T_xojua_12H = T_xojua['STASH_m01s30i204'][11::24,0:10,:,:]


T_xojua_00H = T_xojua_00H.where(T_xojua_00H != 0, other=np.nan)
T_xojua_12H = T_xojua_12H.where(T_xojua_12H != 0, other=np.nan)


##########################################################################################
############# compute time mean of the variables ########################################

T_xojua_00H_mean = T_xojua_00H.mean('time')

T_xojua_12H_mean = T_xojua_12H.mean('time')

###################################################################################
############## Computing potential temperature  ###################################
#lev_po = (np.true_divide(1000, lev))**0.286

def theta(pz, t):
#    ist = t + 273.15
     ist = t 

     try:
        ist =  ist * ((1000 / pz) ** 0.286)
           
     except ValueError:
        ist = (ist.T * ((1000 / pz) ** 0.286)).T

#    return  ist-273.15
#     print(ist.shape)
     return  ist

Theta_xojua_00H_mean = T_xojua_00H_mean.copy(deep=True)*0
Theta_xojua_12H_mean = T_xojua_12H_mean.copy(deep=True)*0



for i in range(0,10):

    Theta_xojua_00H_mean.values[i,:,:] = theta(lev[i],T_xojua_00H_mean.values[i,:,:])
    Theta_xojua_12H_mean.values[i,:,:] =  theta(lev[i],T_xojua_00H_mean.values[i,:,:])




print(Theta_xojua_00H_mean.values - T_xojua_00H_mean.values)

#print(Theta_xojua_00H_mean - T_xojua_00H_mean.values)

