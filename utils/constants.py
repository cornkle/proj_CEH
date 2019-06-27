import os

########## FILE PATHS
if os.path.isdir('/home/ck/DIR/'):
    ext_drive = '/home/ck/DIR/'
    local_data = ext_drive + 'mymachine/'
    network_data = ext_drive + 'cornkle/'
    elements_drive = '/media/ck/Elements/'

else:
    if os.path.isdir('/localscratch/wllf030/cornkle/'):
        local_data = '/localscratch/wllf030/cornkle/'
        network_data = '/users/global/cornkle/shared/'
        elements_drive = '/media/ck/Elements/Obs_Data/'
    else:
        ext_drive = '/media/ck/Seagate/DIR/'
        local_data = ext_drive + 'mymachine/'
        network_data = ext_drive + 'cornkle/'


ANCILS = network_data + 'data/ancils_python/'


ERA_MONTHLY_PL = local_data + 'ERA-I/monthly/monthly_1979-2017_pl_full.nc'
ERA_MONTHLY_SRFC = local_data + 'ERA-I/monthly/monthly_1979-2017_srfc_full.nc'
ERA_MONTHLY_PL_SYNOP = local_data + 'ERA-I/monthly/monthly_synop_1979-2017_pl_full.nc'
ERA_MONTHLY_SRFC_SYNOP = local_data + 'ERA-I/monthly/monthly_synop_1979-2017_srfc_full.nc'
ERA_DAILY_SURFACE =   local_data + 'ERA-I/daily_2006-2010_12UTCsrfc.nc'

ERA_DAILY_PL12UTC = local_data + 'ERA-I/daily_1983-2014_pl.nc'
# ERA_DAILY_PL_NIGHT = '/localscratch/wllf030/cornkle/ERA-I/daily_2006-2010_12UTCpl_night.nc'
# ERA_DAILY_SRFC_ANO = '/localscratch/wllf030/cornkle/ERA-I/daily_2006-2010_12UTCsrfc_anomaly.nc'
ERA5 = local_data + 'ERA5/'

LSTA_TOPO = network_data + 'data/ancils_python/lsta_corr_topo.nc'
WA_TOPO_1MIN = network_data + 'data/ancils_python/gtopo_1min_afr.nc'
TOPO_1MIN = network_data + 'data/ancils_python/gtopo_1min.nc'
MSG5KM_TOPO = network_data + 'data/ancils_python/msg5km_topo.nc'


MCS_POINTS_DOM = network_data + 'MCSfiles/blob_map_MCSs_-40-75000_JJAS_-50-points_dominant.nc' #/users/global/cornkle/MCSfiles/blob_map_MCSs_-40-75000_JJAS_-50-points_dominant.nc'#'/users/global/cornkle/MCSfiles/blob_map_allscales_-50_JJAS_points_dominant.nc'
MCS_TMIN = network_data + 'MCSfiles/blob_map_JJAS_-70minT_15k.nc'
MCS_15K = network_data + 'MCSfiles/blob_map_MCSs_-50_JJAS_gt15k.nc'
MCS_CENTRE70 = network_data + 'MCSfiles/blob_map_JJAS_-70CentreMass_15k.nc'
MCS_CENTRE70_SMALL = network_data + 'MCSfiles/blob_map_JJAS_-70CentreMass_1-5k.nc'
MCS_CENTRE70_GT5 = network_data + 'MCSfiles/blob_map_JJAS_-70CentreMass_GT5k.nc'
MCS_CENTRE40 = network_data + 'MCSfiles/blob_map_JJAS_-40CentreMass_25k.nc'
MCS_BACK = network_data + 'MCSfiles/blob_map_JJAS_25000-40_-70CentreMassBack.nc'

CP4_PATH= network_data + 'data/CP4/'
CP4_OLR = CP4_PATH + 'CLOVER/CP4hist/lw_out_PBLtop/lw_out_PBLtop_MAM_A1hr_mean_ah261_4km_200205040030-200205042330.nc'
CP4_RAIN = CP4_PATH + 'CLOVER/CP4hist/lsRain/lsRain_MAM_A1hr_mean_ah261_4km_200205040030-200205042330.nc'
CP4_LANDSEA = CP4_PATH + 'ANCILS/landseamask_ancil_4km.nc'
CP4_TOPO = CP4_PATH + 'ANCILS/orog_original_ac144_ancil_4km.nc'

CLOVER_SAVES = network_data + '/data/CLOVER/saves/'



########## DIRECTORIES

LSTA = network_data + 'data/OBS/MSG_LSTA/lsta_netcdf/'
LSTA_NEW = network_data + 'data/OBS/MSG_LSTA/lsta_netcdf_new/'
LSTA_NEW_MEAN = network_data + 'data/OBS/MSG_LSTA/lsta_netcdf_new_-mean/'
LSTA_TESTFILE = LSTA_NEW + 'lsta_daily_20060624.nc'
LSTA_DOMSCALE = network_data + 'data/OBS/MSG_LSTA/power_maps/'

AMSRE_DIR = network_data + 'data/OBS/AMSRE/aqua/'
AMSRE_NIGHT = AMSRE_DIR + 'nc_night/'
AMSRE_DAY = AMSRE_DIR + 'nc_day/'
AMSRE_ANO_NIGHT = AMSRE_DIR + 'sma_nc_night/'
AMSRE_ANO_DAY = AMSRE_DIR + 'sma_nc_day/'
AMSRE_NIGHT_TESTFILE = AMSRE_NIGHT + 'AMSR_L3_LPRMv05_A_20060116.nc'




TRMM5KM = network_data + 'TRMMfiles/'
TRMM5KM_FILE = TRMM5KM + 'TRMM5km_2006-2010.nc'

CMORPH = network_data + 'data/OBS/CMORPH/CMORPH_nc/'

GRIDSAT = local_data + 'GRIDSAT/MCS18/'
GRIDSAT_RAW = local_data + 'GRIDSAT/'
CHIRPS = network_data + 'data/OBS/CHIRPS/'

TRMM_DAILY = elements_drive + 'TRMM/data/daily/aggregated/'
TRMM_MONTHLY = elements_drive + 'TRMM/data/monthly/aggregated/'
