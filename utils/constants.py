########## FILE PATHS

ERA_MONTHLY_PL = '/localscratch/wllf030/cornkle/ERA-I/monthly/monthly_1979-2017_pl_full.nc'
ERA_MONTHLY_SRFC = '/localscratch/wllf030/cornkle/ERA-I/monthly/monthly_1979-2017_srfc_full.nc'
ERA_MONTHLY_PL_SYNOP = '/localscratch/wllf030/cornkle/ERA-I/monthly/monthly_synop_1979-2017_pl_full.nc'
ERA_MONTHLY_SRFC_SYNOP = '/localscratch/wllf030/cornkle/ERA-I/monthly/monthly_synop_1979-2017_srfc_full.nc'
ERA_DAILY_SURFACE =   '/localscratch/wllf030/cornkle/ERA-I/daily_2006-2010_12UTCsrfc.nc'
ERA_DAILY_PL = '/localscratch/wllf030/cornkle/ERA-I/daily_2006-2010_12UTCpl.nc'
ERA_DAILY_PL_NIGHT = '/localscratch/wllf030/cornkle/ERA-I/daily_2006-2010_12UTCpl_night.nc'
ERA_DAILY_SRFC_ANO = '/localscratch/wllf030/cornkle/ERA-I/daily_2006-2010_12UTCsrfc_anomaly.nc'
ERA5 = '/localscratch/wllf030/cornkle/ERA5/'



LSTA_TOPO = '/users/global/cornkle/data/pythonWorkspace/proj_CEH/topo/lsta_corr_topo.nc'
WA_TOPO_1MIN = '/users/global/cornkle/data/pythonWorkspace/proj_CEH/topo/gtopo_1min_afr.nc'
TOPO_1MIN = '/users/global/cornkle/data/pythonWorkspace/proj_CEH/topo/gtopo_1min.nc'
MSG5KM_TOPO = '/users/global/cornkle/data/pythonWorkspace/proj_CEH/topo/msg5km_topo.nc'


MCS_POINTS_DOM = '/users/global/cornkle/MCSfiles/blob_map_MCSs_-40-75000_JJAS_-50-points_dominant.nc' #/users/global/cornkle/MCSfiles/blob_map_MCSs_-40-75000_JJAS_-50-points_dominant.nc'#'/users/global/cornkle/MCSfiles/blob_map_allscales_-50_JJAS_points_dominant.nc'
MCS_TMIN = '/users/global/cornkle/MCSfiles/blob_map_JJAS_-70minT_15k.nc'
MCS_15K = '/users/global/cornkle/MCSfiles/blob_map_MCSs_-50_JJAS_gt15k.nc'
MCS_CENTRE70 = '/users/global/cornkle/MCSfiles/blob_map_JJAS_-70CentreMass_15k.nc'
MCS_CENTRE70_SMALL = '/users/global/cornkle/MCSfiles/blob_map_JJAS_-70CentreMass_1-5k.nc'
MCS_CENTRE70_GT5 = '/users/global/cornkle/MCSfiles/blob_map_JJAS_-70CentreMass_GT5k.nc'
MCS_CENTRE40 = '/users/global/cornkle/MCSfiles/blob_map_JJAS_-40CentreMass_25k.nc'
MCS_BACK = '/users/global/cornkle/MCSfiles/blob_map_JJAS_25000-40_-70CentreMassBack.nc'

CP4_PATH= '/users/global/cornkle/data/CP4/'
CP4_OLR = CP4_PATH + 'CLOVER/CP4hist/lw_out_PBLtop/lw_out_PBLtop_MAM_A1hr_mean_ah261_4km_200205040030-200205042330.nc'
CP4_RAIN = CP4_PATH + 'CLOVER/CP4hist/lsRain/lsRain_MAM_A1hr_mean_ah261_4km_200205040030-200205042330.nc'
CP4_LANDSEA = CP4_PATH + 'ANCILS/landseamask_ancil_4km.nc'
CP4_TOPO = CP4_PATH + 'ANCILS/orog_original_ac144_ancil_4km.nc'



########## DIRECTORIES

LSTA = '/users/global/cornkle/data/OBS/MSG_LSTA/lsta_netcdf/'
LSTA_NEW = '/users/global/cornkle/data/OBS/MSG_LSTA/lsta_netcdf_new/'
LSTA_NEW_MEAN = '/users/global/cornkle/data/OBS/MSG_LSTA/lsta_netcdf_new_-mean/'
LSTA_TESTFILE = LSTA_NEW + 'lsta_daily_20060624.nc'
LSTA_DOMSCALE = '/users/global/cornkle/data/OBS/MSG_LSTA/power_maps/'

AMSRE_DIR = '/users/global/cornkle/data/OBS/AMSRE/aqua/'
AMSRE_NIGHT = AMSRE_DIR + 'nc_night/'
AMSRE_ANO_NIGHT = AMSRE_DIR + 'sma_nc_night/'
AMSRE_ANO_DAY = AMSRE_DIR + 'sma_nc_day/'
AMSRE_NIGHT_TESTFILE = AMSRE_NIGHT + 'AMSR_L3_LPRMv05_A_20060116.nc'

TRMM5KM = '/users/global/cornkle/TRMMfiles/'
TRMM5KM_FILE = TRMM5KM + 'TRMM5km_2006-2010.nc'

CMORPH = '/users/global/cornkle/data/OBS/CMORPH/CMORPH_nc/'

GRIDSAT = '/users/global/cornkle/mymachine/GRIDSAT/MCS18/'