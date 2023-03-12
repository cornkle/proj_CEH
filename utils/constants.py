import os

########## FILE PATHS
if os.path.isdir('/home/ck/DIR/'):
    lappi_drive = '/home/ck/DIR/'  # backup on google drive for full drive. This where own data goes
    network_data = lappi_drive + 'cornkle/' # 'shared' on CEH network
    ext_drive = '/media/ck/LStorage/' # most of global_water project goes here / external drive
    lmcs_drive = ext_drive + 'global_water/' # specifically lmcs stuff
    other_drive = ext_drive + 'global_water/other/' # data from other projects
    scratch = ext_drive # scratch does not exist locally, so save on external drive


else:
    if os.path.isdir('/prj/global_water/'):
        scratch = '/scratch/cornkle/'  # scratch in CEH network
        network_data = '/users/global/cornkle/shared/'  # private U drive , equivalent of /cornkle/ on lappi
        lmcs_drive = '/prj/global_water/'  # lmcs stuff
        ext_drive = lmcs_drive # no external drive on CEH network, links to global_water proj
        other_drive = 'prj/global_water/other/' # other project's big files / external drive locally


ANCILS = network_data + 'data/ancils_python/'

DATA = network_data+'data/'
FIGS = network_data+'figs/'


#West Africa ERA5
ERA5_HOURLY_SRFC = DATA + 'ERA5/hourly/surface/'
ERA5_HOURLY_PL = DATA +'ERA5/hourly/pressure_levels/'
ERA5 = DATA + 'ERA5/'
ERA5_MONTHLY_PL_SYNOP = DATA + 'ERA5/monthly/synoptic/pl_1979-2019_monthly_synop_07x07.nc'
ERA5_MONTHLY_SRFC_SYNOP = DATA + 'ERA5/monthly/synoptic/srfc_1979-2019_monthly_synop_07x07.nc'


LSTA_TOPO = ANCILS + 'lsta_corr_topo.nc'
LSTA_TOPO_OLD = ANCILS + 'lsta_topo.nc'
WA_TOPO_1MIN = ANCILS + 'gtopo_1min_afr.nc'
TOPO_1MIN = ANCILS + 'gtopo_1min.nc'
MSG5KM_TOPO = ANCILS + 'msg5km_topo.nc'
WA_TOPO_3KM = ANCILS + 'gtopo_3km_WA.nc'


MCS_POINTS_DOM = network_data + 'MCSfiles/blob_map_allscales_-50_JJAS_points_dominant.nc' #/users/global/cornkle/MCSfiles/blob_map_MCSs_-40-75000_JJAS_-50-points_dominant.nc'#'/users/global/cornkle/MCSfiles/blob_map_allscales_-50_JJAS_points_dominant.nc'
MCS_HOUR_DAILY = network_data + 'MCSfiles/blob_map_MCSs_-50_JJAS_gt15k_daily_14UTC.nc'
MCS_TMIN = network_data + 'MCSfiles/blob_map_JJAS_-70minT_15k.nc'
MCS_15K = network_data + 'MCSfiles/blob_map_MCSs_-50_JJAS_gt15k.nc'
MCS_ALL = network_data + 'MCSfiles/blob_map_MCSs_-50_JJAS.nc'
MCS_CENTRE70 = network_data + 'MCSfiles/blob_map_JJAS_-70CentreMass_15k.nc'
MCS_CENTRE70_SMALL = network_data + 'MCSfiles/blob_map_JJAS_-70CentreMass_1-5k.nc'
MCS_CENTRE70_GT5 = network_data + 'MCSfiles/blob_map_JJAS_-70CentreMass_GT5k.nc'
MCS_CENTRE40 = network_data + 'MCSfiles/blob_map_JJAS_-40CentreMass_25k.nc'
MCS_BACK = network_data + 'MCSfiles/blob_map_JJAS_25000-40_-70CentreMassBack.nc'

CP4_PATH= DATA + '/CP4/'
CP4_OLR = CP4_PATH + 'CLOVER/CP4hist/lw_out_PBLtop/lw_out_PBLtop_MAM_A1hr_mean_ah261_4km_200205040030-200205042330.nc'
CP4_RAIN = CP4_PATH + 'CLOVER/CP4hist/lsRain/lsRain_MAM_A1hr_mean_ah261_4km_200205040030-200205042330.nc'
CP4_LANDSEA = CP4_PATH + 'ANCILS/landseamask_ancil_4km.nc'
CP4_TOPO = CP4_PATH + 'ANCILS/orog_original_ac144_ancil_4km.nc'

CLOVER_SAVES = network_data + '/data/CLOVER/saves/'


########## DIRECTORIES

LSTA = DATA + 'OBS/MSG_LSTA/lsta_netcdf/'
LSTA_NEW = DATA + 'OBS/MSG_LSTA/lsta_netcdf_new/'  # lsta_netcdf_new
LSTA_1330 = DATA + 'OBS/MSG_LSTA/lsta_netcdf_1330/'
LSTA_NEW_MEAN = DATA + 'OBS/MSG_LSTA/lsta_netcdf_new_-mean/'
LSTA_TESTFILE = LSTA_NEW + 'lsta_daily_20060624.nc'
LSTA_DOMSCALE = network_data + 'OBS/MSG_LSTA/power_maps/'

AMSRE_DIR = DATA + 'OBS/AMSRE/aqua/'
AMSRE_NIGHT = AMSRE_DIR + 'nc_night/'
AMSRE_DAY = AMSRE_DIR + 'nc_day/'
AMSRE_ANO_NIGHT = AMSRE_DIR + 'sma_nc_night_new/'
AMSRE_ANO_DAY = AMSRE_DIR + 'sma_nc_day_new/'
AMSRE_ANO_DAY_CORR = AMSRE_DIR + 'sma_nc_day_corr/'
AMSRE_NIGHT_TESTFILE = AMSRE_NIGHT + 'AMSR_L3_LPRMv05_A_20060116.nc'


TRMM5KM = network_data + 'TRMMfiles/'
TRMM5KM_FILE = TRMM5KM + 'TRMM5km_2006-2010.nc'

CMORPH = DATA + 'OBS/CMORPH/CMORPH_nc/'

GRIDSAT = DATA + 'GRIDSAT/MCS18/'
GRIDSAT_PERU = DATA + 'GRIDSAT/MCS18_peru/'
GRIDSAT_RAW = DATA + 'GRIDSAT/'
CHIRPS = DATA + 'OBS/CHIRPS/'
CHIRPS_MONTHLY = DATA + 'OBS/CHIRPS/monthly/'



########### Lost in elements head crash:

# TRMM_DAILY = elements_drive + 'TRMM/data/daily/aggregated/'
# TRMM_MONTHLY = elements_drive + 'TRMM/data/monthly/aggregated/'

# ERA5_MONTHLY_PL_SYNOP_HU = elements_drive + 'SouthAmerica/ERA5/monthly/pressure_levels/synop'
# ERA5_HOURLY_PL_HU = elements_drive + 'SouthAmerica/ERA5/hourly/pressure_levels/'
# ERA5_HOURLY_SRFC_HU = elements_drive + 'SouthAmerica/ERA5/hourly/surface/'
