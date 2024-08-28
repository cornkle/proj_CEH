from MCS_snapshots import MCS_table_create_model_LMCSruns as MCS_table_create
import xarray as xr
import ipdb
import pandas as pd
import numpy as np
from JASMIN import MetUM_variables as mu
import glob

def olr_to_bt(olr):
    #Application of Stefan-Boltzmann law
    sigma = 5.670373e-8
    tf = (olr/sigma)**0.25
    #Convert from bb to empirical BT (degC) - Yang and Slingo, 2001
    a = 1.228
    b = -1.106e-3
    Tb = (-a + np.sqrt(a**2 + 4*b*tf))/(2*b)
    return Tb - 273.15


def load_file(ffile,var):
    try:
        ds=xr.open_dataset(ffile)[var]
    except:
        #Deal with funny landseamask file
        ds=xr.open_dataset(ffile,decode_times=False)[var]
        ds=ds.rename({"rlat":"latitude","rlon":"longitude"})
    ds=ds.assign_coords({"longitude":(ds.longitude-360)})

    try:
        ds=ds.isel(pressure=slice(None,None,-1))
    except:
        pass
    return ds


def make_table():
    """
    Start with scanning image for MCSs as defined in MCS_table_create.process_tir_image.
    :return:
    """
    cp4_h = '/home/users/cornkle/linked_CP4/hist/'
    cp4_f = '/home/users/cornkle/linked_CP4/fut/'
    var = 'lw_out_PBLtop'
    #local = '/media/ck/LStorage/global_water/other/CP4/CP4_WestAfrica/CP4hist/'
    orig_var = mu.create_CP4_filename(var)

    year = [1998]
    month = [8]
    hour = [17]
    box = [-15,20,7,20]

    inrun = cp4_h
    
    for yy in year:
      for mm in month:
        infiles = glob.glob(inrun + '/'+var+'/'+orig_var+'_A1hr_mean_*_4km_'+str(yy)+str(mm).zfill(2)+'*.nc')
        for ff in infiles:
            da = load_file(ff, orig_var).isel(time=hour).sel(latitude=slice(box[2],box[3]), longitude=slice(box[0],box[1]))
            da = olr_to_bt(da) # convert olr to bt
            # Mask out ocean areas - only want continental MCSs
            #mask = load_file('/gws/nopw/j04/impala/shared/CP4A/ncfiles/4km/ANCILS/landseamask_ancil_4km.nc','lsm')
            #mask = mask[0,0,:,:].interp(latitude=da.latitude,longitude=da.longitude)
            #da = da.where(mask>0)
            
            return MCS_table_create.process_tir_image(da, 5)


############
### MCS cutout
basic_tab = make_table()
#vfile = '/media/ck/LStorage/global_water/other/CP4/CP4_WestAfrica/CP4hist/q2/q2_fullPL__A1hr_inst_ad251_4km_199808260100-199808270000.nc'
#tab = MCS_table_create.add_environment_toTable(vfile, basic_tab, envvar_take=['q2'])
#############
### Adding environmental variables to MCS table
#vfile = '/media/ck/LStorage/global_water/other/CP4/CP4_WestAfrica/CP4hist/lsRain/lsRain_fullPL__A1hr_mean_ad251_4km_199808260030-199808262330.nc'
#tab = MCS_table_create.add_environment_toTable(vfile, tab, envvar_take=[], rainvar_name='lsRain')

tab = basic_tab

outpath = '/home/users/cornkle/lmcs/cklein/CP_models/MCS_files/CP4_MCS_table_runOnJasmin/'
outfile = 'MCS_table_test.csv'

tab.pop('cloudMask')
tab.pop('tir')
# for k in tab.keys():
#     print(len(tab[k]))
pd.DataFrame(tab).to_csv(outpath+outfile)



