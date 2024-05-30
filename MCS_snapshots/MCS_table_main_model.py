from LMCS import MCS_table_create
import xarray as xr
import ipdb
import pandas as pd
import numpy as np

def olr_to_bt(olr):
    #Application of Stefan-Boltzmann law
    sigma = 5.670373e-8
    tf = (olr/sigma)**0.25
    #Convert from bb to empirical BT (degC) - Yang and Slingo, 2001
    a = 1.228
    b = -1.106e-3
    Tb = (-a + np.sqrt(a**2 + 4*b*tf))/(2*b)
    return Tb - 273.15


def make_table():
    """
    Start with scanning image for MCSs as defined in MCS_table_create.process_tir_image.
    :return:
    """
    infile = '/media/ck/LStorage/global_water/other/CP4/CP4_WestAfrica/CP4hist/lw_out_PBLtop/lw_out_PBLtop_fullPL_A1hr_mean_ad251_4km_199808260030-199808262330.nc'
    da = (xr.open_dataarray(infile)/100).isel(time=18) # dataset already in degC

    return MCS_table_create.process_tir_image(da, 5)


############
### MCS cutout
basic_tab = make_table()
vfile = '/media/ck/LStorage/global_water/other/CP4/CP4_WestAfrica/CP4hist/q2/q2_fullPL__A1hr_inst_ad251_4km_199808260100-199808270000.nc'
tab = MCS_table_create.add_environment_toTable(vfile, basic_tab, envvar_take=['q2'])
#############
### Adding environmental variables to MCS table
vfile = '/media/ck/LStorage/global_water/other/CP4/CP4_WestAfrica/CP4hist/lsRain/lsRain_fullPL__A1hr_mean_ad251_4km_199808260030-199808262330.nc'
tab = MCS_table_create.add_environment_toTable(vfile, tab, envvar_take=[], rainvar_name='lsRain')

outpath = '/media/ck/LStorage/global_water/'
outfile = 'MCS_table_test.csv'

tab.pop('cloudMask')
tab.pop('tir')
# for k in tab.keys():
#     print(len(tab[k]))
pd.DataFrame(tab).to_csv(outpath+outfile)



