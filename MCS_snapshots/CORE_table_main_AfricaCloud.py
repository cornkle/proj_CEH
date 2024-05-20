from MCS_snapshots import MCS_table_create_01deg
import xarray as xr
import ipdb
import pandas as pd
import glob
import os
from utils import constants as cnst
from GLOBAL import glob_util
import datetime
import numpy as np
import multiprocessing


MREGIONS = {
 'GPlains' : [[-100,-90,32,47], 'nam', -6, (1,7), (5,9), (1,12)], # # 18
 'china' : [[105,115,25,40], 'asia', 8 , (1,7), (5,9), (1,12)], # 4
 'india' : [[70,90, 5,30], 'asia', 5, (1,7), (5,9), (1,12)], # 7
 'WAf' : [[-18,25,4,25], 'spac', 0, (1,7), (5,9), (1,12)], # last is hourly offset to UCT # 12    # [-18,25,4,25]
 'australia' : [[120,140,-23, -11], 'asia', 9, (11,3), (11,3), (1,12)], # 3
 'SAf' : [[20,35, -35,-15], 'spac', 2, (9,12), (11,3), (11,3)], # 10
 'sub_SA' : [[-68,-47, -40, -20.5], 'spac', -4, (11,3), (11,3), (1,12)] , # 16
 'Africa' : [[-18, 52, -36,39], 'spac', 2, (9,12), (11,3), (11,3)],
 'SA_big' : [[9.5, 52, -36,-5.5], 'spac', 2, (9,12), (11,3), (11,3)],
}

#EAf [26,52.5,-5.6, 15.5]

def make_table(yy):
    """
    Start with scanning image for MCSs as defined in MCS_table_create.process_tir_image.
    :return:
    """

    reg = 'WAf'

    lmcs = '/prj/Africa_cloud/ch9_wavelet/'
    out = cnst.lmcs_drive + '/core_tables/' + reg + '/'
    box = MREGIONS[reg][0]

    lpath = glob.glob(cnst.network_data + 'data/nflics/geoloc/lat_lon_2268_2080.npz')[0]  # this is /prj/Africa_cloud/geoloc/*.npz on the Linux system\n",
    msg_latlon = np.load(lpath)
    dlon = msg_latlon['lon']  # [100:-100,100:-100]
    dlat = msg_latlon['lat']  # [100:-100,100:-100]

    infiles = sorted(glob.glob(lmcs + str(yy) + '*/*.nc'))
    full_year = []
    out_dic = {}
    print('Doing', yy)
    outfile = out + str(yy)+'_'+reg+'_CORE_15mins.csv'
    if os.path.isfile(outfile):
       print(outfile, ' exists, continue')
       return

    for infile in infiles:
        print('Doing', infile)
        try:
            da = (xr.open_dataset(infile))
        except:
            print('2d file open error, continue')
            continue

        da['dlat'] = dlat
        da['dlon'] = dlon
        da = da.where((dlon>=box[0]) & (dlon<=box[1]) & (dlat>=box[2]) & (dlat<=box[3]), other=0)
        merge_tab = CORE_table_create_AfricaCloud.process_core_image(da, 3)
        #merge_tab = MCS_table_create_01deg.add_environment_toTable(basic_tab, da,  envvar_take=[], rainvar_name='precipitation')

        merge_tab.pop('cloudMask')
        merge_tab.pop('tir')
        if len(merge_tab['date']) ==0:
           continue

        full_year.append(merge_tab)
        print('Did' , infile)

    for key in full_year[0].keys():
        out_dic[key] = []
    for single_tab in full_year:
        for key in single_tab.keys():
            out_dic[key].extend(single_tab[key])

    pd_out = pd.DataFrame.from_dict(out_dic)
    pd_out.to_csv(outfile)
    del out_dic
    del pd_out
    del merge_tab
    #del basic_tab
    del da


#for yy in np.arange(2004,2024):
#     make_table(reg)

# pool = multiprocessing.Pool(processes=5)
# res = pool.map(make_table, np.arange(2004,2024))
# pool.close()
