import glob
import os
import subprocess
from JASMIN.constants import VARDIC
import numpy as np


class UM_vars(object):

    VARDIC = {
        'sw_net': ('01', '201', 'a', 'W m-2'),
        'sw_in': ('01', '235', 'a', 'W m-2'),
        'lw_net': ('02', '201', 'a', 'W m-2'),
        'lw_in': ('02', '207', 'a', 'W m-2'),
        'SH': ('03', '217', 'a', 'W m-2'),
        'LH': ('03', '234', 'a', 'W m-2'),
        'ConvRain': ('05', '205', 'a', 'kg m-2s-1'),
        'TotalRain': ('05', '216', 'a', 'kg m-2s-1'),
        'LargeScaleRain': ('04', '203', 'a', 'kg m-2s-1'),
        'lst': ('00', '024', 'c', 'K'),
        'q2': ('03', '237', 'c', 'kg kg-1'),
        'T2': ('03', '236', 'c', 'K'),
        'slp': ('16', '222', 'c', 'Pa'),
        'v10': ('03', '226', 'c', 'm s-1'),
        'u10': ('03', '225', 'c', 'm s-1'),
        'w_pl': ('30', '203', 'f', 'm s-1'),
        'u_pl': ('30', '201', 'f', 'm s-1'),
        'v_pl': ('30', '202', 'f', 'm s-1'),
        'T_pl': ('30', '204', 'f', 'K'),
        'q_pl': ('30', '205', 'f', 'kg kg-1'),
        'rh_pl': ('30', '206', 'f', '%'),
        'omega_pl': ('30', '208', 'f', 'Pa s-1'),
        'geoH_pl': ('16', '2', 'f', 'gpm')
    }

    def __init__(self, iterable=VARDIC, **kwargs):
        self.__dict__.update(iterable, **kwargs)



def create_vera_var(var, hourly=False):
    """
    Hourly option fluxes and rain only
    """
    name = 'STASH_m01s' + VARDIC[var][0] + 'i' + VARDIC[var][1]
    if hourly:
        name = name+'_2'
    return name

def create_CP4_filename(var):
    return VARDIC[var][2]+VARDIC[var][0]+VARDIC[var][1]


def read_CP4_filename(str):
    keys = VARDIC.keys()
    for k in keys:

        pos = VARDIC[k]
        if (pos[2]==str[0:1]) & (pos[0]==str[1:3]) & (pos[1]==str[3:6]):

            return k

    print('Var not found')
    return None


def link_folders(inpath, outpath):

    folders = glob.glob(inpath+'/*')

    for f in folders:
        if not os.path.isdir(f):
            continue

        var = read_CP4_filename(os.path.basename(f))

        if var==None:
            continue

        new_f = outpath + os.sep + var

        command = "ln -s "+ f + " " + new_f
        process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()

