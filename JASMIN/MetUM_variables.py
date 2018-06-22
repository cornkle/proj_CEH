import glob
import os
import subprocess
from JASMIN.constants import VARDIC


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

