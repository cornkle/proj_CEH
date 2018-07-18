import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib
import pdb
import pandas as pd
from collections import OrderedDict
import salem
from utils import u_met, u_parallelise, u_gis, u_arrays, constants, u_grid


def cut_kernel(xpos, ypos, arr, date, lon, lat, t, parallax=False, rotate=False, probs=False):

    if parallax:
        km, coords = u_gis.call_parallax_era(date.month, t, lon, lat, 0, 0)
        lx, ly = km

        lx = int(np.round(lx / 3.))
        ly = int(np.round(ly / 3.))  # km into pixels
        xpos = xpos - lx
        ypos = ypos - ly

    dist = 200

    kernel = u_arrays.cut_kernel(arr,xpos, ypos,dist)

    vdic = {}

    for d in probs.data_vars:

        var = u_arrays.cut_kernel(probs[d].values,xpos, ypos,dist)
        vdic[d] = var

    cnt2 = np.zeros_like(vdic[list(vdic.keys())[0]])
    cnt2[np.isfinite(vdic[list(vdic.keys())[0]])] = 1

    if rotate:
        kernel = u_met.era_wind_rotate(kernel,date,lat,lon,level=700, ref_angle=90)

    if (np.sum(np.isfinite(kernel)) < 2):
        return

    kernel3 = kernel - np.nanmean(kernel)

    cnt = np.zeros_like(kernel)
    cnt[np.isfinite(kernel)] = 1



    if kernel.shape != (dist*2+1, dist*2+1):
        pdb.set_trace()



    return kernel, kernel3, cnt, cnt2, vdic



def get_previous_hours(date):


    tdic = {18 : ('36 hours', '15 hours'),
            19 : ('37 hours', '16 hours'),
            20: ('38 hours', '17 hours'),
            21: ('39 hours', '18 hours'),
            0: ('42 hours', '21 hours'),
            1: ('43 hours', '22 hours'),
            2: ('44 hours', '23 hours'),
            3: ('45 hours', '24 hours'),
            6: ('48 hours', '27 hours')}
    before = pd.Timedelta(tdic[date.hour][0])
    before2 = pd.Timedelta(tdic[date.hour][1])

    t1 = date - before
    t2 = date - before2

    file = constants.ERA5
    try:
        cmm = xr.open_dataset(file + 'ERA5_'+str(date.year)+'_12UTCpl.nc')
    except:
        return None
    cmm = cmm.sel(time=slice(t1, t2))
    cm = cmm['t'].sel(level=950).squeeze() * 1000

    cm = cm.to_dataset()

    shear =  (cmm['u'].sel(level=600).squeeze() - cmm['u'].sel(level=925).squeeze() ) #


    vwind_srfc = cmm['v'].sel(level=950).squeeze()
    uwind_srfc = cmm['u'].sel(level=950).squeeze()

    uwind = cmm['u'].sel(level=600).squeeze()
    t = cmm['t'].sel(level=950).squeeze()

    cm['shear'] = shear
    cm['uwind'] = uwind
    cm['u_srfc'] = uwind_srfc
    cm['v_srfc'] = vwind_srfc
    cm['t'] = t