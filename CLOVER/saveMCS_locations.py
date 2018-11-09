import numpy as np
from scipy.ndimage.measurements import label
import xarray as xr
import os
import ipdb
import glob
from scipy.interpolate import griddata
import pandas as pd
import pdb
import itertools


def olr_to_bt(olr):
    sigma = 5.670373e-8
    return ((olr/sigma)**0.25)-273.15


def file_save(cp_dir, out_dir, ancils_dir, datestring, tthresh, storm_dic, box=None):

    v = 'lw_out_PBLtop'

    #load seamask
    landsea_path = glob.glob(ancils_dir+os.sep+'landseamask*.nc')[0]

    landsea = xr.open_dataset(landsea_path, decode_times=False)
    ls = landsea['lsm']

    ls.rlon.values = ls.rlon.values-360
    if box is not None:
        ls = ls.sel(rlon=slice(box[0], box[1]), rlat=slice(box[2], box[3]))
    pos = np.where(ls[0, 0, :, :] == 0)
    lons, lats = np.meshgrid(ls.rlon.values, ls.rlat.values)
    goodinds = 0
    try:
        filepath = glob.glob(cp_dir+os.sep+str(v)+os.sep+'*'+datestring+'*.nc')[0]

    except IndexError:
        print('No file found, return')
        return

    arr = xr.open_dataset(filepath)

    if box is not None:
        darr = arr[v].sel(longitude=slice(box[0],box[1]), latitude=slice(box[2],box[3]))
    else:
        darr = arr[v]

    for dar in darr:

        dar.values[pos[0], pos[1]] = np.nan  # mask sea


        dar.values = olr_to_bt(dar.values)
        dar.values[dar.values >= tthresh] = 0  # T threshold maskout
        dar.values[np.isnan(dar.values)] = 0 # set ocean nans to 0

        try:
            date = dar.time.values[0]
        except IndexError:
            date = dar.time.values

        labels, numL = label(dar.values)

        u, inv = np.unique(labels, return_inverse=True)
        n = np.bincount(inv)

        goodinds = u[n >= 52]  # defines minimum MCS size e.g. 258 pix at 4.4km is 5000km2, 52 pix is 1000km2
        if not sum(goodinds) > 0:
            return

        for gi in goodinds:
            if (gi == 0):  # index 0 is always background, ignore!
                continue
            inds = np.where(labels == gi)
            #mask = np.where(labels!=gi)

            if np.nanmin(dar.values[inds]) > -70:
                continue

            indir = {}
            pdb.set_trace()
            index = np.datetime_as_string(dar.time.values)[0:-13] + '_' + str(gi)

            indir['hour'] =  int(dar['time.hour'])
            indir['month'] = int(dar['time.month'])
            indir['year'] = int(dar['time.year'])
            indir['date'] = date

            # cut a box for every single blob from msg - get min max lat lon of the blob
            latmax, latmin = np.nanmax(lats[inds]), np.nanmin(lats[inds])
            lonmax, lonmin = np.nanmax(lons[inds]), np.nanmin(lons[inds])

            indir['latmid'] = latmin + (latmax-latmin)/2
            indir['lonmid'] = lonmin + (lonmax-lonmin)/2
            indir['latmax'] = latmax
            indir['latmin'] = latmin
            indir['lonmax'] = lonmax
            indir['lonmin'] = lonmin
            indir['lons'] = lons[inds]
            indir['lats'] = lats[inds]
            indir['pixel'] = len(inds[0])

            storm_dic[index] = indir

            print('Wrote ' + index)

def run():
    data_path = '/users/global/cornkle/shared/data/CP4/CLOVER/CP4hist'  # CP4 data directory
    ancils_path = '/users/global/cornkle/shared/data/CP4/ANCILS' # directory with seamask file inside
    out_path = '/users/global/cornkle/shared/data/CP4/CLOVER/test'  # out directory to save loc files
    box = [-13, 13, 4.5, 9]  # W- E , S - N geographical coordinates box

    years = np.array(np.arange(1998,2007), dtype=str)
    months = ['03', '04', '05', '06', '07', '08', '09', '10'] #np.array(['10'])
    days = np.array(np.arange(1,32), dtype=str)

    tthresh = -50 # chosen temperature threshold, e.g. -50, -60, -70

    datelist = []
    for y,m,d in itertools.product(years, months, days):
        datelist.append(y+m+str(d).zfill(2))

    storm_dic = {}

    for d in datelist:
        file_save(data_path, out_path, ancils_path, d, tthresh, storm_dic, box=box)

    return storm_dic
