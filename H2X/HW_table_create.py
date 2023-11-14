import numpy as np
import datetime as dt
import xarray as xr
import os
import pandas as pd
import glob
import pickle as pkl
import ipdb
from scipy.ndimage import label


def dictionary():

    dic = {}
    vars = ['date', 'month', 'hour', 'minute', 'year', 'day', 'area', 'tmax',
            'minlon', 'minlat', 'maxlon', 'maxlat', 'clon', 'clat', 'tmaxlon', 'tmaxlat',
            'tmean', 'tp1', 'tp99', 'HW_ID']
           #, 'hwMask', 'hw']


    for v in vars:
        dic[v] = []
    return dic

def hw_define_wloop(array, thresh, min_area=None, max_area=None, minmax_area=None):
    """

    :param array: 3d input array, assuming input as X day rolling slices i.e. time dimension reflects minimum length of considered HW
    :param thresh: heat wave threshold
    :param min_area: minimum area of considered HW
    :param max_area: maximum area of considered HW
    :param minmax_area: tuple indicating only HW bigger than tuple[0] and smaller than tuple[1]
    :return: 2d array with labelled HW blobs
    """
    array[array <= thresh] = 0  # T threshold maskout
    array[np.isnan(array)] = 0  # set nans to 0

    labels, numL = label(array)

    u, inv = np.unique(labels, return_inverse=True)
    goodinds_size = u[u != 0]

    goodinds = []
    badinds = []

    for uni in goodinds_size:

        l_dummy = labels.copy()
        l_dummy[l_dummy!=uni] = 0
        dummy = np.mean(l_dummy, axis=0)

        if min_area != None:
            if np.sum(dummy == uni) >= min_area:
                goodinds.append(uni)
            else:
                badinds.append(uni)

        if max_area != None:
            if np.sum(dummy == uni) <= max_area:
                goodinds.append(uni)
            else:
                badinds.append(uni)

        if minmax_area != None:
            if (np.sum(dummy == uni) <= max_area) & (np.sum(dummy==uni) >= min_area):
                goodinds.append(uni)
            else:
                badinds.append(uni)

        del l_dummy

    if (min_area is not None) | (max_area is not None) | (minmax_area is not None):
        for b in badinds:
            pos = np.where(labels==b)
            labels[pos]=0

    return labels, np.array(goodinds)


def hw_define(array, thresh, min_area=None, max_area=None, minmax_area=None):
    """

    :param array: 3d input array, assuming input as X day rolling slices i.e. time dimension reflects minimum length of considered HW
    :param thresh: heat wave threshold
    :param min_area: minimum area of considered HW
    :param max_area: maximum area of considered HW
    :param minmax_area: tuple indicating only HW bigger than tuple[0] and smaller than tuple[1]
    :return: 2d array with labelled HW blobs
    """
    array[array <= thresh] = 0  # T threshold maskout
    array[np.isnan(array)] = 0  # set nans to 0

    labels, numL = label(array)
    same_count = np.sum(labels == labels[0,:,:], axis=0)

    labels = np.mean(labels, axis=0)
    labels[same_count != array.shape[0]] = 0

    # lmask = np.broadcast_to(same_count != array.shape[0], array.shape) for 3d masking
    # labels[lmask] = 0    # determine pixels affected over all time steps.

    u, inv = np.unique(labels, return_inverse=True)
    n = np.bincount(inv)

    goodinds = u[u != 0]

    if min_area != None:
        goodinds = u[(n >= min_area) & (u != 0)]
        badinds = u[n < min_area]

    if max_area != None:
        goodinds = u[(n <= max_area) & (u != 0)]
        badinds = u[n > max_area]

    if minmax_area != None:
        goodinds = u[(n <= minmax_area[1]) & (u != 0) & (n >= minmax_area[0])]
        badinds = u[(n > minmax_area[1]) | (n < minmax_area[0])]

    if (min_area is not None) | (max_area is not None) | (minmax_area is not None):
        for b in badinds:
            pos = np.where(labels == b)
            labels[pos] = 0

    return labels, goodinds


def process_hw_image(hw, data_res, t_thresh=35, min_hw_size=60, max_hw_size=None):
    """
    This function cuts out HWs. By default, a HW is defined as contiguous wet bulb temperature area at >= 35 degC over at least 60 km2 (3 CP4 pixels).
    :param hw: hw temperature image (in degC)
    :param data_res: spatial resolution of input image (approximately, in km - this defines how many pixel are needed to define a heat wave)
    :param t_thresh: temperature threshold (in degC) for considered contiguous HW - minimum threshold to avoid counting one-pixel spikes.
    :param min_mcs_size: minimum size of contiguous HW to be considered a HW (in km2)
    :return: dictionary with dates and HW characteristics from temperature information only.
    """

    min_pix_nb = min_hw_size / data_res**2

    if max_hw_size!=None:
        max_pix_nb = max_hw_size / data_res**2  # this is to capture satellite artefacts that come in large contiguous stripes.
    else:
        max_pix_nb = 5000000 # some large area to allow most upper bounds

    labels, goodinds = hw_define(hw.values, t_thresh, minmax_area=[min_pix_nb, max_pix_nb]) # 4.4x4.4km = 20km2 per pix in CP4. Min standard HW size = 4 contiguous CP4 pixels (80km2)
    dic = dictionary()

    for g in goodinds:

        if g==0:
            continue

        pos = np.where(labels==g)
        npos = np.where(labels!=g)
        datestr = str(hw['time.year'].values)+'-'+str(hw['time.month'].values).zfill(2)+'-'+str(hw['time.day'].values).zfill(2)+'_'+\
                      str(hw['time.hour'].values).zfill(2)+':'+str(hw['time.minute'].values).zfill(2)

        dic['date'].append(datestr)
        dic['month'].append(int(hw['time.month']))
        dic['hour'].append(int(hw['time.hour']))
        dic['year'].append(int(hw['time.year']))
        dic['day'].append(int(hw['time.day']))
        dic['minute'].append(int(hw['time.minute']))

        hw_obj = hw.copy()
        hw_obj.values[npos] = np.nan
        ipdb.set_trace()
        tmin_pos = np.nanargmin(hw_obj.values)
        tpos_2d = np.unravel_index(tmin_pos, hw_obj.shape)

        latmin = np.nanmin(hw.latitude.values[pos[0]])
        latmax = np.nanmax(hw.latitude.values[pos[0]])
        lonmin = np.nanmin(hw.longitude.values[pos[1]])
        lonmax = np.nanmax(hw.longitude.values[pos[1]])

        dic['area'].append(np.sum(np.isfinite(hw_obj.values)))
        dic['minlon'].append(lonmin)
        dic['minlat'].append(latmin)
        dic['maxlon'].append(lonmax)
        dic['maxlat'].append(latmax)
        dic['clon'].append(lonmin + (lonmax - lonmin)/2)  # hw centre points as per box
        dic['clat'].append(latmin + (latmax - latmin)/2) # hw centre points as per box
        dic['tmax'].append(np.nanmax(hw_obj))
        dic['tmaxlat'].append(hw.latitude[tpos_2d[0]].values)  # alternative hw centre point: tmax
        dic['tmaxlon'].append(hw.longitude[tpos_2d[1]].values) # alternative hw centre point: tmax
        dic['tmean'].append(np.nanmean(hw_obj))
        dic['tp1'].append(np.nanpercentile(hw_obj, 1)) # percentiles within object
        dic['tp99'].append(np.nanpercentile(hw_obj, 99)) # percentiles within object
        dic['HW_ID'].append(datestr + '_' + str(g)) # some unique ID including date
        #dic['hwMask'].append(labels==g)
        #dic['hw'].append(hw_obj.values) # full heatwave array, can only be saved in pickle not csv

    # for k in dic.keys():
    #     print(k, len(dic[k]))

    return dic



############
### Heat wave cutout

infile = '/media/ck/LStorage/global_water/other/CP4/CP4_WestAfrica/CP4hist/'
var = 't2'

cp4_files = sorted(glob.glob(infile + var + '/' + var + '*_200005*.nc'))  # needs more clever way to bring CP4 dates in order
chunks = []

for idx, ff in enumerate(cp4_files[1::]):

    da = xr.open_mfdataset(cp4_files[idx:idx+3], concat_dim="time", combine="nested", decode_times=False)

    basic_tab = process_hw_image(da[var].load()-273.15, 4.4)
    del da

    outpath = '/outpath/outpath/'
    outfile = 'test_table'
    pd.DataFrame(basic_tab).to_csv(outpath+outfile+'.csv')

    #pkl.dump(basic_tab, open(outpath+outfile+".p", "wb")) to alternatively save original dictionary including full arrays