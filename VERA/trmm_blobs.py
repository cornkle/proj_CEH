import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from utils import u_arrays as ua
from collections import OrderedDict
import pandas as pd
import multiprocessing
import pickle as pkl
from scipy.ndimage.measurements import label
import pdb

def run():
    files = ua.locate(".nc", '/users/global/cornkle/VERA')
    pool = multiprocessing.Pool(processes=4)
    out = '/users/global/cornkle/VERA/blobs/'


    res = pool.map(file_loop, files)
    pool.close()
    res = [x for x in res if x is not None]

    nb_sys = len(res)

    print('Number systems: ', nb_sys)

    res = [item for sublist in res for item in sublist]  # flatten list of lists


    # blon, blat, blon_c, blat_c, p ,pmax, pmean

    dic = OrderedDict([('lons', []) , ('lats', []), ('lon_c', []), ('lat_c', []),
                       ('p', []), ('pmax', []), ('pmean', []), ('hour', []), ('area', []) ])

    dic2 = OrderedDict([('lon_centre', []), ('lat_centre', []),
                        ('pmax', []), ('pmean', []), ('hour', []), ('area', []) ])

    dic3 = {'precip_pixel' : []}

    pick = [2,3,5,6,7,8]

    for v in res:

        for cnt, kk in enumerate(dic.keys()):

            dic[kk].append(v[cnt])  # omit kernel and kernelt

    for v in res:

        for cnt, kk in zip(pick,dic2.keys()):
            if (cnt == 0) or (cnt ==1) or (cnt ==4):
                continue
            dic2[kk].append(v[cnt])  # omit kernel and kernelt
        dic3['precip_pixel'].extend(v[4])


    pkl.dump(dic, open(out + 'trmm_blobs_1000km2.p', 'wb'))

    df = pd.DataFrame.from_dict(dic2)

    df.to_csv(out + 'trmm_cluster.csv')

    df = pd.DataFrame.from_dict(dic3)
    df.to_csv(out + 'trmm_pixel.csv')




def file_loop(f):

    ff = xr.open_dataset(f)
    print('Doing file: ' + f)

    larray = ff['p'].values.copy()
    precip = ff['p'].values
    lon = ff['lon'].values
    lat = ff['lat'].values

    hour = ff['time.hour'].values.tolist()

    larray[larray < 1] = 0
    larray[np.isnan(larray)] = 0

    labels, numL = label(larray)

    ret = []

    for l in np.unique(labels):
        if l == 0:
            continue

        blob = np.where(labels == l)

        if np.sum(len(blob[0]))< 7:  # at least 1000m2 for TRMM on 12km grid
            continue

        blon = lon[blob]
        blat = lat[blob]

        area = blon.size*144

        blon_c = np.min(blon) + (np.max(blon)-np.min(blon))
        blat_c = np.min(blat) + (np.max(blat)- np.min(blat))

        p = precip[blob]
        pmax = np.max(p)

        pmean = np.nanmean(p)

        ret.append((blon, blat, blon_c, blat_c, p ,pmax, pmean, hour, area))

    return ret
