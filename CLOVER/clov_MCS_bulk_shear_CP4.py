import numpy as np
import xarray as xr
from utils import u_arrays as ua
from scipy import ndimage
import matplotlib.pyplot as plt
import multiprocessing
import ipdb
import pickle as pkl
from collections import defaultdict
import cartopy.crs as ccrs
from utils import constants
from utils import u_met
import cartopy
import pdb

def dictionary():

    dic = {}
    vars = ['hour', 'month', 'year', 'area',
            'lon', 'lat', 'clon', 'clat',
            'tmin', 'tmean',
            'pmax', 'pmean',
            'qmax' , 'qmean',
            'umax_srfc', 'umean_srfc',
            'umin_mid', 'umean_mid',
            'shearmin', 'shearmean',
            'pgt30', 'pgt01'
            'isvalid',
             't', 'p', 'q', 'u_srfc', 'u_mid', 'shear' ]

    for v in vars:
        dic[v] = []
    return dic

def perSys():

    pool = multiprocessing.Pool(processes=5)
    tthresh = '-40'
    files = ua.locate(".nc", '/users/global/cornkle/data/CP4/CLOVER/MCS_-40_5000km2_JJAS')
    print('Nb files', len(files))
    mdic = dictionary() #defaultdict(list)
    res = pool.map(file_loop, files)
    pool.close()
    #
    #res = [item for sublist in res for item in sublist]  # flatten list of lists

    keys = mdic.keys()
    for v in res:
        for k in keys:
            try:
                mdic[k].append(v[k])
            except TypeError:
                continue


        # if v[2]*25 > 1000000:
        #     tplt = v[9]
        #     tplt[np.where(tplt==np.nan)]=0
            # f = plt.figure()
            # ax = plt.axes(projection=ccrs.PlateCarree())
            # plt.contourf(v[10], v[11], tplt, transform=ccrs.PlateCarree())
            # ax.coastlines()
            # plt.colorbar()
            # ax.add_feature(cartopy.feature.BORDERS, linestyle='--')


    # f = plt.figure()
    # siz = 3
    #
    # ax = f.add_subplot(1, 1, 1)
    # plt.scatter(mdic['tmin'], mdic['pmax'])
    # plt.title('bulk', fontsize=9)

    pkl.dump(mdic, open('/users/global/cornkle/data/CLOVER/saves/bulk_'+tthresh+'_zeroRain_gt1k_shear_CP4_JJASNORTH.p',
                           'wb'))


def file_loop(f):
    print('Doing file: ' + f)
    dic = xr.open_dataset(f)
    out = dictionary()

    outt = dic['lw_out_PBLtop'].values
    outp = dic['lsRain'].values
    outu_srfc = dic['u_srfc'].values
    outu_mid = dic['u_mid'].values
    outshear = dic['shear'].values
    outq = dic['q_pl'].values

    out['lon'] = dic['longitude'].values
    out['lat'] = dic['latitude'].values
    out['hour'] = dic['time.hour'].item()
    out['month'] = dic['time.month'].item()
    out['year'] = dic['time.year'].item()
    out['date'] = dic['time'].values

    t_thresh = -40  # -40C ~ 167 W m-2
    mask = np.isfinite(outp) & (outt<=t_thresh) & np.isfinite(outq) & np.isfinite(outshear)

    if np.sum(mask) < 3:
        return

    out['area'] = np.sum(mask)*(4.4**2)

    out['clat'] = np.min(out['lat'])+((np.max(out['lat'])-np.min(out['lat']))*0.5)
    out['clon'] = np.min(out['lon']) + ((np.max(out['lon']) - np.min(out['lon'])) * 0.5)

    out['tmin'] = np.min(outt[mask])
    out['tmean'] = np.mean(outt[mask])
    out['pmax'] = np.max(outp[mask])
    out['pmean'] = np.mean(outp[mask])
    out['qmax'] = np.max(outq[mask])
    out['qmean'] = np.mean(outq[mask])
    out['umax_srfc'] = np.max(outu_srfc[mask])
    out['umean_srfc'] = np.mean(outu_srfc[mask])
    out['umin_mid'] = np.min(outu_mid[mask])
    out['umean_mid'] = np.mean(outu_mid[mask])
    out['shearmin'] = np.min(outshear[mask])
    out['shearmean'] = np.mean(outshear[mask])

    out['pgt30'] = np.sum(outp[mask]>30)
    out['isvalid'] = np.sum(mask)
    out['pgt01'] = np.sum(outp[mask]>0.1)

    out['p'] = outp[mask]
    out['t'] = outt[mask]
    out['q'] = outq[mask]
    out['u_mid'] = outu_mid[mask]
    out['u_srfc'] = outu_srfc[mask]
    out['shear'] = outshear[mask]

    dic.close()

    return out