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

    return {
        'hour' : [],
        'month' : [],
        'year' : [],

    }

def perSys():

    pool = multiprocessing.Pool(processes=5)
    tthresh = '-60'
    files = ua.locate(".nc", '/users/global/cornkle/data/CP4/CLOVER/MCS_-60_1000km2')
    print('Nb files', len(files))
    mdic = defaultdict(list)
    res = pool.map(file_loop, files)
    pool.close()
    #
    #res = [item for sublist in res for item in sublist]  # flatten list of lists

    #
    p=[]
    t=[]

    for v in res:
        try:
            mdic['tmin'].append(v[0])
            mdic['pmax'].append(v[1])
            mdic['area'].append(v[2])
            mdic['ao60'].append(v[3])
            mdic['tmean'].append(v[4])
            mdic['pperc'].extend(v[5])
            mdic['clat'].append(v[6])
            mdic['po30'].append(v[7])
            mdic['isfin'].append(v[8])
            mdic['t'].append(v[9])
            mdic['lon30'].extend(v[10])
            mdic['lat30'].extend(v[11])
            mdic['lonisfin'].extend(v[12])
            mdic['latisfin'].extend(v[13])
            mdic['hour'].append(v[14])
            mdic['month'].append(v[15])
            mdic['latmin'].append(v[16])
            mdic['latmax'].append(v[17])
            mdic['isnz'].append(v[18])
            mdic['clon'].append(v[19])
            mdic['p'].append(v[20])
            mdic['year'].append(v[21])
            mdic['date'].append(v[22])
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

    pkl.dump(mdic, open('/users/global/cornkle/data/CLOVER/saves/bulk_'+tthresh+'_zeroRain_gt1k_shear_CP4.p',
                           'wb'))


def file_loop(f):
    #print('Doing file: ' + f)
    dic = xr.open_dataset(f)
    res = []
    outt = dic['lw_out_PBLtop'].values
    outp = dic['lsRain'].values
    outu_srfc = dic['u_srfc'].values
    outu_mid = dic['u_mid'].values
    outshear = dic['shear'].values
    outq = dic['q_pl'].values

    lon = dic['longitude'].values
    lat = dic['latitude'].values
    h = dic['time.hour'].values
    m = dic['time.month'].values
    y = dic['time.year'].values

    lons, lats = np.meshgrid(lon,lat)

    date = dic['time']

    t_thresh = -60  # -40C ~ 167 W m-2
    mask = np.isfinite(outp) & (outt<=t_thresh) & np.isfinite(outq)

    if np.sum(mask) < 3:
        return

    area = np.sum(mask)*(4.4**2)

    clat = np.min(lat)+((np.max(lat)-np.min(lat))*0.5)
    clon = np.min(lon) + ((np.max(lon) - np.min(lon)) * 0.5)

    tmin = np.min(outt[mask])
    tmean = np.mean(outt[mask])
    pmax = np.max(outp[mask])
    pmean = np.mean(outp[mask])
    qmax = np.max(outq[mask])
    qmean = np.mean(outq[mask])
    umax_srfc = np.max(outu_srfc[mask])
    umean_srfc = np.mean(outu_srfc[mask])
    umin_mid = np.min(outu_mid[mask])
    umean_mid = np.mean(outu_mid[mask])
    shear_min = np.min(outshear[mask])
    shear_mean = np.mean(outshear[mask])


    ao60 = np.sum(outt[mask])
    po30 = np.sum(outp[mask]>30)
    isfin = np.sum(mask)
    isnz = np.sum(outp[mask]>0.1)

    latmin = lat.min()
    latmax = lat.max()

    p = outp[mask]
    t = outt[mask]
    q = outq[mask]
    u_mid = outu_mid[mask]
    u_srfc = outu_srfc[mask]
    shear = outshear[mask]

    dic.close()
    return (h, m, y, lons, lats, date, area, clat, clon, tmin)