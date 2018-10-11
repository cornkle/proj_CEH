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

def perSys():

    pool = multiprocessing.Pool(processes=5)
    tthresh = '-50'
    files = ua.locate(".nc", '/users/global/cornkle/data/CP4/CLOVER/MCS_-50_1000km2_JA_sahel')
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

    pkl.dump(mdic, open('/users/global/cornkle/data/CLOVER/saves/bulk_'+tthresh+'_zeroRain_gt1k_shear_CP4_JA_sahel.p',
                           'wb'))


def file_loop(f):
    #print('Doing file: ' + f)
    dic = xr.open_dataset(f)
    res = []
    out = dic['lw_out_PBLtop'].values
    outt = u_met.OLR_to_Tb(out)
    outt = outt-273.15
    outp = dic['lsRain'].values * 3600
    lon = dic['longitude'].values
    lat = dic['latitude'].values
    h = dic['time.hour'].values
    m = dic['time.month'].values
    y = dic['time.year'].values

    lons, lats = np.meshgrid(lon,lat)

    date = dic['time']

    t_thresh = -50  # -40C ~ 167 W m-2


    clat = np.min(lat)+((np.max(lat)-np.min(lat))*0.5)
    clon = np.min(lon) + ((np.max(lon) - np.min(lon)) * 0.5)
    pdb.set_trace()
    tt = np.min(outt[(np.isfinite(outp))&((outt<=t_thresh))])
    pp = np.max(outp[(np.isfinite(outp))&((outt<=t_thresh))])

    isfin = np.sum((np.isfinite(outp)) & ((outt<=t_thresh)))

    if isfin < 3:
        return

    try:
        pperc = outp[(outt<=t_thresh) & (outp>0.1)]
    except IndexError:
        pperc = np.nan
    tmean = np.mean(outt[(np.isfinite(outp)) & (outt<=t_thresh)])

    print(np.nanmax(outt))

    area = np.sum(outt<=t_thresh)


    ao40 = np.sum(outt<=t_thresh)
    po30 = np.sum(outp[(np.isfinite(outp))&((outt<=t_thresh))]>30)
    isfin = np.sum((np.isfinite(outp)) & ((outt<=t_thresh)))
    isnz = np.sum((outp>0.1) & ((outt<=t_thresh)))

    lon30 = lons[(np.isfinite(outp) & ((outt<=t_thresh)) & (outp > 30))]
    lat30 = lats[(np.isfinite(outp) & ((outt<=t_thresh)) & (outp > 30))]
    lonisfin = lons[(np.isfinite(outp)) & ((outt<=t_thresh))]
    latisfin = lats[(np.isfinite(outp)) & ((outt<=t_thresh))]

    latmin = lat.min()
    latmax = lat.max()

    p = outp[(np.isfinite(outp))&((outt<=t_thresh))]
    outt = outt[(np.isfinite(outp)) & ((outt<=t_thresh))]

    dic.close()
    return (tt,pp, area, ao40, tmean, pperc, clat, po30, isfin, outt, lon30, lat30, lonisfin, latisfin, h, m, latmin, latmax, isnz, clon, p, y, date)