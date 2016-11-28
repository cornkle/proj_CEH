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
import cartopy

def perSys():

    pool = multiprocessing.Pool(processes=5)
    files = ua.locate(".nc", '/users/global/cornkle/MCSfiles/WA15_big_-40_15W-20E/')
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

    pkl.dump(mdic, open('/users/global/cornkle/C_paper/wavelet/saves/bulk_40big.p',
                           'wb'))




def file_loop(f):
    print('Doing file: ' + f)
    dic = xr.open_dataset(f)
    res = []
    outt = dic['tc_lag0'].values
    outp = dic['p'].values
    lon = dic['lon'].values
    lat = dic['lat'].values
    h = dic['time.hour'].values
    m = dic['time.month'].values
    clat = np.min(dic.lat)+((np.max(dic.lat)-np.min(dic.lat))*0.5)
    clon = np.min(dic.lon) + ((np.max(dic.lon) - np.min(dic.lon)) * 0.5)

    tt = np.min(outt[(np.isfinite(outp))&(np.isfinite(outt))])
    pp = np.max(outp[(np.isfinite(outp))&(np.isfinite(outt))])
    try:
        pperc = outp[(outt<-40) & (outp>0.1)]
    except IndexError:
        pperc = np.nan
    tmean = np.mean(outt[(np.isfinite(outp)) & (np.isfinite(outt))])
    area = np.sum(outt<=-40)

    if (area*25 < 300) or (pp<0.1) or (pp>200):
        return

    ao40 = np.sum(outt<=-40)
    po30 = np.sum(outp[(np.isfinite(outp))&(np.isfinite(outt))]>30)
    isfin = np.size(outp[(np.isfinite(outp)) & (np.isfinite(outt))])
    isnz = np.size(outp[(outp>0) & (np.isfinite(outt))])

    lon30 = lon[(np.isfinite(outp) & (np.isfinite(outt)) & (outp > 30))]
    lat30 = lat[(np.isfinite(outp) & (np.isfinite(outt)) & (outp > 30))]
    lonisfin = lon[(np.isfinite(outp)) & (np.isfinite(outt))]
    latisfin = lat[(np.isfinite(outp)) & (np.isfinite(outt))]

    latmin = lat.min()
    latmax = lat.max()



    dic.close()
    return (tt,pp, area, ao40, tmean, pperc, clat, po30, isfin, outt, lon30, lat30, lonisfin, latisfin, h, m, latmin, latmax, isnz)


def file_loop_wavelet(f):
    print('Doing file: ' + f)
    dic = xr.open_dataset(f)
    res = []
    outt = dic['tc_lag0'].values
    outp = dic['p'].values
    lon = dic['lon'].values
    lat = dic['lat'].values
    h = dic['time.hour'].values
    m = dic['time.month'].values
    clat = np.min(dic.lat)+((np.max(dic.lat)-np.min(dic.lat))*0.5)

    tt = np.min(outt[(np.isfinite(outp))&(np.isfinite(outt))])
    pp = np.max(outp[(np.isfinite(outp))&(np.isfinite(outt))])
    try:
        pperc = outp[(outt<-40) & (outp>0.1)]
    except IndexError:
        pperc = np.nan
    tmean = np.mean(outt[(np.isfinite(outp)) & (np.isfinite(outt))])
    area = np.sum(outt<=-40)

    if (area*25 < 15000) or (pp<0.1) or (pp>200):
        return

    ao40 = np.sum(outt<=-40)
    po30 = np.sum(outp[(np.isfinite(outp))&(np.isfinite(outt))]>30)
    isfin = np.size(outp[(np.isfinite(outp)) & (np.isfinite(outt))])

    lon30 = lon[(np.isfinite(outp) & (np.isfinite(outt)) & (outp > 30))]
    lat30 = lat[(np.isfinite(outp) & (np.isfinite(outt)) & (outp > 30))]
    lonisfin = lon[(np.isfinite(outp)) & (np.isfinite(outt))]
    latisfin = lat[(np.isfinite(outp)) & (np.isfinite(outt))]

    latmin = lat.min()
    latmax = lat.max()



    dic.close()
    return (tt,pp, area, ao40, tmean, pperc, clat, po30, isfin, outt, lon30, lat30, lonisfin, latisfin, h, m, latmin, latmax)