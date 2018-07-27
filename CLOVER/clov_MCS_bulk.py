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
import cartopy
import pdb

def perSys():

    pool = multiprocessing.Pool(processes=5)
    tthresh = '-40'
    files = ua.locate(".nc", '/users/global/cornkle/MCSfiles/WA350_4-8N_14W-10E_'+tthresh+'/')
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
            mdic['pc'].append(v[21])
            mdic['year'].append(v[22])
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

    pkl.dump(mdic, open('/users/global/cornkle/data/CLOVER/saves/bulk_'+tthresh+'_zeroRain_gt15k_shear.p',
                           'wb'))


def file_loop(f):
    #print('Doing file: ' + f)
    dic = xr.open_dataset(f)
    res = []
    outt = dic['tc_lag0'].values
    outp = dic['p'].values
    outpc = dic['pconv'].values
    lon = dic['lon'].values
    lat = dic['lat'].values
    h = dic['time.hour'].values
    m = dic['time.month'].values
    y = dic['time.year'].values

    date = dic['time']

    t_thresh = -40
    pos = (outt<=t_thresh)

    clat = np.min(lat)+((np.max(lat)-np.min(lat))*0.5)
    clon = np.min(lon) + ((np.max(lon) - np.min(lon)) * 0.5)


    era = xr.open_dataset(constants.)


    tt = np.min(outt[(np.isfinite(outp))&((outt<=t_thresh))])
    pp = np.max(outp[(np.isfinite(outp))&((outt<=t_thresh))])
    ppmin = np.min(outp[(np.isfinite(outp)) & ((outt<=t_thresh))])

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

    if  (pp>200) | (ppmin<0) | (area*25 < 15000) : #(pp<0.1)  | (area*25>800000)
        return

    ao40 = np.sum(outt<=t_thresh)
    po30 = np.sum(outp[(np.isfinite(outp))&((outt<=t_thresh))]>30)
    isfin = np.sum((np.isfinite(outp)) & ((outt<=t_thresh)))
    isnz = np.sum((outp>0.1) & ((outt<=t_thresh)))

    lon30 = lon[(np.isfinite(outp) & ((outt<=t_thresh)) & (outp > 30))]
    lat30 = lat[(np.isfinite(outp) & ((outt<=t_thresh)) & (outp > 30))]
    lonisfin = lon[(np.isfinite(outp)) & ((outt<=t_thresh))]
    latisfin = lat[(np.isfinite(outp)) & ((outt<=t_thresh))]

    latmin = lat.min()
    latmax = lat.max()

    p = outp[(np.isfinite(outp))&((outt<=t_thresh))]
    pc = outpc[(np.isfinite(outpc)) & ((outt<=t_thresh))]
    outt = outt[(np.isfinite(outp)) & ((outt<=t_thresh))]

    dic.close()
    return (tt,pp, area, ao40, tmean, pperc, clat, po30, isfin, outt, lon30, lat30, lonisfin, latisfin, h, m, latmin, latmax, isnz, clon, p, pc, y)