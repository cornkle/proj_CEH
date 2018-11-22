import numpy as np
import xarray as xr
from utils import u_arrays as ua
import matplotlib.pyplot as plt
import multiprocessing
import pickle as pkl
from collections import defaultdict
from utils import constants as cnst
import pdb


def dictionary():

    dic = {}
    vars = ['hour', 'month', 'year', 'area',
            'lon', 'lat', 'clon', 'clat',
            'tmin', 'tmean',
            'pmax', 'pmean',
            'q925', 'q650',
            'u925', 'u650',
            'v925', 'v650',
            'w925', 'w650',
            'rh925','rh650',
            't925', 't650',
            'div925','div650',
            'pv925', 'pv650',
            'shear',
            'pgt30', 'pgt01'
            'isvalid',
             't', 'p' ]

    for v in vars:
        dic[v] = []
    return dic

def perSys():

    pool = multiprocessing.Pool(processes=4)
    tthresh = '-40'
    files = ua.locate(".nc", cnst.network_data + 'MCSfiles/WA5000_4-8N_13W-13E_'+tthresh+'_18UTC/')

    print('Nb files', len(files))
    mdic = dictionary() #defaultdict(list)
    res = pool.map(file_loop, files)
    pool.close()
    #
    #
    # for f in files:
    #     file_loop(f)

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


    pkl.dump(mdic, open(cnst.network_data + 'data/CLOVER/saves/bulk_'+tthresh+'_zeroRain_gt5k_-40thresh_OBSera_thicklayer.p',
                           'wb'))


def file_loop(f):
    print('Doing file: ' + f)
    dic = xr.open_dataset(f)
    era = xr.open_dataset(cnst.ERA_DAILY_PL12UTC)

    getera =np.where((era['time.day']==dic['time.day']) & (era['time.month']==dic['time.month']) & (era['time.year']==dic['time.year']))
    try:
        era_day = era.isel(time=int(getera[0]))
    except TypeError:
        print('Era missing')
        return

    out = dictionary()
    res = []
    outt = dic['tc_lag0'].values
    outp = dic['p'].values
    outpc = dic['pconv'].values

    tminpos = np.where(dic['tc_lag0'].values == np.nanmin(dic['tc_lag0'].values)) # era position close to min temp
    if len(tminpos[0])>1:
        ptmax = np.nanmax((dic['p'].values)[tminpos])
        if ptmax > 0:
            prpos = np.where((dic['p'].values)[tminpos] == ptmax)
            tminpos = ((tminpos[0])[prpos], (tminpos[1])[prpos] )
        else:
            tminpos = ((tminpos[0])[0], (tminpos[1])[0])

    elon = dic['lon'].values[tminpos]
    elat = dic['lat'].values[tminpos]

    #e925 = era_day.sel(latitude=elat, longitude=elon, level=925, method='nearest')
    e925 = era_day.sel(level=slice(925,900)).mean('level').sel(latitude=elat, longitude=elon , method='nearest')
    #e650 = era_day.sel(latitude=elat, longitude=elon, level=650, method='nearest')
    e650 = era_day.sel(level=slice(600,601)).mean('level').sel(latitude=elat, longitude=elon , method='nearest')


    out['lon'] = dic['lon'].values
    out['lat'] = dic['lat'].values
    out['hour'] = dic['time.hour'].item()
    out['month'] = dic['time.month'].item()
    out['year'] = dic['time.year'].item()
    out['date'] = dic['time'].values

    t_thresh = -40  # -40C ~ 167 W m-2
    mask = np.isfinite(outp) & (outt<=t_thresh) & np.isfinite(outt)

    if np.sum(mask) < 3:
        return

    out['clat'] = np.min(out['lat'])+((np.max(out['lat'])-np.min(out['lat']))*0.5)
    out['clon'] = np.min(out['lon']) + ((np.max(out['lon']) - np.min(out['lon'])) * 0.5)

    isfin = np.sum((np.isfinite(outp)) & ((outt<=t_thresh)))

    if isfin < 3:
        return

    print(np.nanmax(outt[mask]))   # can be bigger than cutout threshold because of interpolation to 5km grid after cutout

    out['area'] = np.sum(mask)*(4.4**2)

    out['clat'] = np.min(out['lat'])+((np.max(out['lat'])-np.min(out['lat']))*0.5)
    out['clon'] = np.min(out['lon']) + ((np.max(out['lon']) - np.min(out['lon'])) * 0.5)

    out['tmin'] = np.min(outt[mask])
    out['tmean'] = np.mean(outt[mask])
    out['pmax'] = np.max(outp[mask])
    out['pmean'] = np.mean(outp[mask])
    try:
        out['q925'] =float(e925['q'])
    except TypeError:
        pdb.set_trace()
    out['q650'] = float(e650['q'])
    out['v925'] = float(e925['v'])
    out['v650'] = float(e925['v'])
    out['u925'] = float(e925['u'])
    out['u650'] = float(e650['u'])
    out['w925'] = float(e925['w'])
    out['w650'] = float(e650['w'])
    out['rh925'] = float(e925['r'])
    out['rh650'] = float(e650['r'])
    out['t925'] = float(e925['t'])
    out['t650'] = float(e650['t'])
    out['pv925'] = float(e925['pv'])
    out['pv650'] = float(e650['pv'])
    out['div925'] = float(e925['d'])
    out['div650'] = float(e650['d'])


    out['shear'] = float(e650['u']-e925['u'])

    out['pgt30'] = np.sum(outp[mask]>30)
    out['isvalid'] = np.sum(mask)
    out['pgt01'] = np.sum(outp[mask]>0.1)
    #
    out['p'] = outp[mask]
    out['t'] = outt[mask]

    dic.close()

    return out
