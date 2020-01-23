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
from utils import constants as cnst
import metpy
from metpy import calc
from metpy.units import units


def dictionary():

    dic = {}
    vars = ['hour', 'month', 'year', 'area',
            'lon', 'lat', 'clon', 'clat',
            'tmin', 'tmean', 'thetamean', 'thetamax', 'tmidmax', 'tmidmean', 'tsrfcmax', 'tsrfcmean',
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

    pool = multiprocessing.Pool(processes=4)
    tthresh = '-50'
    files = ua.locate(".nc", cnst.network_data +'data/CP4/CLOVER/CP4_18UTC_5000km2_-50_5-20N_new')  #CP25_-50C_5000km2
    print('Nb files', len(files))
    mdic = dictionary() #defaultdict(list)
    res = pool.map(file_loop, files)
    pool.close()
    # res=[]
    # for f in files:
    #     res.append(file_loop(f))

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

    pkl.dump(mdic, open(cnst.network_data +'data/CLOVER/saves/bulk_'+tthresh+'_5000km2_CP4_ERA5_30km_WA_5-20N_new.p',
                           'wb'))


def file_loop(f):
    print('Doing file: ' + f)
    dic = xr.open_dataset(f)
    out = dictionary()

    outt = dic['lw_out_PBLtop'].values
    try:
        outp = dic['lsRain'].values
    except KeyError:
        outp = dic['totRain'].values
    outu_srfc = dic['u_srfc'].values
    outu_mid = dic['u_mid'].values
    outshear = dic['shear'].values
    outq = dic['q_srfc'].values
    #theta = dic['theta'].values
    tmid = dic['t_mid'].values
    tsrfc = dic['t_srfc'].values

    pes = units.Quantity(650, 'hPa')
    tes = units.Quantity(tmid+273.15, 'K')

    theta_low = u_met.theta_e(925, tsrfc, outq/1000.)

    theta_es = calc.saturation_equivalent_potential_temperature(pes, tes)

    theta = theta_low-(np.array(theta_es)-273.15)

    print(theta)

    #ipdb.set_trace()

    out['lon'] = dic['longitude'].values
    out['lat'] = dic['latitude'].values
    out['hour'] = dic['time.hour'].item()
    out['month'] = dic['time.month'].item()
    out['year'] = dic['time.year'].item()
    out['date'] = dic['time'].values

    t_thresh = -50  # -40C ~ 167 W m-2
    mask = np.isfinite(outp) & (outt<=t_thresh) & np.isfinite(outq) & np.isfinite(outshear)
    antimask = ~mask

    for var in [outt,outp,outu_srfc,outu_mid,outshear,outq,theta,tmid]:
        var[antimask] = np.nan


    if np.sum(mask) < 3:
        return

    maxpos = np.unravel_index(np.nanargmax(outp), outp.shape)
    minpos = np.unravel_index(np.nanargmin(outt), outt.shape)
    minshear = np.unravel_index(np.nanargmin(outshear), outt.shape)
    maxq = np.unravel_index(np.nanargmax(outq), outt.shape)

    out['area'] = np.sum((outt<=t_thresh))

    out['clat'] = np.min(out['lat'])+((np.max(out['lat'])-np.min(out['lat']))*0.5)
    out['clon'] = np.min(out['lon']) + ((np.max(out['lon']) - np.min(out['lon'])) * 0.5)

    out['tmin'] = np.nanmin(outt)
    out['tmean'] = np.mean(outt[mask])


    out['pmax'] = np.nanmean(ua.cut_kernel(outp,maxpos[1], maxpos[0],2)) # degrade rainfall to 13km

    #out['pmax'] = np.max(outp[mask]) # degrade rainfall to 20km
    out['pmean'] = np.mean(outp[mask])
    out['qmax'] = np.nanmean(ua.cut_kernel(outq,maxq[1], maxq[0],2)) #np.max(outq[mask])

    out['qmean'] = np.nanmean(outq[mask])
    out['umax_srfc'] = np.max(outu_srfc[mask])
    out['umean_srfc'] = np.nanmean(outu_srfc[mask])
    out['umin_mid'] = np.min(outu_mid[mask])
    out['umean_mid'] = np.mean(outu_mid[mask])


    out['shearmin'] =  np.nanmean(ua.cut_kernel(outshear,minshear[1], minshear[0],2))

    out['shearmean'] = np.mean(outshear[mask])
    out['thetamax'] = np.max(theta[mask])
    out['thetamean'] = np.mean(theta[mask])
    out['tmidmax'] = np.max(tmid[mask])
    out['tmidmean'] = np.mean(tmid[mask])
    out['tsrfcmax'] = np.max(tsrfc[mask])
    out['tsrfcmean'] = np.mean(tsrfc[mask])

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
