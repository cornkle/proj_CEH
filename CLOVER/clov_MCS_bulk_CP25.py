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
from metpy import calc
from metpy.units import units

def dictionary():

    dic = {}
    vars = ['hour', 'month', 'year', 'area',
            'lon', 'lat', 'clon', 'clat',
            'tmin', 'tmean', 'thetamean', 'thetamax', 'tmidmax', 'tmidmean', 'tsrfcmax', 'tsrfcmean',
            'pmax', 'pmean',
            'qmax' , 'qmean',
            'tcwv', 'tcwvmean',
            'tgrad', 'tdiff',
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

    tthresh = '-50'
    files = ua.locate(".nc", '/media/ck/Elements/Africa/WestAfrica/CP4/CP25_16-19UTC_future_5000km2_-50C_TCWV')  #CP25_-50C_5000km2
    print('Nb files', len(files))
    for y in range(1999,2007):

        yfiles = []
        for f in files:
            if str(y) in f:
                yfiles.append(f)
        pool = multiprocessing.Pool(processes=4)

        mdic = dictionary() #defaultdict(list)
        print('Yearly files', len(yfiles))
        # ipdb.set_trace()
        # res = pool.map(file_loop, yfiles)
        # pool.close()

        res=[]
        for f in yfiles:
            res.append(file_loop(f))

        ipdb.set_trace()
        keys = mdic.keys()
        for v in res:
            for k in keys:
                try:
                    mdic[k].append(v[k])
                except TypeError:
                    continue


        pkl.dump(mdic, open(cnst.network_data +'data/CLOVER/saves/bulk_'+tthresh+'_5000km2_P25means_hourly_SAHEL_15kmprecip_WA_5-20N_-50C_TCWV_fut_'+str(y)+'.p',
                           'wb'))
        print('Saved file')


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
    tmid = dic['t_mid'].values
    tsrfc = dic['t_srfc'].values

    pes = units.Quantity(650, 'hPa')
    tes = units.Quantity(tmid+273.15, 'K')

    theta_low = u_met.theta_e(925, tsrfc, outq/1000.)

    theta_es = calc.saturation_equivalent_potential_temperature(pes, tes)

    theta = theta_low-(np.array(theta_es)-273.15)
    tcwv = dic['tcwv'].values
    tgrad = dic.attrs['Tgrad']
    tdiff = dic.attrs['Tgradbox']

    #print(theta)

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

    for var in [outt,outp,outu_srfc,outu_mid,outshear,outq,theta,tmid,tsrfc]:
        var[antimask] = np.nan


    if np.sum(mask) < 3:
        print('NOT ENOUGH IN MASK')
        return

    maxpos = np.unravel_index(np.nanargmax(outp), outp.shape)
    minpos = np.unravel_index(np.nanargmin(outt), outt.shape)

    out['area'] = np.sum((outt<=t_thresh))

    out['clat'] = np.min(out['lat'])+((np.max(out['lat'])-np.min(out['lat']))*0.5)
    out['clon'] = np.min(out['lon']) + ((np.max(out['lon']) - np.min(out['lon'])) * 0.5)

    out['tmin'] = np.nanmin(outt)
    out['tmean'] = np.mean(outt[mask])
    out['tgrad'] = tgrad
    out['tdiff'] = tdiff

    out['pmax'] = np.nanmean(ua.cut_kernel(outp,maxpos[1], maxpos[0],1)) # degrade rainfall to 30km

    #out['pmax'] = np.max(outp[mask]) # rain at 4.4km
    out['pmean'] = np.mean(outp[mask])
    out['qmax'] = np.nanmean(ua.cut_kernel(outq,minpos[1], minpos[0],3)) #np.max(outq[mask]) #30km q and shear

    out['qmean'] = np.nanmean(outq[mask])
    out['umax_srfc'] = np.max(outu_srfc[mask])
    out['umean_srfc'] = np.nanmean(outu_srfc[mask])
    out['umin_mid'] = np.min(outu_mid[mask])
    out['umean_mid'] = np.mean(outu_mid[mask])

    out['shearmin'] =  np.nanmean(ua.cut_kernel(outshear,minpos[1], minpos[0],3))

    out['shearmean'] = np.mean(outu_mid[mask]) - np.nanmean(outu_srfc[mask]) #np.mean(outshear[mask])
    #out['thetamax'] = np.max(theta[mask])

    out['thetamax'] = np.nanmean(ua.cut_kernel(theta,minpos[1], minpos[0],3))
    out['tcwv'] = np.nanmean(ua.cut_kernel(tcwv, minpos[1], minpos[0], 3))
    out['tcwvmean'] = np.mean(tcwv[mask])

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

    print(out['pmax'])

    return out
