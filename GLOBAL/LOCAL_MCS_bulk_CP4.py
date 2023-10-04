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
import glob


def dictionary():

    dic = {}
    vars = ['hour', 'month', 'year', 'area',
            'lon', 'lat', 'clon', 'clat',
            'tmin', 'tmean', 'thetamean', 'thetamax', 'tmidmax', 'tmidmean', 'tsrfcmax', 'tsrfcmean',
            'pmax', 'pmean',
            'qmax' , 'qmean',
            'tcwv', 'tcwvmean',
            'tgrad',
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

    pool = multiprocessing.Pool(processes=3)
    tthresh = '-50'
    files = glob.glob('/media/ck/LStorage/global_water/CP_models/MCS_files/MODELS/CP4_box/CP4_allHours_historical_5000km2_-50_WAf_box_v2/*.nc') # CP4_18UTC_5000km2_-50_5-20N_new')  #CP25_-50C_5000km2
    #files = glob.glob('/media/ck/LStorage/global_water/CP_models/MCS_files/MODELS/CP4_box/CP4_allHours_future_5000km2_-50_WAf_box_v2/*.nc')
    print('Nb files', len(files))
    for y in range(1998,2007):

        yfiles = []
        for f in files:
            if str(y) in f:
                yfiles.append(f)
        pool = multiprocessing.Pool(processes=4)

        mdic = dictionary() #defaultdict(list)
        print('Yearly files', len(yfiles))
        # res = pool.map(file_loop, yfiles)
        # pool.close()

        res=[]
        for f in files:
            res.append(file_loop(f))
            ipdb.set_trace()

        keys = mdic.keys()
        for v in res:
            for k in keys:
                try:
                    mdic[k].append(v[k])
                except TypeError:
                    continue


        #ipdb.set_trace()
        pkl.dump(mdic, open(cnst.network_data +'data/LMCS/CP4_study_saves/bulk_CP4/bulk_'+tthresh+'_5000km2_CP4means_hourly_SAHEL_15kmprecip_WA_18W-25E_9-25N_-50C_LMCSfiles_17h_hist_'+str(y)+'.p',
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
    tmid = dic['t_mid'].values
    tsrfc = dic['t_srfc'].values

    pes = units.Quantity(650, 'hPa')
    pes_down = units.Quantity(925, 'hPa')
    tes_up = units.Quantity(tmid+273.15, 'K')
    qes = units.Quantity(outq/1000, 'kg/kg')
    tes_down = units.Quantity(tsrfc+273.15, 'K')

    theta_low = u_met.theta_e(925, tsrfc, outq/1000.)
    rh = calc.relative_humidity_from_specific_humidity(pes_down, tes_down, qes)
    dew = calc.dewpoint_from_specific_humidity(pes_down, tes_down, qes)
    theta_e = calc.equivalent_potential_temperature(pes_down, tes_down,dew)
    theta_es = calc.saturation_equivalent_potential_temperature(pes, tes_up)
    ipdb.set_trace()
    theta = theta_low-(np.array(theta_es)-273.15)
    tcwv = dic['tcwv'].values


    out['lon'] = dic.minlon
    out['lat'] = dic.minlat
    out['hour'] = dic['time.hour'].item()
    out['month'] = dic['time.month'].item()
    out['year'] = dic['time.year'].item()
    out['date'] = dic['time'].values

    out['rh'] = rh

    t_thresh = -50  # -40C ~ 167 W m-2
    mask = np.isfinite(outp) & (outt<=t_thresh) & np.isfinite(outq) & np.isfinite(outshear)
    antimask = ~mask

    for var in [outt,outp,outu_srfc,outu_mid,outshear,outq,theta,tmid,tsrfc]:
        var[antimask] = np.nan


    if np.sum(mask) < 3:
        return

    maxpos = np.unravel_index(np.nanargmax(outp), outp.shape)
    minpos = np.unravel_index(np.nanargmin(outt), outt.shape)

    tgrad_lat = dic.isel(longitude=slice(minpos[1] - 3, minpos[1] + 3)).mean('longitude').squeeze()
    ipdb.set_trace()
    tgrad = tgrad_lat.polyfit(dim='latitude', deg=1)

    out['area'] = np.sum((outt<=t_thresh))


    out['tmin'] = np.nanmin(outt)
    out['tmean'] = np.mean(outt[mask])
    out['tgrad'] = tgrad

    out['pmax'] = np.nanmean(ua.cut_kernel(outp,maxpos[1], maxpos[0],1)) # degrade rainfall to 15km

    out['pmax_native'] = np.max(outp[mask]) # rain at 4.4km
    out['pmean'] = np.mean(outp[mask])
    out['qmax30'] = np.nanmean(ua.cut_kernel(outq,minpos[1], minpos[0],3)) #np.max(outq[mask]) #30km q and shear

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
    out['tcwvmean'] = np.nanmean(tcwv[mask])
    out['rh'] = np.nanmean(rh[mask])

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
