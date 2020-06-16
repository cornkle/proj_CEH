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
            'tmin', 'tmean', 'thetamax', 'tmidmax',  'tsrfcmax',
            'pmax', 'pmean',
            'qmax' ,
            'tcwv',
            'tgrad', 'tdiff',
            'umax_srfc', 'umean_srfc',
            'umin_mid',
            'shearmin',
            'pgt30', 'pgt01'
            'isvalid',
             't', 'p' ]

    for v in vars:
        dic[v] = []
    return dic

def perSys():


    tthresh = '-50'
    files = ua.locate(".nc", '/media/ck/Elements/Africa/WestAfrica/CP4/CP4_16-19UTC_future_5000km2_-50C_TCWV') # CP4_18UTC_5000km2_-50_5-20N_new')  #CP25_-50C_5000km2
    print('Nb files', len(files))

    for y in range(1998,2007):

        yfiles = []
        for f in files:
            if str(y) in f:
                yfiles.append(f)
        pool = multiprocessing.Pool(processes=4)

        mdic = dictionary() #defaultdict(list)
        print('Yearly files', len(yfiles))
        res = pool.map(file_loop, yfiles)
        pool.close()

        # res=[]
        # for f in files:
        #     res.append(file_loop(f))


        keys = mdic.keys()
        for v in res:
            for k in keys:
                try:
                    mdic[k].append(v[k])
                except TypeError:
                    continue


        #ipdb.set_trace()
        pkl.dump(mdic, open(cnst.network_data +'data/CLOVER/saves/bulk_'+tthresh+'_5000km2_CP40.7deg_monthly_SAHEL_15kmprecip_WA_5-20N_-50C_TCWV_fut_'+str(y)+'.p',
                           'wb'))
        print('Written file for ', y)


def file_loop(f):
    print('Doing file: ' + f)
    dic = xr.open_dataset(f)

    tag = 'fut'

    mfiles = '/home/ck/DIR/cornkle/data/CP4/CLOVER/CP4_monthly_'+tag+'/'

    y = dic['time.year'].values
    m = dic['time.month'].values

    out = dictionary()

    outt = dic['lw_out_PBLtop'].values
    try:
        outp = dic['lsRain'].values
    except KeyError:
        outp = dic['totRain'].values


    out['lon'] = dic['longitude'].values
    out['lat'] = dic['latitude'].values
    out['hour'] = dic['time.hour'].item()
    out['month'] = dic['time.month'].item()
    out['year'] = dic['time.year'].item()
    out['date'] = dic['time'].values

    t_thresh = -50  # -40C ~ 167 W m-2

    maxpos = np.unravel_index(np.nanargmax(outp), outp.shape)
    minpos = np.unravel_index(np.nanargmin(outt), outt.shape)
    #ipdb.set_trace()
    latpick = out['lat'][minpos[0]]
    lonpick = out['lon'][minpos[1]]

    out['area'] = np.sum((outt<=t_thresh))

    out['clat'] = np.min(out['lat'])+((np.max(out['lat'])-np.min(out['lat']))*0.5)
    out['clon'] = np.min(out['lon']) + ((np.max(out['lon']) - np.min(out['lon'])) * 0.5)

    if (out['clat']<9) | (out['clon']<-15) | (out['clon']>15):
        print('MCS out of box')
        return

    tcwv = xr.open_dataarray(mfiles+'tcwv/tcwv_monthly_12UTC_0.7deg_'+tag+'_4km_'+str(y)+str(m).zfill(2)+'.nc')
    tpl = xr.open_dataarray(mfiles + 't_pl/t_pl_monthly_12UTC_0.7deg_'+tag+'_4km_' + str(y) + str(m).zfill(2) + '.nc')
    uu = xr.open_dataarray(mfiles + 'u_pl/u_pl_monthly_12UTC_0.7deg_'+tag+'_4km_' + str(y) + str(m).zfill(2) + '.nc')
    qq = xr.open_dataarray(mfiles + 'q_pl/q_pl_monthly_12UTC_0.7deg_'+tag+'_4km_' + str(y) + str(m).zfill(2) + '.nc')

    tcwvout = tcwv.sel(latitude=latpick, longitude=lonpick, method='nearest', tolerance=0.7)
    tt = tpl.sel(latitude=latpick, longitude=lonpick, method='nearest', tolerance=0.7)
    uu = uu.sel(latitude=latpick, longitude=lonpick, method='nearest', tolerance=0.7)
    qq = qq.sel(latitude=latpick, longitude=lonpick, method='nearest', tolerance=0.7)

    outu_srfc = uu.sel(pressure=925).values
    outu_mid = uu.sel(pressure=650).values
    outshear = uu.sel(pressure=650).values-uu.sel(pressure=925).values
    outq = qq.sel(pressure=925).values/100
    tmid = tt.sel(pressure=650).values/100
    tsrfc = tt.sel(pressure=925).values/100

    if np.isnan(tsrfc):
        print('TSRFC IS NAN, RETURN!!!')
        return


    pes = units.Quantity(650, 'hPa')
    tes = units.Quantity(tmid+273.15, 'K')
    #ipdb.set_trace()
    theta_low = u_met.theta_e(925, tsrfc, outq/1000.)

    theta_es = calc.saturation_equivalent_potential_temperature(pes, tes)

    theta = theta_low-(np.array(theta_es)-273.15)

    tgrad = dic.attrs['Tgrad']
    tdiff = dic.attrs['Tgradbox']

    print(theta)

    #ipdb.set_trace()

    mask = np.isfinite(outp) & (outt<=t_thresh)
    antimask = ~mask

    if np.sum(mask) < 3:
        return

    for var in [outt,outp]:
        var[antimask] = np.nan

    out['tmin'] = np.nanmin(outt)
    out['tmean'] = np.mean(outt[mask])
    out['tgrad'] = tgrad
    out['tdiff'] = tdiff

    out['pmax'] = np.nanmean(ua.cut_kernel(outp,maxpos[1], maxpos[0],1)) # degrade rainfall to 30km

    out['pmean'] = np.mean(outp[mask])
    out['qmax'] = outq #np.max(outq[mask]) #30km q and shear

    out['umax_srfc'] = outu_srfc
    out['umin_mid'] = outu_mid
    out['shearmin'] =  outshear


    out['thetamax'] = theta
    out['tcwv'] = tcwvout.values

    out['tmidmax'] = tmid
    out['tsrfcmax'] = tsrfc

    out['pgt30'] = np.sum(outp[mask]>30)
    out['isvalid'] = np.sum(mask)
    out['pgt01'] = np.sum(outp[mask]>0.1)

    out['p'] = outp[mask]
    out['t'] = outt[mask]

    dic.close()

    #ipdb.set_trace()

    return out
