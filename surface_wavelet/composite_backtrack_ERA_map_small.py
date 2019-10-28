# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 10:15:40 2016

@author: cornkle
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib
import ipdb
import pandas as pd
from collections import OrderedDict
import salem
from utils import u_met, u_parallelise, u_gis, u_arrays, constants as cnst, u_grid, u_darrays
from scipy.interpolate import griddata
import multiprocessing

import pickle as pkl


matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)


def diurnal_loop():

    for l in np.arange(17,23):
        print('Doing '+str(l))
        composite(l)


def composite(h, eh):

    file = cnst.MCS_CENTRE70

    hour = h

    msg = xr.open_dataarray(file)

    #for year in np.arange(2008, 2009):
    year='all'

    msgo = msg[((msg['time.hour'] >= 17) & (msg['time.hour']<=19) ) & (msg['time.minute'] == 0) & (
        msg['time.year'] >=2006) & (msg['time.month'] >=7) ]

    msgo = msgo.sel(lat=slice(10.2,19), lon=slice(-9.9,9.9))
    msgo.attrs['eh'] = eh
    msgo.attrs['refhour'] = h

    dic = u_parallelise.era_run_arrays(2, file_loop, msgo)

    # res = []
    # for mm in msgo:
    #     out =file_loop(mm)
    #     res.append(out)
    #
    # res = [x for x in res if x is not None]
    # dic = {}
    #
    # rres = []
    # dic_names = (res[0])[1]
    # for r in res:
    #     rres.append(np.array(r[0]))
    #
    # vars = np.array(rres)
    # for id, l in enumerate(dic_names):
    #     try:
    #         dic[l] = dic[l] + np.nansum(np.squeeze(vars[:, id, ...]), axis=0)
    #     except KeyError:
    #         dic[l] = np.nansum(np.squeeze(vars[:, id, ...]), axis=0)

    # return
    # print('Returned from parallel')
    #

    pkl.dump(dic, open(cnst.network_data + "figs/LSTA/corrected_LSTA/new/composite_backtrack"+str(eh) + "UTCERA"+str(hour).zfill(2)+'_'+str(year)+"_stormCentre.p", "wb"))
    print('Dumped file')



def cut_kernel(xpos, ypos, arr, dist, probs=False, probs2=False, probs3=False):

    kernel = u_arrays.cut_kernel(arr,xpos, ypos,dist)

    vdic = {}

    for d in probs.data_vars:

        var = u_arrays.cut_kernel(probs[d].values,xpos, ypos,dist)
        vdic[d] = var

    cnt2 = np.zeros_like(kernel)
    cnt2[np.isfinite(vdic[list(vdic.keys())[0]])] = 1


    if (np.sum(np.isfinite(kernel)) < 2):
        return

    cnt = np.zeros_like(kernel)
    cnt[np.isfinite(kernel)] = 1


    if kernel.shape != (dist*2+1, dist*2+1):
        print('Kernels shape wrong!')
        ipdb.set_trace()

    kernel = kernel - np.nanmean(kernel)


    if np.nansum(probs2) > 0:
        prob = u_arrays.cut_kernel(probs2,xpos, ypos,dist)
        cnt3 = np.zeros_like(kernel)
        cnt3[np.isfinite(prob)] = 1

    else:
        prob = np.zeros_like(kernel)
        cnt3 = np.zeros_like(kernel)


    if np.nansum(probs3) > 0:
        probcm = u_arrays.cut_kernel(probs3,xpos, ypos,dist)
        cnt4 = np.zeros_like(kernel)
        cnt4[np.isfinite(probcm)] = 1

    else:
        probcm = np.zeros_like(kernel)
        cnt4 = np.zeros_like(kernel)


    return kernel,  cnt, cnt2, vdic, prob, cnt3, probcm, cnt4




def get_previous_hours(date, ehour, refhour):

    #ehour = 3
    #
    # if (date.hour) <= 16:
    #     print('Nighttime')
    #     edate = date - pd.Timedelta('1 days')
    #
    #     edate = edate.replace(hour=ehour)
    #
    # else:

    date = date.replace(hour=refhour)
    cm = xr.Dataset()

    if ehour > 0:
        edate = date + pd.Timedelta(str(ehour) + ' hours')
    else:
        edate = date - pd.Timedelta(str(np.abs(ehour)) + ' hours')
    #edate = edate.replace(hour=ehour)


    t1 = edate

    file = cnst.ERA5

    try:
        cmm = xr.open_dataset(file + 'hourly/pressure_levels/ERA5_'+str(edate.year)+'_' + str(edate.month).zfill(2) + '_pl.nc')
        cmm = u_darrays.flip_lat(cmm)
    except:
        return None

    try:
        css = xr.open_dataset(file + 'hourly/surface/ERA5_'+str(edate.year)+'_' + str(edate.month).zfill(2) + '_srfc.nc')
        css = u_darrays.flip_lat(css)
    except:
        return None


    pl_clim = xr.open_dataset(file + 'monthly/synop_selfmade/CLIM_2006-2010_new/ERA5_2006-2010_CLIM_'+str(edate.month).zfill(2)+'-'+str(edate.hour).zfill(2)+'_pl_rw.nc').load()
    pl_clim = u_darrays.flip_lat(pl_clim)
    srfc_clim = xr.open_dataset(file + 'monthly/synop_selfmade/CLIM_2006-2010_new/ERA5_2006-2010_CLIM_'+str(edate.month).zfill(2)+'-'+str(edate.hour).zfill(2)+'_srfc_rw.nc').load()
    #srfc_clim = u_darrays.flip_lat(srfc_clim)

    ## latitude in surface is already flipped, not for pressure levels though... ?!


    cmm = cmm.sel(longitude=slice(-13, 13), latitude=slice(8, 21))
    css = css.sel(longitude=slice(-13, 13), latitude=slice(8, 21))
    pl_clim = pl_clim.sel(longitude=slice(-13, 13), latitude=slice(8, 21))
    srfc_clim = srfc_clim.sel(longitude=slice(-13, 13), latitude=slice(8, 21))


    cmm = cmm.sel(time=t1)
    pl_clim = pl_clim.squeeze()

    css = css.sel(time=t1)
    srfc_clim = srfc_clim.squeeze()

    t = cmm['t'].sel(level=925).squeeze() - pl_clim['t'].sel(level=925).squeeze() #* 1000

    shear =  (cmm['u'].sel(level=650).squeeze() - cmm['u'].sel(level=925).squeeze() ) #- (pl_clim['u'].sel(level=600).squeeze() - pl_clim['u'].sel(level=925).squeeze() ) #

    vwind_srfc = cmm['v'].sel(level=925).squeeze() - pl_clim['v'].sel(level=925).squeeze()
    uwind_srfc = cmm['u'].sel(level=925).squeeze() - pl_clim['u'].sel(level=925).squeeze()
    #wwind_srfc = cmm['w'].sel(level=400).squeeze() - pl_clim['w'].sel(level=400).squeeze()
    div = cmm['d'].sel(level=925).squeeze()


    surface_pressure = css['sp'].squeeze() / 100 - srfc_clim['sp'].squeeze() / 100
    # #sl_pressure = css['msl'].squeeze() / 100 - srfc_clim['msl'].squeeze()/ 100
    # lv = 2.26 * 1e6  # energy of evaporation in J / kg of water
    sh = css['ishf'].squeeze() * (-1)
    #lh = css['ie'].squeeze()* lv * (-1)
    #ef = lh / (sh + lh)
    #ef_clim = (srfc_clim['ie'].squeeze()*lv*(-1)) / (srfc_clim['ishf'].squeeze()*(-1) + srfc_clim['ie'].squeeze()*lv*(-1))

    sh_f = sh -  srfc_clim['ishf'].squeeze()*(-1)

    #ef_ano = ef - ef_clim

    q = cmm['q'].sel(level=925).squeeze() - pl_clim['q'].sel(level=925).squeeze()

    cm['shear'] = shear
    cm['u925'] = uwind_srfc
    cm['v925'] = vwind_srfc

    cm['u925raw'] = cmm['u'].sel(level=925).squeeze()
    cm['v925raw'] = cmm['v'].sel(level=925).squeeze()
    cm['slp'] = surface_pressure

    cm['tciw'] = css['tciw'].squeeze() - srfc_clim['tciw'].squeeze()
    cm['sh'] = sh_f
    cm['div'] = div *1000
    cm['q'] = q
    cm['t'] = t

    del srfc_clim
    del pl_clim
    del css
    del cmm
    return cm

def get_previous_hours_msg(date, ehour, refhour):

    # tdic = {18 : ('36 hours', '15 hours'),
    #         19 : ('37 hours', '16 hours'),
    #         20: ('38 hours', '17 hours'),
    #         21: ('39 hours', '18 hours'),
    #         22: ('40 hours', '19 hours'),
    #         23: ('41 hours', '20 hours'),
    #         0: ('42 hours', '21 hours'),
    #         3: ('45 hours', '24 hours'),
    #         6: ('48 hours', '27 hours')}
    # before = pd.Timedelta(tdic[date.hour][0])
    # before2 = pd.Timedelta(tdic[date.hour][1])
    date = date.replace(hour=refhour)

    if ehour > 0:
        edate = date + pd.Timedelta(str(ehour) + ' hours')
    else:
        edate = date - pd.Timedelta(str(np.abs(ehour)) + ' hours')
        #edate = edate.replace(hour=ehour)

    t1 = edate - pd.Timedelta('1 hours')
    t2 = edate + pd.Timedelta('1 hours')

    file = cnst.MCS_15K# MCS_15K #_POINTS_DOM
    msg = xr.open_dataarray(file)
    try:
        msg = msg.sel(time=slice(t1.strftime("%Y-%m-%dT%H"), t2.strftime("%Y-%m-%dT%H")))
    except OverflowError:
        return None

    #print(prev_time.strftime("%Y-%m-%dT%H"), date.strftime("%Y-%m-%dT%H"))
    pos = np.where((msg.values <= -40) ) #(msg.values >= 5) & (msg.values < 65)) # #

    out = np.zeros_like(msg)
    out[pos] = 1
    out = np.sum(out, axis=0)
    out[out>0]=1
    # if np.sum(out>1) != 0:
    #     'Stop!!!'
    #     pdb.set_trace()

    msg = msg.sum(axis=0)*0

    xout = msg.copy()
    del msg
    xout.name = 'probs'
    xout.values = out

    return xout


def get_previous_hours_CMORPH(date):


    tdic = {16 : ('34 hours', '6 hours'),  # 6am prev - 9am storm day
            17: ('35 hours', '7 hours'),
            18 : ('36 hours', '8 hours'),
            19 : ('37 hours', '9 hours'),
            20: ('38 hours', '10 hours'),
            21: ('39 hours', '11 hours'),
            22: ('40 hours', '12 hours'),
            23: ('41 hours', '13 hours'),
            0: ('42 hours', '14 hours'),
            1: ('43 hours', '15 hours'),
            2: ('44 hours', '16 hours'),
            3: ('45 hours', '17 hours'),
            6: ('48 hours', '18 hours')}
    before = pd.Timedelta(tdic[date.hour][0])
    before2 = pd.Timedelta(tdic[date.hour][1])

    t1 = date - before
    t2 = date - before2

    # before2 = pd.Timedelta('15 minutes')
    #
    # t1 = date #- before
    # t2 = date + before2

    file = cnst.CMORPH
    try:
        cmm = xr.open_dataarray(file + 'CMORPH_WA_' + str(date.year) + '.nc')
    except:
        return None
    cmm = cmm.sel( time=slice(t1, t2)).sum(dim='time') # lat=slice(10.9,19), lon=slice(-9.8,9.8),
    # cmm = cmm.sel(lat=slice(10.9, 19), lon=slice(-9.8, 9.8))
    # posi = np.where(pd.to_datetime(cmm.time.values) == date)
    #
    # cmm = cmm.isel(time=slice(posi[0][0]-2,posi[0][0]))
    cm = cmm
    pos = np.where(cm.values>=10) #(msg.values >= 5) & (msg.values < 65)) # #

    out = np.zeros_like(cm)
    out[pos] = 1
    #out = np.sum(out, axis=0) / out.shape[0]

    xout = cm.copy()
    xout.name = 'probs'
    xout.values = out

    del cm

    return xout



def file_loop(fi):

    print('Doing day: ', fi.time.values)

    date = pd.Timestamp(fi.time.values)

    dayd = pd.Timedelta('1 days')
    if fi['time.hour'].values.size != 1:
        'hour array too big, problem!'

    if (fi['time.hour'].values) <= 16:
        print('Nighttime')
        daybefore = date - dayd
    else:
        print('Daytime')
        daybefore = date

    fdate = str(daybefore.year) + str(daybefore.month).zfill(2) + str(daybefore.day).zfill(2)

    try:
        lsta = xr.open_dataset(cnst.LSTA_NEW + 'lsta_daily_' + fdate + '.nc')
    except OSError:
        return None
    print('Doing '+ 'lsta_daily_' + fdate + '.nc')

    topo = xr.open_dataset(cnst.LSTA_TOPO)
    ttopo = topo['h']
    #
    grad = np.gradient(ttopo.values)
    gradsum = abs(grad[0])+abs(grad[1])


    lsta_da = lsta['LSTA'].squeeze()
    if (np.sum(np.isfinite(lsta_da)) / lsta_da.size) < 0.05:
        print('Not enough valid')
        return None

    lsta_da.values[ttopo.values>=450] = np.nan
    lsta_da.values[gradsum>30] = np.nan
    pos = np.where((fi.values == 1)) #(fi.values >= 5) & (fi.values < 65)

    del topo

    if (np.sum(pos) == 0):
        print('No blobs found')
        return None


    dist = 200

    kernel2_sum = np.zeros((dist*2+1, dist*2+1))
    cnt_sum = np.zeros((dist*2+1, dist*2+1))
    cntp_sum = np.zeros((dist*2+1, dist*2+1))
    cntm_sum = np.zeros((dist*2+1, dist*2+1))
    probm_sum = np.zeros((dist*2+1, dist*2+1))

    edic = {}

    probs = get_previous_hours(date, fi.attrs['eh'], fi.attrs['refhour'])
    print('Era5 collect')

    try:
        probs_on_lsta = lsta.salem.transform(probs)
    except RuntimeError:
        print('Era5 on LSTA interpolation problem')
        return None
    del probs
    probs_msg = get_previous_hours_msg(date, fi.attrs['eh'], fi.attrs['refhour'])

    probsm_on_lsta = lsta.salem.transform(probs_msg, interp='nearest')
    del probs_msg

    probs_cm = get_previous_hours_CMORPH(date)
    try:
        probscm_on_lsta = lsta.salem.transform(probs_cm)
    except RuntimeError:
        return None
    del probs_cm
    del lsta

    mcs_hour = xr.open_dataarray(cnst.MCS_HOUR_DAILY)
    for y, x in zip(pos[0], pos[1]):

        lat = fi['lat'][y]
        lon = fi['lon'][x]

        mhour = mcs_hour.sel(time=pd.datetime(int(fdate[0:4]), int(fdate[4:6]), int(fdate[6:8]), 16), lon=lon,
                             lat=lat).values
        if mhour < 16:
            mhour += 24

        if mhour == 0:
            mhour += 24

        chour = fi['time.hour'].values

        # ipdb.set_trace()

        if (chour >= 0) & (chour <= 15):
            chour += 24
        if (mhour < fi['time.hour'].values) | (np.isnan(mhour)):
            print('Core overlaps: earliest:', mhour, ' core: ', chour)
            continue

        point = lsta_da.sel(lat=lat, lon=lon, method='nearest')
        plat = point['lat'].values
        plon = point['lon'].values

        xpos = np.where(lsta_da['lon'].values == plon)
        xpos = int(xpos[0])
        ypos = np.where(lsta_da['lat'].values == plat)
        ypos = int(ypos[0])
        try:
            kernel2, cnt, cntp, vdic, probm, cntm, probcm, cntc = cut_kernel(xpos, ypos, lsta_da, dist, probs=probs_on_lsta, probs2=probsm_on_lsta, probs3=probscm_on_lsta)
        except TypeError:
            continue

        if np.nansum(probm[:,dist::])>=2:   # filter out cases with MCSs at 12
            print('MCS continue')
            continue

        if np.nanmax(vdic['tciw'][:,dist::])>=0.05:   # filter out cases with MCSs at 12
            print('MCS continue')
            continue

        if np.nansum(probcm[100:200,:])>=100:   # filter out cases with rainfall to the south
            print('Southward rainfall, continue')
            continue



        kernel2_sum = np.nansum(np.stack([kernel2_sum, kernel2]), axis=0)
        cntp_sum = np.nansum(np.stack([cntp_sum, cntp]), axis=0)
        cnt_sum = np.nansum(np.stack([cnt_sum, cnt]), axis=0)

        cntm_sum = np.nansum(np.stack([cntm_sum, cntm]), axis=0)
        probm_sum = np.nansum(np.stack([probm_sum, probm]), axis=0)

        for ks in vdic.keys():
            if ks in edic:
                edic[ks] = np.nansum(np.stack([edic[ks], vdic[ks]]), axis=0)
            else:
                edic[ks] = vdic[ks]

    if np.sum(cnt_sum)==0:
        return None

    outlist = [kernel2_sum, cnt_sum, cntp_sum, cntm_sum, probm_sum]
    outnames = ['lsta',  'cnt', 'cntp', 'cntm', 'probmsg']
    for ek in edic.keys():
        outnames.append(ek)
        outlist.append(edic[ek])

    print('Returning')
    #lsta.close()

    return outlist, outnames



def plot_doug(h, eh):

    # dic = pkl.load(open("/users/global/cornkle/figs/LSTA-bullshit/corrected_LSTA/system_scale/doug/large_scale/composite_backtrack_ERA_SEP"+str(hour).zfill(2)+".p", "rb"))
    # dic2 = pkl.load(open("/users/global/cornkle/figs/LSTA-bullshit/corrected_LSTA/system_scale/doug/CMORPH_6-6UTC_antecedent/composite_backtrack_CMORPH_SEP_"+str(chour).zfill(2)+".p", "rb"))
    dic = pkl.load(open(
        cnst.network_data + "figs/LSTA/corrected_LSTA/new/composite_backtrack"+str(eh) + "UTCERA"+str(h).zfill(2)+"_all_small.p", "rb"))
    dic2 = pkl.load(open(
        cnst.network_data + "figs/LSTA/corrected_LSTA/new/composite_backtrack_CMORPH_" + str(h).zfill(2) + ".p", "rb"))

    extent = (dic['lsta'].shape[1]-1)/2
    xlen = dic['lsta'].shape[1]
    ylen = dic['lsta'].shape[0]

    xv, yv = np.meshgrid(np.arange(ylen), np.arange(xlen))
    st=30
    xquiv = xv[4::st, 4::st]
    yquiv = yv[4::st, 4::st]

    u = (dic['u925']/ dic['cntp'])[4::st, 4::st]
    v = (dic['v925']/ dic['cntp'])[4::st, 4::st]


    f = plt.figure(figsize=(10,8))
    ax = f.add_subplot(221)

    plt.contourf((dic['lsta'] / dic['cnt']), cmap='RdBu_r', levels=[-1, -0.8, -0.6,-0.4,-0.2,0.2,0.4,0.6,  0.8, 1], extend='both') #-(rkernel2_sum / rcnt_sum)
    #plt.plot(extent, extent, 'bo')
    plt.colorbar(label='K')
    #pdb.set_trace()

    contours = plt.contour((dic2['prob']/ dic2['cntp'])*100, extend='both', levels=np.arange(10,70,10), cmap='jet') # #, levels=np.arange(1,5, 0.5)
    plt.clabel(contours, inline=True, fontsize=11, fmt='%1.0f')
    plt.plot(extent, extent, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6 , dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.title('23-01UTC | '+str(np.max(dic['cnt']))+' cores, LSTA & 06-06UTC antecedent rain', fontsize=9)


    ax1 = f.add_subplot(222)
    plt.contourf(((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r',levels=[-1, -0.8, -0.6,-0.4,-0.2,0.2,0.4,0.6,  0.8, 1]) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='K')
    contours = plt.contour((dic['t'] / dic['cntp']), extend='both',levels=[-1, -0.8, -0.6,-0.4,-0.2,0.2,0.4,0.6,  0.8, 1], cmap='PuOr_r') #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.1f')
    #qu = ax1.quiver(xquiv, yquiv, u, v, scale=50)
    plt.plot(extent, extent, 'bo')
    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('Shading: MSG LSTA, Contours: ERA5 T2m', fontsize=9)

    ax1 = f.add_subplot(223)
    plt.contourf(((dic['q'])*1000/ dic['cntp']), extend='both',  cmap='RdBu',levels=np.arange(-1,1.1,0.1)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='g kg-1')
    contours = plt.contour((dic['shear'] / dic['cntp']), extend='both',levels=np.arange(-18,-12,0.5), cmap='viridis') #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    plt.plot(extent, extent, 'bo')
    #qu = ax1.quiver(xquiv, yquiv, u, v, scale=50)
    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('Shading: 950hPa q anomaly, Contours: 600hPa-925hPa wind shear ', fontsize=9)

    ax1 = f.add_subplot(224)
 #   plt.contourf(((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r', vmin=-1.5, vmax=1.5) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.contourf(((dic['div'])/ dic['cntp'])*100, extend='both',  cmap='RdBu', levels=np.linspace(-1,1,10)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='K')
    plt.plot(extent, extent, 'bo')
    # contours = plt.contour((dic['shear'] / dic['cntp']), extend='both',levels=np.arange(-17,-12,0.5), cmap='viridis') #np.arange(-15,-10,0.5)
    # plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    qu = ax1.quiver(xquiv, yquiv, u, v, scale=15)
    qk = plt.quiverkey(qu, 0.9, 0.02,1, '1 m s$^{-1}$',
                       labelpos='E', coordinates='figure')

    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('Shading: MSG LSTA, vectors: 950hPa wind anomaly', fontsize=9)


    plt.tight_layout()
    plt.show()
    #plt.savefig('/users/global/cornkle/figs/LSTA-bullshit/corrected_LSTA/system_scale/doug/large_scale/'+str(hour).zfill(2)+'_JUN.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png')
    #plt.close()



def plot_all():

    hours = [16,17,18,19,20,21,22,23,0,1,2,3,4,5,6,7]

    for h in hours:
        plot_doug(h)
