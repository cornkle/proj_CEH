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

import pickle as pkl


matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)


def diurnal_loop():

    for l in np.arange(17,23):
        print('Doing '+str(l))
        composite(l)

def composite(h, eh):
    #pool = multiprocessing.Pool(processes=8)


    file = cnst.MCS_CENTRE70

    hour = h

    msg = xr.open_dataarray(file)

    msg = msg[(msg['time.hour'] == hour ) & (msg['time.minute'] == 0) & (
        msg['time.year'] >= 2006) & (msg['time.year'] <= 2010) & (msg['time.month'] >=6) ]

    msg = msg.sel(lat=slice(10.9,19), lon=slice(-9.8,9.8))
    msg.attrs['eh'] = eh
    msg.attrs['refhour'] = h
    dic = {}
    #ipdb.set_trace()
    # for ids in range(0,len(msg), 50):
    #     dic = u_parallelise.era_run_arrays(1,file_loop,msg[ids:ids+50], dic) #'rano', 'rregional', 'rcnt',


    res = []
    for mm in msg:
        out =file_loop(mm)
        res.append(out)

    print('Returned from parallel')

    res = [x for x in res if x is not None]

    rres = []
    dic_names = (res[0])[1]
    for r in res:
        rres.append(np.array(r[0]))


    vars = np.array(rres)
    for id, l in enumerate(dic_names):
            dic[l] = np.nansum(np.squeeze(vars[:,id,...]), axis=0)
    # for k in dic.keys():
    #    dic[k] = np.nansum(dic[k], axis=0)


    pkl.dump(dic, open(cnst.network_data + "figs/LSTA-bullshit/AGU/composite_backtrack"+str(eh) + "UTCERA"+str(hour).zfill(2)+".p", "wb"))
    print('Dumped file')



def cut_kernel(xpos, ypos, arr, dist, probs=False, probs2=False):

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


    return kernel,  cnt, cnt2, vdic, prob, cnt3

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

    if ehour > 0:
        edate = date + pd.Timedelta(str(ehour) + ' hours')
    else:
        edate = date - pd.Timedelta(str(np.abs(ehour)) + ' hours')
    #edate = edate.replace(hour=ehour)


    t1 = edate

    file = cnst.ERA5

    try:
        cmm = xr.open_dataset(file + 'pressure_levels/ERA5_'+str(edate.year)+'_' + str(edate.month).zfill(2) + '_pl.nc')
    except:
        return None

    try:
        css = xr.open_dataset(file + 'surface/ERA5_'+str(edate.year)+'_' + str(edate.month).zfill(2) + '_srfc.nc')
    except:
        return None


    pl_clim = xr.open_dataset(file + 'CLIM/ERA5_2006-2010_CLIM_'+str(edate.month)+'-'+str(edate.hour)+'_pl.nc')
    srfc_clim = xr.open_dataset(file + 'CLIM/ERA5_2006-2010_CLIM_'+str(edate.month)+'-'+str(edate.hour)+'_srfc.nc')

    cmm = cmm.sel(time=t1)
    css = css.sel(time=t1)
    # pl_clim = pl_clim.sel(month=cmm['time.month'].values)
    # srfc_clim = srfc_clim.sel(month=cmm['time.month'].values)

    cm = cmm['t'].sel(level=925).squeeze() - pl_clim['t'].sel(level=925).squeeze() #* 1000

    cm = cm.to_dataset()

    shear =  (cmm['u'].sel(level=600).squeeze() - cmm['u'].sel(level=925).squeeze() ) #- (pl_clim['u'].sel(level=600).squeeze() - pl_clim['u'].sel(level=925).squeeze() ) #
    rh = (cmm['r'].sel(level=925).squeeze() - pl_clim['r'].sel(level=925).squeeze() )

    vwind_srfc = cmm['v'].sel(level=925).squeeze() - pl_clim['v'].sel(level=925).squeeze()
    uwind_srfc = cmm['u'].sel(level=925).squeeze() - pl_clim['u'].sel(level=925).squeeze()
    wwind_srfc = cmm['w'].sel(level=400).squeeze() - pl_clim['w'].sel(level=400).squeeze()
    div = cmm['d'].sel(level=925).squeeze()

    cape = css['cape'].squeeze() - srfc_clim['cape'].squeeze()
    surface_pressure = css['sp'].squeeze() / 100 - srfc_clim['sp'].squeeze() / 100
    #sl_pressure = css['msl'].squeeze() / 100 - srfc_clim['msl'].squeeze()/ 100
    lv = 2.26 * 1e6  # energy of evaporation in J / kg of water
    sh = css['ishf'].squeeze() * (-1)
    #lh = css['ie'].squeeze()* lv * (-1)
    #ef = lh / (sh + lh)
    #ef_clim = (srfc_clim['ie'].squeeze()*lv*(-1)) / (srfc_clim['ishf'].squeeze()*(-1) + srfc_clim['ie'].squeeze()*lv*(-1))

    sh_f = sh -  srfc_clim['ishf'].squeeze()*(-1)

    #ef_ano = ef - ef_clim

    #t2 = css['t2m'].squeeze() - srfc_clim['t2m'].squeeze()
    q = cmm['q'].sel(level=925).squeeze() - pl_clim['q'].sel(level=925).squeeze()

    cm['shear'] = shear
    cm['u925'] = uwind_srfc
    cm['v925'] = vwind_srfc

    cm['u925raw'] = cmm['u'].sel(level=925).squeeze()
    cm['v925raw'] = cmm['v'].sel(level=925).squeeze()
    cm['cape'] = cape
    cm['rh'] = rh
    cm['sf'] = surface_pressure

    cm['tciw'] = css['tciw'].squeeze() - srfc_clim['tciw'].squeeze()
    cm['sh'] = sh_f
    #cm['t2'] = t2
    cm['div'] = div *1000
    cm['q'] = q
    cm['w400'] = wwind_srfc
    srfc_clim.close()
    pl_clim.close()
    css.close()
    cmm.close()
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
    xout.name = 'probs'
    xout.values = out

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

    grad = np.gradient(ttopo.values)
    gradsum = abs(grad[0])+abs(grad[1])


    lsta_da = lsta['LSTA'].squeeze()
    if (np.sum(np.isfinite(lsta_da)) / lsta_da.size) < 0.05:
        print('Not enough valid')
        return None

    lsta_da.values[ttopo.values>=500] = np.nan
    lsta_da.values[gradsum>30] = np.nan
    pos = np.where((fi.values == 1)) #(fi.values >= 5) & (fi.values < 65)

    topo.close()

    if (np.sum(pos) == 0):
        print('No blobs found')
        return None


    dist = 300

    kernel2_sum = np.zeros((dist*2+1, dist*2+1))
    cnt_sum = np.zeros((dist*2+1, dist*2+1))
    cntp_sum = np.zeros((dist*2+1, dist*2+1))
    cntm_sum = np.zeros((dist*2+1, dist*2+1))
    probm_sum = np.zeros((dist*2+1, dist*2+1))

    edic = {}

    probs = get_previous_hours(date, fi.attrs['eh'], fi.attrs['refhour'])
    print(probs)
    try:
        probs_on_lsta = lsta.salem.transform(probs)
    except RuntimeError:
        return None

    probs_msg = get_previous_hours_msg(date, fi.attrs['eh'], fi.attrs['refhour'])
    probsm_on_lsta = lsta.salem.transform(probs_msg, interp='nearest')
    # if np.sum(probsm_on_lsta>1) != 0:
    #     'Stopp!!'
    #     pdb.set_trace()

    for y, x in zip(pos[0], pos[1]):

        lat = fi['lat'][y]
        lon = fi['lon'][x]

        point = lsta_da.sel(lat=lat, lon=lon, method='nearest')
        plat = point['lat'].values
        plon = point['lon'].values

        xpos = np.where(lsta_da['lon'].values == plon)
        xpos = int(xpos[0])
        ypos = np.where(lsta_da['lat'].values == plat)
        ypos = int(ypos[0])
        try:
            kernel2, cnt, cntp, vdic, probm, cntm = cut_kernel(xpos, ypos, lsta_da, dist, probs=probs_on_lsta, probs2=probsm_on_lsta)
        except TypeError:
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

    outlist = [kernel2_sum, cnt_sum, cntp_sum, cntm_sum, probm_sum]
    outnames = ['lsta',  'cnt', 'cntp', 'cntm', 'probmsg']
    for ek in edic.keys():
        outnames.append(ek)
        outlist.append(edic[ek])

    print('Returning')
    lsta.close()

    return outlist, outnames


def plot_gewex(h, ehour):
    hour=h

    tag = str(ehour)


    dic = pkl.load(open(cnst.network_data + "figs/LSTA-bullshit/AGU/composite_backtrack"+tag+"UTCERA"+str(hour).zfill(2)+".p", "rb"))
    #dic2 = pkl.load(open(cnst.network_data + "figs/LSTA-bullshit/AGU/CMORPH_6-6UTC_antecedent/composite_backtrack_CMORPH_"+str(hour).zfill(2)+".p", "rb"))
    extent = (dic['lsta'].shape[1]-1)/2
    xlen = dic['lsta'].shape[1]
    ylen = dic['lsta'].shape[0]

    xv, yv = np.meshgrid(np.arange(ylen), np.arange(xlen))
    st=30
    xquiv = xv[4::st, 4::st]
    yquiv = yv[4::st, 4::st]

    u = (dic['u925']/ dic['cntp'])[4::st, 4::st]
    v = (dic['v925']/ dic['cntp'])[4::st, 4::st]

    ur = (dic['u925raw']/ dic['cntp'])[4::st, 4::st]
    vr = (dic['v925raw']/ dic['cntp'])[4::st, 4::st]


    f = plt.figure(figsize=(15, 8))
    ax = f.add_subplot(231)

    plt.contourf((dic['lsta'] / dic['cnt']), cmap='RdBu_r', levels=[-1, -0.8, -0.6,-0.4,-0.2,0.2,0.4,0.6,  0.8, 1], extend='both') #-(rkernel2_sum / rcnt_sum)
    #plt.plot(extent, extent, 'bo')
    plt.colorbar(label='K')
    #pdb.set_trace()

    contours = plt.contour((dic['probmsg']/ dic['cntm'])*100, extend='both', levels=np.arange(10,70,10), cmap='jet') # #, levels=np.arange(1,5, 0.5)
    plt.clabel(contours, inline=True, fontsize=11, fmt='%1.0f')
    plt.plot(extent, extent, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 9) ) * 3, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.title('17-19UTC | '+str(np.max(dic['cnt']))+' cores, LSTA & 06-06UTC antecedent rain', fontsize=9)


    ax1 = f.add_subplot(232)
    #ipdb.set_trace()
    plt.contourf(((dic['sh'])/ dic['cntp']), extend='both',  cmap='RdBu_r',levels=[-10, -5, -2.5, -1, 1,2.5, 5,10]) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='s-1')
    contours = plt.contour((dic['t'] / dic['cntp']), extend='both',levels=[ -1, -0.8, -0.6,-0.4,-0.2,0.2,0.4,0.6,  0.8, 1], cmap='PuOr_r') #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.1f')
    #qu = ax1.quiver(xquiv, yquiv, u, v, scale=50)
    plt.plot(extent, extent, 'bo')
    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('Shading: MSG LSTA, Contours: ERA5 T2m', fontsize=9)

    ax1 = f.add_subplot(233)
    plt.contourf(((dic['q'])/ dic['cntp'])*1000, extend='both',  cmap='RdBu', levels=np.linspace(-1,1,7 )) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='kg kg-1')
    contours = plt.contour((dic['shear'] / dic['cntp']), extend='both',levels=np.arange(-17,-12,0.5), cmap='viridis') #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    plt.plot(extent, extent, 'bo')
    #qu = ax1.quiver(xquiv, yquiv, u, v, scale=50)
    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('Shading: 950hPa q anomaly, Contours: 600hPa-925hPa wind shear ', fontsize=9)

    ax1 = f.add_subplot(234)
    plt.contourf(((dic['cape'])/ dic['cnt']), extend='both',  cmap='RdBu_r', levels=np.linspace(-500,500,10)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='s-1')
    plt.plot(extent, extent, 'bo')
    contours = plt.contour((dic['sf'] / dic['cntp'])*1000, extend='both', cmap='viridis') #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    qu = ax1.quiver(xquiv, yquiv, u, v, scale=15)
    qk = plt.quiverkey(qu, 0.7, 0.95,2, '2 m s$^{-1}$',
                       labelpos='E', coordinates='figure')

    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('Shading: MSG LSTA, vectors: 900hPa wind anomaly', fontsize=9)

    ax1 = f.add_subplot(235)
    plt.contourf(((dic['rh'])/ dic['cntp']), extend='both',  cmap='RdBu', levels=np.arange(-10,11, 2)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='m s-1')

    contours = plt.contour((dic['v925']/ dic['cntp']), extend='both', cmap='viridis', levels=np.arange(-1.5,1.51, 0.5)) #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    qu = ax1.quiver(xquiv, yquiv, ur, vr, scale=30)

    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')

    # ax1 = f.add_subplot(236)
    # plt.contourf(((dic['div'])/ dic['cntp']), extend='both',  cmap='RdBu') #  tciw  #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    # plt.colorbar(label='m s-1')
    # contours = plt.contour((dic['sf'] *1000/ dic['cntp']), extend='both', cmap='viridis') #np.arange(-15,-10,0.5)
    # plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    # qu = ax1.quiver(xquiv, yquiv, u, v, scale=30)
    #
    # ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    # ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    # ax1.set_xlabel('km')
    # ax1.set_ylabel('km')

    ax1 = f.add_subplot(236)
    plt.contourf(((dic['div'])/ dic['cntp'])*100, extend='both',  cmap='RdBu', levels=np.linspace(-1,1,10)) #  tciw  #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='m s-1')
    contours = plt.contour((dic['tciw'] *100/ dic['cntp']), extend='both', cmap='viridis') #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    qu = ax1.quiver(xquiv, yquiv, u, v, scale=30)

    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')

    plt.title('ERA5 950hPa divergence (shading) & 700-950hPa wind shear', fontsize=9)


    plt.tight_layout()
    plt.show()
    plt.savefig(cnst.network_data + 'figs/LSTA-bullshit/AGU/test_' + str(ehour).zfill(2)+ 'hours_' +str(hour).zfill(2)+'UTC_single.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png')





def plot_doug(h):
    hour=h
    chour=hour
    # dic = pkl.load(open("/users/global/cornkle/figs/LSTA-bullshit/corrected_LSTA/system_scale/doug/large_scale/composite_backtrack_ERA_SEP"+str(hour).zfill(2)+".p", "rb"))
    # dic2 = pkl.load(open("/users/global/cornkle/figs/LSTA-bullshit/corrected_LSTA/system_scale/doug/CMORPH_6-6UTC_antecedent/composite_backtrack_CMORPH_SEP_"+str(chour).zfill(2)+".p", "rb"))
    dic = pkl.load(open(
        "/users/global/cornkle/figs/LSTA-bullshit/corrected_LSTA/system_scale/doug/composite_backtrack_ERA_JUN" + str(
            hour).zfill(2) + ".p", "rb"))
    dic2 = pkl.load(open(
        "/users/global/cornkle/figs/LSTA-bullshit/corrected_LSTA/system_scale/doug/composite_backtrack_CMORPH_JUN_" + str(
            chour).zfill(2) + ".p", "rb"))

    extent = (dic['lsta'].shape[1]-1)/2
    xlen = dic['lsta'].shape[1]
    ylen = dic['lsta'].shape[0]

    xv, yv = np.meshgrid(np.arange(ylen), np.arange(xlen))
    st=30
    xquiv = xv[4::st, 4::st]
    yquiv = yv[4::st, 4::st]

    u = (dic['u950']/ dic['cntp'])[4::st, 4::st]
    v = (dic['v950']/ dic['cntp'])[4::st, 4::st]


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
    contours = plt.contour((dic['t2'] / dic['cntp']), extend='both',levels=[-1, -0.8, -0.6,-0.4,-0.2,0.2,0.4,0.6,  0.8, 1], cmap='PuOr_r') #np.arange(-15,-10,0.5)
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
    plt.contourf(((dic['div'])/ dic['cntp']), extend='both',  cmap='RdBu_r') # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
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
    plt.savefig('/users/global/cornkle/figs/LSTA-bullshit/corrected_LSTA/system_scale/doug/large_scale/'+str(hour).zfill(2)+'_JUN.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png')
    plt.close()


def plot_chris(h):
    hour=h
    chour=hour
    # dic = pkl.load(open("/users/global/cornkle/figs/LSTA-bullshit/corrected_LSTA/system_scale/doug/large_scale/composite_backtrack_ERA_SEP"+str(hour).zfill(2)+".p", "rb"))
    # dic2 = pkl.load(open("/users/global/cornkle/figs/LSTA-bullshit/corrected_LSTA/system_scale/doug/CMORPH_6-6UTC_antecedent/composite_backtrack_CMORPH_SEP_"+str(chour).zfill(2)+".p", "rb"))
    dic = pkl.load(open(
        "/users/global/cornkle/figs/LSTA-bullshit/corrected_LSTA/system_scale/doug/composite_backtrack_ERA" + str(
            hour).zfill(2) + ".p", "rb"))
    dic2 = pkl.load(open(
        "/users/global/cornkle/figs/LSTA-bullshit/corrected_LSTA/system_scale/doug/CMORPH_6-6UTC_antecedent/composite_backtrack_CMORPH_" + str(
            chour).zfill(2) + ".p", "rb"))

    extent = (dic['lsta'].shape[1]-1)/2
    xlen = dic['lsta'].shape[1]
    ylen = dic['lsta'].shape[0]

    xv, yv = np.meshgrid(np.arange(ylen), np.arange(xlen))
    st=30
    xquiv = xv[4::st, 4::st]
    yquiv = yv[4::st, 4::st]

    u = (dic['u950']/ dic['cntp'])[4::st, 4::st]
    v = (dic['v950']/ dic['cntp'])[4::st, 4::st]


    f = plt.figure(figsize=(12,5))



    ax1 = f.add_subplot(121)
 #   plt.contourf(((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r', vmin=-1.5, vmax=1.5) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
 #    plt.contourf(((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r', levels=[-1.5,-1, -0.8, -0.6,-0.4,-0.2,0.2,0.4,0.6,  0.8, 1,1.5]) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
 #    plt.colorbar(label='K')
    plt.contourf(((dic['w400']) / dic['cntp']), extend='both', cmap='RdBu_r', levels=np.arange(-0.1,0.11,0.01))
    plt.colorbar(label='K')
    plt.contour(((dic['div']) / dic['cntp']), extend='both', cmap='RdBu', levels=[-0.0075,-0.005, -0.0025, 0.0025,0.005,0.0075])

    plt.plot(extent, extent, 'bo')
    # contours = plt.contour((dic['shear'] / dic['cntp']), extend='both',levels=np.arange(-17,-12,0.5), cmap='viridis') #np.arange(-15,-10,0.5)
    # plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    qu = ax1.quiver(xquiv, yquiv, u, v, scale=15)
    qk = plt.quiverkey(qu, 0.9, 0.02,1, '1 m s$^{-1}$',
                       labelpos='E', coordinates='figure')

    plt.title('Shading: LSTA, contours & vectors: 950hPa divergence & wind anomaly', fontsize=11)

    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9)  - extent) * 3, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')

    ax1 = f.add_subplot(122)
    plt.contourf(((dic['q'])*1000/ dic['cntp']), extend='both',  cmap='RdBu',levels=np.arange(-1,1.1,0.1)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='g kg-1')
    contours = plt.contour((dic['shear'] / dic['cntp']), extend='both', cmap='viridis') #np.arange(-15,-10,0.5)     contours = plt.contour((dic['shear'] / dic['cntp']), extend='both', cmap='viridis') #np.arange(-15,-10,0.5)

    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    plt.plot(extent, extent, 'bo')
    #qu = ax1.quiver(xquiv, yquiv, u, v, scale=50)
    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -extent) * 3, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('Shading: 950hPa q anomaly, Contours: 600hPa-925hPa wind shear ', fontsize=11)



    plt.tight_layout()
    plt.savefig('/users/global/cornkle/figs/LSTA-bullshit/corrected_LSTA/system_scale/doug/'+str(hour).zfill(2)+'_test.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png')
    plt.close()


def plot_all():

    hours = [16,17,18,19,20,21,22,23,0,1,2,3,4,5,6,7]

    for h in hours:
        plot_gewex(h)
