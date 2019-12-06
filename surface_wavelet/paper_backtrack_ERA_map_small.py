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
from utils import u_met, u_parallelise, u_arrays as ua, constants as cnst, u_darrays
from scipy.interpolate import griddata
import multiprocessing
import glob

import pickle as pkl


matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)


def diurnal_loop():
    afternoon = list(range(14,24))
    night = list(range(0,8))
    all = night #afternoon + night

    hlist = []
    for hh in all:
        if hh >= 14:
            hlist.append((hh,12-hh))
        else:
            hlist.append((hh, 12-(hh+24)))

    for l in hlist:
        print('Doing '+str(l))
        composite(l[0], l[1])


def eh_loop():

    #all = np.arange(-38,1,3)
    all= np.arange(4,30,3)

    #hlist = []
    # for hh in all:
    #     if hh >= 14:
    #         hlist.append((hh,12-hh))
    #     else:
    #         hlist.append((hh, 12-(hh+24)))

    for l in all:
        print('Doing '+str(l))
        composite(17, l)


def composite(h, eh):

    file = cnst.MCS_POINTS_DOM
    #file = cnst.MCS_CENTRE70

    hour = h

    msg = xr.open_dataarray(file)

    for year in np.arange(2006, 2011):
        #year='all'

        msgo = msg[((msg['time.hour'] == hour) ) & (msg['time.minute'] == 0) & (
            msg['time.year'] ==year) & ((msg['time.month'] >=6) & (msg['time.month'] <=9))  ]


        msgo = msgo.sel(lat=slice(10.2,19), lon=slice(-9.9,9.9))
        msgo.attrs['eh'] = eh
        msgo.attrs['refhour'] = h
        print('MCS dataset length:', len(msgo))
        dic = u_parallelise.era_run_arrays(4, file_loop, msgo)
        #res = []
        # for m in msgo:
        #     out = file_loop(m)
        #     res.append(out)
        #
        # ipdb.set_trace()
        # return

        pkl.dump(dic, open(cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/ERA5_composite_cores_LSTA_500w04_15k_WET_"+str(eh) + "UTCERA"+str(hour).zfill(2)+'_'+str(year)+"_small_cores.p", "wb"))
        del dic
        print('Dumped file')



def cut_kernel(xpos, ypos, arr, dist, probs=False, probs2=False, probs3=False, lsta_prev=False):

    kernel = ua.cut_kernel(arr,xpos, ypos,dist)

    vdic = {}

    for d in probs.data_vars:

        var = ua.cut_kernel(probs[d].values,xpos, ypos,dist)
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

    kernel = kernel #- np.nanmean(kernel)


    if np.nansum(probs2) > 0:
        prob = ua.cut_kernel(probs2,xpos, ypos,dist)
        cnt3 = np.zeros_like(kernel)
        cnt3[np.isfinite(prob)] = 1

    else:
        prob = np.zeros_like(kernel)
        cnt3 = np.zeros_like(kernel)

    if np.nansum(probs3) > 0:
        probcm = ua.cut_kernel(probs3,xpos, ypos,dist)
        cnt4 = np.zeros_like(kernel)
        cnt4[np.isfinite(probcm)] = 1

    else:
        probcm = np.zeros_like(kernel)
        cnt4 = np.zeros_like(kernel)

    if np.nansum(lsta_prev) > 0:
        plsta = ua.cut_kernel(lsta_prev,xpos, ypos,dist)
        lcnt = np.zeros_like(kernel)
        lcnt[np.isfinite(plsta)] = 1

    else:
        plsta = np.zeros_like(kernel)
        lcnt = np.zeros_like(kernel)

    kmean = kernel - np.nanmean(kernel)
    ycirc100e, xcirc100e = ua.draw_circle(dist+51, dist+1, 17)  # at - 150km, draw 50km radius circle
    e100 = np.nanmean(kmean[ycirc100e,xcirc100e])

    # if e100 >= -3.5:   ### random LSTA p10
    #     return
    #
    # if e100 <= 2.7:   ### random LSTA p90
    #     return

    # if e100 >= -3:  ### core LSTA p10
    #     return

    # if e100 <= 2.88:  ### random LSTA p90
    #     return

    if e100 >= -1.36:  ### core LSTA p25
        return
    #
    # if e100 <= 1.45:  ### random LSTA p75
    #     return

    if (e100 >= 1.1) | (e100<=1):
        return

    return kernel,  cnt, cnt2, vdic, prob, cnt3, probcm, cnt4,plsta, lcnt


def get_previous_hours(storm_date, lsta_date, ehour, refhour):


    date = storm_date.replace(hour=refhour)
    cm = xr.Dataset()

    if ehour > 0:
        edate = date + pd.Timedelta(str(ehour) + ' hours')
    else:
        edate = date - pd.Timedelta(str(np.abs(ehour)) + ' hours')
    #edate = edate.replace(hour=ehour)


    t1 = edate
    t1_MCScheck = lsta_date.replace(hour=12)  # no MCS at 12, otherwise allowed

    file = cnst.ERA5

    #ipdb.set_trace()

    try:
        cmp = xr.open_dataset(file + 'hourly/pressure_levels/ERA5_'+str(edate.year)+'_' + str(edate.month).zfill(2) + '_pl.nc')
        cmp = u_darrays.flip_lat(cmp)
    except:
        return None

    # try:
    #     css = xr.open_dataset(file + 'hourly/surface/ERA5_'+str(t1_MCScheck.year)+'_' + str(t1_MCScheck.month).zfill(2) + '_srfc.nc')
    #     css = u_darrays.flip_lat(css)
    # except:
    #     return None

    csm = xr.open_dataset(
        file + 'hourly/surface/ERA5_' + str(edate.year) + '_' + str(edate.month).zfill(2) + '_srfc.nc')
    csm = u_darrays.flip_lat(csm)


    pl_clim = xr.open_dataset(file + 'monthly/synop_selfmade/CLIM_2006-2010_new/ERA5_2006-2010_CLIM_'+str(edate.month).zfill(2)+'-'+str(edate.hour).zfill(2)+'_pl_rw.nc').load()
    pl_clim = u_darrays.flip_lat(pl_clim)

    esrfc_clim = xr.open_dataset(
        file + 'monthly/synop_selfmade/CLIM_2006-2010_new/ERA5_2006-2010_CLIM_' + str(edate.month).zfill(
            2) + '-' + str(edate.hour).zfill(2) + '_srfc_rw.nc').load()

    ## latitude in surface is already flipped, not for pressure levels though... ?!


    cmp = cmp.sel(longitude=slice(-13, 13), latitude=slice(8, 21))
    #css = css.sel(longitude=slice(-13, 13), latitude=slice(8, 21))
    pl_clim = pl_clim.sel(longitude=slice(-13, 13), latitude=slice(8, 21))

    #srfc_clim = srfc_clim.sel(longitude=slice(-13, 13), latitude=slice(8, 21))
    esrfc_clim = esrfc_clim.sel(longitude=slice(-13, 13), latitude=slice(8, 21))

    cmm = cmp.sel(time=t1)
    pl_clim = pl_clim.squeeze()

    css = cmp.sel(time=t1_MCScheck)
    scm = csm.sel(time=t1)

    #srfc_clim = srfc_clim.squeeze()

    sh = scm['ishf'].squeeze() - esrfc_clim['ishf'].squeeze()
    ev = scm['ie'].squeeze() - esrfc_clim['ie'].squeeze()
    skt = scm['skt'].squeeze() - esrfc_clim['skt'].squeeze()

    t = cmm['t'].sel(level=925).squeeze() - pl_clim['t'].sel(level=925).squeeze() #* 1000

    shear =  (cmm['u'].sel(level=650).squeeze() - cmm['u'].sel(level=925).squeeze() ) #- (pl_clim['u'].sel(level=600).squeeze() - pl_clim['u'].sel(level=925).squeeze() ) #

    vwind_srfc = cmm['v'].sel(level=925).squeeze() - pl_clim['v'].sel(level=925).squeeze()
    uwind_srfc = cmm['u'].sel(level=925).squeeze() - pl_clim['u'].sel(level=925).squeeze()

    uwind_up = cmm['u'].sel(level=650).squeeze() #- pl_clim['u'].sel(level=650).squeeze()
    vwind_up = cmm['v'].sel(level=650).squeeze() #- pl_clim['v'].sel(level=650).squeeze()
    #wwind_up = cmm['w'].sel(level=650).squeeze()

    uwind_up_ano = cmm['u'].sel(level=650).squeeze() - pl_clim['u'].sel(level=650).squeeze()
    vwind_up_ano = cmm['v'].sel(level=650).squeeze() - pl_clim['v'].sel(level=650).squeeze()

    div = cmm['d'].sel(level=925).squeeze()

    theta_e_diff = u_met.theta_e(925, cmm['t'].sel(level=925).squeeze().values - 273.15,
                            cmm['q'].sel(level=925).squeeze()) - \
                   u_met.theta_e(650, cmm['t'].sel(level=650).squeeze().values - 273.15,
                                 cmm['q'].sel(level=650).squeeze())

    theta_e_diff_clim = u_met.theta_e(925, pl_clim['t'].sel(level=925).squeeze().values - 273.15,
                                 pl_clim['q'].sel(level=925).squeeze()) - \
                   u_met.theta_e(650, pl_clim['t'].sel(level=650).squeeze().values - 273.15,
                                 pl_clim['q'].sel(level=650).squeeze())

    theta_e = theta_e_diff - theta_e_diff_clim

    q = cmm['q'].sel(level=925).squeeze() - pl_clim['q'].sel(level=925).squeeze()

    cm['shear'] = shear
    cm['u925'] = uwind_srfc
    cm['v925'] = vwind_srfc
    cm['v925_orig'] = cmm['v'].sel(level=925).squeeze()
    cm['u925_orig'] = cmm['u'].sel(level=925).squeeze()

    cm['tciw'] = css['w'].sel(level=350).squeeze() #- srfc_clim['tciw'].squeeze()
    cm['tciwlow'] = css['w'].sel(level=850).squeeze()
    cm['tciwmid'] = css['w'].sel(level=500).squeeze()
    #ipdb.set_trace()
    #cm['sp'] = scm['sp'].squeeze() - esrfc_clim['sp'].squeeze()

    cm['u650_orig'] = uwind_up
    cm['v650_orig'] = vwind_up

    cm['u650'] = uwind_up_ano
    cm['v650'] = vwind_up_ano
    cm['sh'] = sh
    cm['ev'] = ev
    cm['skt'] = skt

    cm['div'] = div *1000
    cm['q'] = q
    cm['t'] = t
    cm['theta_e'] = theta_e

    #del srfc_clim
    del pl_clim
    #del esrfc_clim
    del css
    #del csm
    del cmm
    return cm

def get_previous_hours_msg(storm_date, lsta_date, ehour, refhour):


    date = storm_date.replace(hour=refhour)

    if ehour > 0:
        edate = date + pd.Timedelta(str(ehour) + ' hours')
    else:
        edate = date - pd.Timedelta(str(np.abs(ehour)) + ' hours')

    edate = lsta_date.replace(hour=12)  # make 12 reference hour for MCS filter

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


    tdic = {14 : ('32 hours', '4 hours'),
            15: ('33 hours', '5 hours'),
            16 : ('34 hours', '6 hours'),  # 6am prev - 10am storm day
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
            4: ('46 hours', '17 hours'),
            5: ('47 hours', '18 hours'),
            6: ('48 hours', '19 hours'),
            7: ('49 hours', '20 hours')}

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
    cmm = cmm.sel( time=slice(t1, t2)).sum(dim='time')
    # cmm = cmm.sel(lat=slice(10.9, 19), lon=slice(-9.8, 9.8))

    cm = cmm
    pos = np.where(cm.values>=5)

    out = np.zeros_like(cm)
    out[pos] = 1
    #out = np.sum(out, axis=0) / out.shape[0]

    xout = cm.copy()
    xout.name = 'probs'
    xout.values = out

    return xout


def file_loop(fi):

    print('Doing day: ', fi.time.values)

    storm_date = pd.Timestamp(fi.time.values)

    dayd = pd.Timedelta('1 days')
    if fi['time.hour'].values.size != 1:
        'hour array too big, problem!'

    if (fi['time.hour'].values) <= 13:
        print('Nighttime')
        lsta_date = storm_date - dayd
    else:
        print('Daytime')
        lsta_date = storm_date

    prev_lsta_date = lsta_date - dayd

    fdate = str(lsta_date.year) + str(lsta_date.month).zfill(2) + str(lsta_date.day).zfill(2)
    pfdate = str(prev_lsta_date.year) + str(prev_lsta_date.month).zfill(2) + str(prev_lsta_date.day).zfill(2)

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
    pos = np.where((fi.values >= 5) & (fi.values <= 65)) #np.where((fi.values == 1)) #np.where((fi.values >= 5) & (fi.values <= 65)) #np.where((fi.values == 1)) #(fi.values >= 5) & (fi.values < 65)

    del topo

    if (np.sum(pos) == 0):
        print('No blobs found')
        return None

    try:
        lsta_prev = xr.open_dataset(cnst.LSTA_NEW + 'lsta_daily_' + pfdate + '.nc')
    except OSError:
        return None
    prev_lsta_da = lsta_prev['LSTA'].squeeze()


    dist = 200

    kernel2_sum = np.zeros((dist*2+1, dist*2+1))
    kernelp_sum = np.zeros((dist * 2 + 1, dist * 2 + 1))
    cnt_sum = np.zeros((dist*2+1, dist*2+1))
    lcnt_sum = np.zeros((dist * 2 + 1, dist * 2 + 1))
    cntp_sum = np.zeros((dist*2+1, dist*2+1))
    cntm_sum = np.zeros((dist*2+1, dist*2+1))
    probm_sum = np.zeros((dist*2+1, dist*2+1))
    probc_sum = np.zeros((dist * 2 + 1, dist * 2 + 1))
    cntc_sum = np.zeros((dist * 2 + 1, dist * 2 + 1))

    edic = {}

    probs = get_previous_hours(storm_date,lsta_date, fi.attrs['eh'], fi.attrs['refhour'])
    print('Era5 collect')

    try:
        probs_on_lsta = lsta.salem.transform(probs)
    except (RuntimeError, ValueError):
        print('Era5 on LSTA interpolation problem or ERA5 missing')
        return None
    del probs


    probs_msg = get_previous_hours_msg(storm_date,lsta_date, fi.attrs['eh'], fi.attrs['refhour'])

    probsm_on_lsta = lsta.salem.transform(probs_msg, interp='nearest')
    del probs_msg

    probs_cm = get_previous_hours_CMORPH(storm_date)   # get previous rain to storm
    try:
        probscm_on_lsta = lsta.salem.transform(probs_cm)
    except RuntimeError:
        return None
    del probs_cm
    del lsta

    mcs_hour = xr.open_dataarray(cnst.MCS_HOUR_DAILY)
    mcsimage = xr.open_dataarray(cnst.MCS_15K)
    mcsimage = mcsimage.sel(time=fi.time, lat=slice(10.2,19), lon=slice(-9.9,9.9))
    counter = 0

    labels, goodinds = ua.blob_define(mcsimage.values, -50, minmax_area=[600, 50000], max_area=None)

    for y, x in zip(pos[0], pos[1]):

        lat = fi['lat'][y]
        lon = fi['lon'][x]

        if (labels[y,x] not in goodinds) | (labels[y,x] == 0):
            print('MCS too small!!')
            continue

        if (mcsimage.values[y,x] > -60):
            print('Core too warm!!')
            continue

        # if ((mcsimage.values[y,x] < -73) | (mcsimage.values[y,x] > -67)):
        #     print('Core too cold!!')
        #     continue

        mhour = mcs_hour.sel(time=pd.datetime(int(fdate[0:4]), int(fdate[4:6]), int(fdate[6:8]), 14), lon=lon,
                             lat=lat).values
        if mhour <= 13:
            mhour += 24

        if mhour == 0:
            mhour += 24

        chour = fi['time.hour'].values

        if (chour >= 0) & (chour <= 13):
            chour += 24
        if (mhour < chour) | (np.isnan(mhour)):
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
            kernel2, cnt, cntp, vdic, probm, cntm, probcm, cntc ,plsta, lcnt = cut_kernel(xpos, ypos, lsta_da, dist, probs=probs_on_lsta, probs2=probsm_on_lsta, probs3=probscm_on_lsta, lsta_prev=prev_lsta_da)
        except TypeError:
            continue

        if np.nansum(probm[dist-50:dist+50,dist-30:dist+100])>=2:   # filter out cases with MCSs at 12 [dist-50:dist+50,dist-30:dist+100]
            print('Meteosat MCS continue')
            continue

        if np.nanmin((vdic['tciwmid'][dist-50:dist+50,dist-30:dist+100]))<=-0.4:   # 0.03 for tciw, -0.3 for w [dist-50:dist+50,dist-30:dist+100] ,filter out cases with MCSs at 12
            print('ERA MCS continue')
            continue

        # if np.nansum(probcm[100:200,:])>=100:   # filter out cases with rainfall to the south
        #     print('Southward rainfall, continue')
        #     continue

        kernel2_sum = np.nansum(np.stack([kernel2_sum, kernel2]), axis=0)
        cntp_sum = np.nansum(np.stack([cntp_sum, cntp]), axis=0)
        cnt_sum = np.nansum(np.stack([cnt_sum, cnt]), axis=0)

        cntm_sum = np.nansum(np.stack([cntm_sum, cntm]), axis=0)
        probm_sum = np.nansum(np.stack([probm_sum, probm]), axis=0)

        probc_sum = np.nansum(np.stack([probc_sum, probcm]), axis=0)
        cntc_sum = np.nansum(np.stack([cntc_sum, cntc]), axis=0)

        kernelp_sum = np.nansum(np.stack([kernelp_sum, plsta]), axis=0)
        lcnt_sum = np.nansum(np.stack([lcnt_sum, lcnt]), axis=0)

        for ks in vdic.keys():
            if ks in edic:
                edic[ks] = np.nansum(np.stack([edic[ks], vdic[ks]]), axis=0)
            else:
                edic[ks] = vdic[ks]
        counter += 1
        print('Saved core ', counter)

    del mcs_hour

    if np.sum(cnt_sum)==0:
        return None

    outlist = [kernel2_sum, cnt_sum, cntp_sum, cntm_sum, probm_sum, cntc_sum,probc_sum,kernelp_sum,lcnt_sum]
    outnames = ['lsta',  'cnt', 'cntp', 'cntm', 'probmsg', 'cntc', 'probc','plsta', 'plcnt']
    for ek in edic.keys():
        outnames.append(ek)
        outlist.append(edic[ek])

    print('Returning')
    #lsta.close()

    return outlist, outnames


def plot_eh():
    for eh in range(-38,1,3):
        plot_doug_all(17,eh)


def plot_noon():
    afternoon = list(range(15,24))
    night = list(range(0,8))
    all = afternoon + night

    hlist = []
    for hh in all:
        if hh >= 15:
            hlist.append((hh,12-hh))
        else:
            hlist.append((hh, 12-(hh+24)))

    for h in hlist:
        plot_doug_all(h[0],h[1])



def plot_doug_timeseries():

    x = len(list(range(-38, 4, 3)))
    y = 401

    # outticks = list(range(-30, 1, 5))
    # ranges = np.arange(-30,1,3)
    #
    # outticks = [12,17,22,3,8,13,18]
    #outticks = [1, 6, 11, 16, 21, 2, 7, 12, 17]
    ranges = np.arange(-38, 4, 3)

    h=17

    outdic = {}
    dummyfile = glob.glob(cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/ERA5_composite_cores_*_small_cores.p")
    dummy = pkl.load(open(dummyfile[0], "rb"))

    for k in dummy.keys():
        outdic[k] = np.zeros((y, x))

    for ids, eh in enumerate(range(-38,4,3)):

        dic = {}

        def coll(dic, h, eh, year):
            print(h)
            core = pkl.load(open(
                cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/ERA5_composite_cores_LSTA_500w04_15k_"+str(eh) + "UTCERA"+str(h).zfill(2)+'_'+str(year)+"_small_cores.p", "rb"))
            for id, k in enumerate(core.keys()):
                try:
                    dic[k] = dic[k] + core[k]
                except KeyError:
                    dic[k] = core[k]


        for y in range(2006,2011):
            coll(dic, h, eh, y)

        for k in dic.keys():
            outdic[k][:,ids] = dic[k][:,190:211].mean(axis=1)

    print(outdic.keys())
    f = plt.figure(figsize=(15,10))
    ax = f.add_subplot(321)

    plt.contourf(ranges,np.arange(401), (outdic['q'])*1000/ outdic['cntp'],levels=np.linspace(-0.5,0.5,14), cmap='RdBu', extend='both')
    plt.colorbar(label=r'g kg$^{-1}$')
    contours = plt.contour(ranges, np.arange(401), (outdic['v650']) / outdic['cntp'], levels=np.linspace(-0.8,0.8,11), cmap='RdBu')
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    ax.set_yticklabels(np.array((np.linspace(0, 200 * 2, 9) - 200) * 3, dtype=int))
    #ax.set_xticklabels(outticks)

    plt.title('Shading: q-anomaly, contours: 650hpa v-wind anomaly')
    plt.xlabel('Hour relative to convective core at 1700UTC')
    plt.ylabel('North-South distance from core (km)')

    plt.hlines(200,xmin=ranges[0], xmax=ranges[-1], linestyle='dashed')


    ax = f.add_subplot(324)

    plt.contourf(ranges,np.arange(401), (outdic['t'] / outdic['cntp']), extend='both',levels=[ -0.8,-0.7, -0.6,-0.5,-0.4,-0.2,-0.1, 0.1,0.2,0.4,0.5,0.6, 0.7, 0.8], cmap='RdBu_r')
    plt.colorbar(label=r'K')
    contours = plt.contour(ranges, np.arange(401), (outdic['u650_orig'] / outdic['cntp']), extend='both', cmap='viridis', levels=np.arange(-15,-8,0.5)) #levels=np.arange(-15,-8,0.5)) #

    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    ax.set_yticklabels(np.array((np.linspace(0, 200 * 2, 9) - 200) * 3, dtype=int))
    #ax.set_xticklabels(outticks)

    plt.title('Shading: t-anomaly, contours: 650hPa u-wind')
    plt.hlines(200,xmin=ranges[0], xmax=ranges[-1], linestyle='dashed')
    plt.xlabel('Hour of day')
    plt.ylabel('North-South distance from core (km)')


    ax = f.add_subplot(323)
    plt.contourf(ranges,np.arange(401), (outdic['t'])/ outdic['cntp'], extend='both',  cmap='RdBu_r',levels=[ -0.8,-0.7, -0.6,-0.5,-0.4,-0.2,-0.1, 0.1,0.2,0.4,0.5,0.6, 0.7, 0.8])
    plt.colorbar(label=r'K')
    contours = plt.contour(ranges, np.arange(401), (outdic['u650'] / outdic['cntp']), extend='both', cmap='RdBu', levels=np.linspace(-0.7,0.7,9))
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    ax.set_yticklabels(np.array((np.linspace(0, 200 * 2, 9) - 200) * 3, dtype=int))
    #ax.set_xticklabels(outticks)

    plt.title('Shading:t-anomaly, contours: 650hpa u-wind anomaly')
    plt.hlines(200,xmin=ranges[0], xmax=ranges[-1], linestyle='dashed')
    plt.xlabel('Hour of day')
    plt.ylabel('North-South distance from core (km)')



    ax = f.add_subplot(322)
    plt.contourf(ranges,np.arange(401), (outdic['q'])*1000/ outdic['cntp'],levels=np.linspace(-0.5,0.5,14), cmap='RdBu', extend='both')
    plt.colorbar(label=r'g kg$^{-1}$')
    contours = plt.contour(ranges, np.arange(401), (outdic['u925'] / outdic['cntp']), extend='both', cmap='RdBu', levels=np.linspace(-0.7,0.7,9))
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    ax.set_yticklabels(np.array((np.linspace(0, 200 * 2, 9) - 200) * 3, dtype=int))
    #ax.set_xticklabels(outticks)
    plt.title('Shading: q-anomaly, contours: 925hpa u-wind anomaly')
    plt.hlines(200,xmin=ranges[0], xmax=ranges[-1], linestyle='dashed')
    plt.xlabel('Hour of day')
    plt.ylabel('North-South distance from core (km)')

    ax = f.add_subplot(325)
    plt.contourf(ranges,np.arange(401), (outdic['theta_e'] / outdic['cntp']), extend='both', levels=np.linspace(-2,2,12), cmap='PuOr_r')
    plt.colorbar(label=r'K')
    contours = plt.contour(ranges, np.arange(401), (outdic['v925'] / outdic['cntp']), extend='both', cmap='RdBu', levels=np.linspace(-0.7,0.7,9))
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    ax.set_yticklabels(np.array((np.linspace(0, 200 * 2, 9) - 200) * 3, dtype=int))
    #ax.set_xticklabels(outticks)

    plt.title(r'Shading:$\Delta \theta_{e} anomaly$, contours: 925hpa v-wind anomaly')
    plt.hlines(200,xmin=ranges[0], xmax=ranges[-1], linestyle='dashed')
    plt.xlabel('Hour of day')
    plt.ylabel('North-South distance from core (km)')

    ax = f.add_subplot(326)
    plt.contourf(ranges,np.arange(401), (outdic['div'])/ outdic['cntp']*100, extend='both',  cmap='RdBu', levels=np.linspace(-0.5,0.5,10))
    plt.colorbar(label=r'10$^{-2}$ s$^{-1}$')
    contours = plt.contour(ranges, np.arange(401),(outdic['v925_orig'] / outdic['cntp']), extend='both',levels=np.linspace(-3,3,11), cmap='RdBu')
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    ax.set_yticklabels(np.array((np.linspace(0, 200 * 2, 9) - 200) * 3, dtype=int))
    #ax.set_xticklabels(outticks)
    plt.title(r'Shading:Divergence, contours: 925hpa v-wind')
    plt.hlines(200,xmin=ranges[0], xmax=ranges[-1], linestyle='dashed')
    plt.xlabel('Hour of day')


    plt.tight_layout()
    plt.show()
    plt.savefig(cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/plots/ERA5_"+str(h).zfill(2)+'_timeseries_short.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png)
    #plt.close()



def plot_doug_all_diff(h, eh):

    dic = {}
    dic2 = {}

    name = "ERA5_composite_cores_LSTA_500w04_15k_LSTAcp90"#"ERA5_composite_cores_AMSRE_w1_15k_minusMean"
    name2 = "ERA5_composite_cores_LSTA_500w04_15k_LSTAcp10"

    def coll(dic, h, eh, year,name):
        print(h)
        core = pkl.load(open(
            cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/"+name+str(eh) + "UTCERA"+str(h).zfill(2)+'_'+str(year)+"_small_cores.p", "rb"))
        for id, k in enumerate(core.keys()):
            try:
                dic[k] = dic[k] + core[k]
            except KeyError:
                dic[k] = core[k]


    for y in range(2006,2011):
        coll(dic, h, eh, y,name)

    for y in range(2006,2011):
        coll(dic2, h, eh, y,name2)

    extent = (dic['lsta'].shape[1]-1)/2
    xlen = dic['lsta'].shape[1]
    ylen = dic['lsta'].shape[0]

    xv, yv = np.meshgrid(np.arange(ylen), np.arange(xlen))
    st=30
    xquiv = xv[4::st, 4::st]
    yquiv = yv[4::st, 4::st]

    u = (dic['u925']/ dic['cntp'])[4::st, 4::st]
    v = (dic['v925']/ dic['cntp'])[4::st, 4::st]

    u600 = (dic['u650_orig']/ dic['cntp'])[4::st, 4::st]
    v600 = (dic['v650_orig']/ dic['cntp'])[4::st, 4::st]

    u_orig = (dic['v925_orig']/ dic['cntp'])[4::st, 4::st]
    v_orig = (dic['v925_orig']/ dic['cntp'])[4::st, 4::st]


    f = plt.figure(figsize=(15,7))
    ax = f.add_subplot(231)


    plt.contourf((dic['plsta'] / dic['plcnt'])-((dic2['plsta'] / dic2['plcnt'])), cmap='RdBu_r', extend='both', levels=np.linspace(-5,5,12)) #-(rkernel2_sum / rcnt_sum)
    plt.plot(extent, extent, 'bo')
    plt.colorbar(label='K', format='%1.2f')
    plt.text(0.02,0.08, 'ITD 0-line', color='turquoise', fontsize=12, transform=ax.transAxes)
    plt.text(0.02, 0.03, 'ITD anomaly 0-line', color='k', fontsize=12, transform=ax.transAxes)

    # plt.annotate('ITD 0-line', xy=(0.04, 0.1), xytext=(0, 4), size=15, color='turquoise', xycoords=('figure fraction', 'figure fraction'))
    #              #             textcoords='offset points')   #transform=ax.transAxes
    #pdb.set_trace()

    contours = plt.contour((dic['t'] / dic['cntp'])-(dic2['t'] / dic2['cntp']), extend='both', cmap='PuOr_r', linewidths=2, levels=np.linspace(-1,1,8)) # #, levels=np.arange(1,5, 0.5)
    plt.clabel(contours, inline=True, fontsize=11, fmt='%1.1f')
    #contours2 = plt.contour((dic['v925']) / dic['cntp'], extend='both', cmap='RdBu', levels=np.linspace(-1, 1,9))  # , levels=np.linspace(-1,1,10)#(dic['probc']/ dic['cntc'])*100, extend='both', levels=np.arange(15,70,12), cmap='jet') # #, levels=np.arange(1,5, 0.5)
    #plt.clabel(contours2, inline=True, fontsize=11, fmt='%1.0f')

    contours = plt.contour((dic['v925']/ dic['cntc'])-(dic2['v925']/ dic2['cntc']), extend='both',colors='k', linewidths=4, levels=[-50,0,50])  #np.arange(-15,-10,0.5)
    #plt.clabel(contours, inline=True, fontsize=9, fmt='%1.1f')

    contours = plt.contour((dic['v925_orig']/ dic['cntc'])-(dic2['v925_orig']/ dic2['cntc']), extend='both',colors='turquoise', linewidths=4, levels=[-50,0,50])  #np.arange(-15,-10,0.5)
    #plt.clabel(contours, inline=True, fontsize=9, fmt='%1.1f')

    plt.plot(extent, extent, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6 , dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.title(str(h).zfill(2)+'00UTC | '+str(np.max(dic['cnt']))+' cores, LSTA day-1, ERA5$_{noon}$ 925hPa T anomaly', fontsize=9)


    ax1 = f.add_subplot(232)
    plt.contourf(((dic['lsta'])/ dic['cnt'])-((dic2['lsta'])/ dic2['cnt']), extend='both',  cmap='RdBu_r', levels=np.linspace(-5,5,12)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='K', format='%1.2f')
    contours = plt.contour((dic['probc']/ dic['cntc'])*100 - (dic2['probc']/ dic2['cntc'])*100, extend='both', cmap='jet', linewidths=2)  #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.1f')



    #qu = ax1.quiver(xquiv, yquiv, u, v, scale=50)
    plt.plot(extent, extent, 'bo')
    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('LSTA day0, Contours: CMORPH rainP>5mm [6am|day-1 to 10am|day0]', fontsize=9)

    ax1 = f.add_subplot(233)
    plt.contourf(((dic['q'])*1000/ dic['cntp'])-((dic2['q'])*1000/ dic2['cntp']), extend='both',  cmap='RdBu', levels=np.linspace(-0.8,0.8,10)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='g kg-1')
    contours = plt.contour((dic['shear'] / dic['cntp'])-(dic2['shear'] / dic2['cntp']), extend='both', cmap='viridis_r', linewidths=2) #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    plt.plot(extent, extent, 'bo')
    #qu = ax1.quiver(xquiv, yquiv, u, v, scale=50)
    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('925hPa q anomaly, contours: 650hPa-925hPa zonal shear', fontsize=9)

    ax1 = f.add_subplot(234)
 #   plt.contourf(((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r', vmin=-1.5, vmax=1.5) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.contourf((((dic['div'])/ dic['cntp'])-((dic2['div'])/ dic2['cntp']))*100, extend='both',  cmap='PuOr', levels=np.linspace(-0.7, 0.7, 10)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='10$^{-2}$ s$^{-1}$', format='%1.3f')
    plt.plot(extent, extent, 'bo')
    # contours = plt.contour((dic['v650_orig'] / dic['cntp']), extend='both', cmap='viridis') #np.arange(-15,-10,0.5)
    # plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    # qu = ax1.quiver(xquiv, yquiv, u, v, scale=15)
    #qk = plt.quiverkey(qu, 0.2, 0.02,1, '1 m s$^{-1}$',
    #                   labelpos='E', coordinates='figure')

    ax1.streamplot(xv, yv, (dic['u925']/ dic['cntp'])-(dic2['u925']/ dic2['cntp']), (dic['v925']/ dic['cntp'])-(dic2['v925']/ dic2['cntp']), density=[0.5, 1])

    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('925hPa divergence, vectors: 925hPa wind anomaly', fontsize=9)


    ax1 = f.add_subplot(235)
    #   plt.contourf(((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r', vmin=-1.5, vmax=1.5) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
     # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.contourf((dic['tciwmid'] / dic['cntp'])-(dic2['tciwmid'] / dic2['cntp']), extend='both', cmap='PuOr', levels=np.linspace(-0.02, 0.02, 10))
    plt.colorbar(label=r'Pa s$^{-1}$', format='%1.3f')
    plt.plot(extent, extent, 'bo')
    # contours = plt.contour(((dic['v925_orig']) / dic['cntp']) , extend='both', cmap='RdBu', levels=np.linspace(-2,2,11), linewidths=2) #np.arange(-15,-10,0.5)
    # plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    ax1.streamplot(xv, yv, dic['u925_orig']/ dic['cntp']-dic2['u925_orig']/ dic2['cntp'], dic['v925_orig']/ dic['cntp']-dic2['v925_orig']/ dic2['cntp'], density=[0.5, 1])

    # qu = ax1.quiver(xquiv, yquiv, u_orig, v_orig, scale=30)
    # qk = plt.quiverkey(qu, 0.55, 0.02, 1, '1 m s$^{-1}$',
    #                    labelpos='E', coordinates='figure')

    ax1.set_xticklabels(np.array((np.linspace(0, extent * 2, 9) - 100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent * 2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title(r'500hpa omega, vectors: 925hPa wind', fontsize=9)

    ax1 = f.add_subplot(236)
    #   plt.contourf(((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r', vmin=-1.5, vmax=1.5) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.contourf(((dic['theta_e']) / dic['cntp'])-((dic2['theta_e']) / dic2['cntp']), extend='both', cmap='RdBu', levels=np.linspace(-2.5,2.5,14)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label=r'K')
    plt.plot(extent, extent, 'bo')
    contours = plt.contour((dic['u650'] / dic['cntp'])-(dic2['u650'] / dic2['cntp']), extend='both', cmap='PuOr', levels=[-2,-1.5,-1,-0.5,-0.2,-0.1, 0,0.1,0.2,0.5,1,1.5,2], linewidths=2) #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    # qu = ax1.quiver(xquiv, yquiv, u600, v600, scale=20)
    # qk = plt.quiverkey(qu, 0.9, 0.02, 1, '1 m s$^{-1}$',
    #                    labelpos='E', coordinates='figure')

    ax1.set_xticklabels(np.array((np.linspace(0, extent * 2, 9) - 100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent * 2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title(r'$\Delta \theta_{e}$ anomaly, contours: 650hPa u-wind anomaly', fontsize=9)

    plt.tight_layout()
    plt.show()
    plt.savefig(cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/plots/"+name+str(h).zfill(2)+'_'+str(eh).zfill(2)+'_DIFF_9010.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png)
    #plt.close()



def plot_doug_all(h, eh):

    dic = {}
    dic2 = {}

    name = "ERA5_composite_cores_LSTA_500w04_15k_WET_"#"ERA5_composite_cores_AMSRE_w1_15k_minusMean"


    def coll(dic, h, eh, year):
        print(h)
        core = pkl.load(open(
            cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/"+name+str(eh) + "UTCERA"+str(h).zfill(2)+'_'+str(year)+"_small_cores.p", "rb"))
        for id, k in enumerate(core.keys()):
            try:
                dic[k] = dic[k] + core[k]
            except KeyError:
                dic[k] = core[k]


    for y in range(2006,2011):
        coll(dic, h, eh, y)

    extent = (dic['lsta'].shape[1]-1)/2
    xlen = dic['lsta'].shape[1]
    ylen = dic['lsta'].shape[0]

    xv, yv = np.meshgrid(np.arange(ylen), np.arange(xlen))
    st=30
    xquiv = xv[4::st, 4::st]
    yquiv = yv[4::st, 4::st]

    u = (dic['u925']/ dic['cntp'])[4::st, 4::st]
    v = (dic['v925']/ dic['cntp'])[4::st, 4::st]

    u600 = (dic['u650_orig']/ dic['cntp'])[4::st, 4::st]
    v600 = (dic['v650_orig']/ dic['cntp'])[4::st, 4::st]

    u_orig = (dic['v925_orig']/ dic['cntp'])[4::st, 4::st]
    v_orig = (dic['v925_orig']/ dic['cntp'])[4::st, 4::st]


    f = plt.figure(figsize=(15,7))
    ax = f.add_subplot(231)


    plt.contourf((dic['plsta'] / dic['plcnt'])-np.mean((dic['plsta'] / dic['plcnt'])), cmap='RdBu_r', levels=np.linspace(-2,2,16), extend='both') #-(rkernel2_sum / rcnt_sum)
    plt.plot(extent, extent, 'bo')
    plt.colorbar(label='K', format='%1.2f')
    plt.text(0.02,0.08, 'ITD 0-line', color='turquoise', fontsize=12, transform=ax.transAxes)
    plt.text(0.02, 0.03, 'ITD anomaly 0-line', color='k', fontsize=12, transform=ax.transAxes)

    # plt.annotate('ITD 0-line', xy=(0.04, 0.1), xytext=(0, 4), size=15, color='turquoise', xycoords=('figure fraction', 'figure fraction'))
    #              #             textcoords='offset points')   #transform=ax.transAxes
    #pdb.set_trace()

    contours = plt.contour((dic['t'] / dic['cntp']), extend='both',levels=[ -0.8,-0.7, -0.6,-0.5,-0.4,-0.2,0, 0.2,0.4,0.5,0.6, 0.7, 0.8], cmap='PuOr_r', linewidths=2) # #, levels=np.arange(1,5, 0.5)
    plt.clabel(contours, inline=True, fontsize=11, fmt='%1.1f')
    #contours2 = plt.contour((dic['v925']) / dic['cntp'], extend='both', cmap='RdBu', levels=np.linspace(-1, 1,9))  # , levels=np.linspace(-1,1,10)#(dic['probc']/ dic['cntc'])*100, extend='both', levels=np.arange(15,70,12), cmap='jet') # #, levels=np.arange(1,5, 0.5)
    #plt.clabel(contours2, inline=True, fontsize=11, fmt='%1.0f')

    contours = plt.contour((dic['v925']/ dic['cntc']), extend='both',colors='k', linewidths=4, levels=[-50,0,50])  #np.arange(-15,-10,0.5)
    #plt.clabel(contours, inline=True, fontsize=9, fmt='%1.1f')

    contours = plt.contour((dic['v925_orig']/ dic['cntc']), extend='both',colors='turquoise', linewidths=4, levels=[-50,0,50])  #np.arange(-15,-10,0.5)
    #plt.clabel(contours, inline=True, fontsize=9, fmt='%1.1f')

    plt.plot(extent, extent, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6 , dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.title(str(h).zfill(2)+'00UTC | '+str(np.max(dic['cnt']))+' cores, LSTA day-1, ERA5$_{noon}$ 925hPa T anomaly', fontsize=9)


    ax1 = f.add_subplot(232)
    plt.contourf(((dic['lsta'])/ dic['cnt'])-np.mean((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r', levels=np.linspace(-2,2,16)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='K', format='%1.2f')
    contours = plt.contour((dic['probc']/ dic['cntc'])*100, extend='both', levels=np.arange(15,70,12), cmap='jet', linewidths=2)  #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.1f')



    #qu = ax1.quiver(xquiv, yquiv, u, v, scale=50)
    plt.plot(extent, extent, 'bo')
    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('LSTA day0, Contours: CMORPH rainP>5mm [6am|day-1 to 10am|day0]', fontsize=9)

    ax1 = f.add_subplot(233)
    plt.contourf(((dic['q'])*1000/ dic['cntp']), extend='both',  cmap='RdBu',levels=np.linspace(-0.9,0.9,16)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='g kg-1')
    contours = plt.contour((dic['shear'] / dic['cntp']), extend='both',levels=np.arange(-17,-12,0.5), cmap='viridis_r', linewidths=2) #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    plt.plot(extent, extent, 'bo')
    #qu = ax1.quiver(xquiv, yquiv, u, v, scale=50)
    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('925hPa q anomaly, contours: 650hPa-925hPa zonal shear', fontsize=9)

    ax1 = f.add_subplot(234)
 #   plt.contourf(((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r', vmin=-1.5, vmax=1.5) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.contourf(((dic['div'])/ dic['cntp'])*100, extend='both',  cmap='PuOr', levels=np.linspace(-0.7,0.7,10)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='10$^{-2}$ s$^{-1}$', format='%1.3f')
    plt.plot(extent, extent, 'bo')
    # contours = plt.contour((dic['v650_orig'] / dic['cntp']), extend='both', cmap='viridis') #np.arange(-15,-10,0.5)
    # plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    # qu = ax1.quiver(xquiv, yquiv, u, v, scale=15)
    #qk = plt.quiverkey(qu, 0.2, 0.02,1, '1 m s$^{-1}$',
    #                   labelpos='E', coordinates='figure')

    ax1.streamplot(xv, yv, (dic['u925']/ dic['cntp']), (dic['v925']/ dic['cntp']), density=[0.5, 1])

    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('925hPa divergence, vectors: 925hPa wind anomaly', fontsize=9)


    ax1 = f.add_subplot(235)
    #   plt.contourf(((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r', vmin=-1.5, vmax=1.5) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
     # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.contourf((dic['tciwmid'] / dic['cntp']), extend='both', cmap='PuOr', levels=np.linspace(-0.02, 0.02, 10))
    plt.colorbar(label=r'Pa s$^{-1}$', format='%1.3f')
    plt.plot(extent, extent, 'bo')
    # contours = plt.contour(((dic['v925_orig']) / dic['cntp']) , extend='both', cmap='RdBu', levels=np.linspace(-2,2,11), linewidths=2) #np.arange(-15,-10,0.5)
    # plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    ax1.streamplot(xv, yv, dic['u925_orig']/ dic['cntp'], dic['v925_orig']/ dic['cntp'], density=[0.5, 1])

    # qu = ax1.quiver(xquiv, yquiv, u_orig, v_orig, scale=30)
    # qk = plt.quiverkey(qu, 0.55, 0.02, 1, '1 m s$^{-1}$',
    #                    labelpos='E', coordinates='figure')

    ax1.set_xticklabels(np.array((np.linspace(0, extent * 2, 9) - 100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent * 2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title(r'500hpa omega, vectors: 925hPa wind', fontsize=9)

    ax1 = f.add_subplot(236)
    #   plt.contourf(((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r', vmin=-1.5, vmax=1.5) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.contourf(((dic['theta_e']) / dic['cntp']), extend='both', cmap='RdBu', levels=np.linspace(-3,3,18)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label=r'K')
    plt.plot(extent, extent, 'bo')
    contours = plt.contour((dic['u650'] / dic['cntp']), extend='both', cmap='PuOr', levels=[-2,-1.5,-1,-0.5,-0.2,-0.1, 0,0.1,0.2,0.5,1,1.5,2], linewidths=2) #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    # qu = ax1.quiver(xquiv, yquiv, u600, v600, scale=20)
    # qk = plt.quiverkey(qu, 0.9, 0.02, 1, '1 m s$^{-1}$',
    #                    labelpos='E', coordinates='figure')

    ax1.set_xticklabels(np.array((np.linspace(0, extent * 2, 9) - 100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent * 2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title(r'$\Delta \theta_{e}$ anomaly, contours: 650hPa u-wind anomaly', fontsize=9)

    plt.tight_layout()
    plt.show()
    plt.savefig(cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/plots/"+name+str(h).zfill(2)+'_'+str(eh).zfill(2)+'.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png)
    #plt.close()


def plot_doug_CLIM(h, eh):

    dic = {}
    dic2 = {}

    name = "ERA5_composite_cores_LSTA_500w04_15k_-80_clim"#"ERA5_composite_cores_AMSRE_w1_15k_minusMean"


    def coll(dic, h, eh, year):
        print(h)
        core = pkl.load(open(
            cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/"+name+str(eh) + "UTCERA"+str(h).zfill(2)+'_'+str(year)+"_small_cores.p", "rb"))
        for id, k in enumerate(core.keys()):
            try:
                dic[k] = dic[k] + core[k]
            except KeyError:
                dic[k] = core[k]


    for y in range(2006,2011):
        coll(dic, h, eh, y)

    extent = (dic['lsta'].shape[1]-1)/2
    xlen = dic['lsta'].shape[1]
    ylen = dic['lsta'].shape[0]

    xv, yv = np.meshgrid(np.arange(ylen), np.arange(xlen))
    st=30
    xquiv = xv[4::st, 4::st]
    yquiv = yv[4::st, 4::st]

    u = (dic['u925']/ dic['cntp'])[4::st, 4::st]
    v = (dic['v925']/ dic['cntp'])[4::st, 4::st]

    u600 = (dic['u650_orig']/ dic['cntp'])[4::st, 4::st]
    v600 = (dic['v650_orig']/ dic['cntp'])[4::st, 4::st]

    u_orig = (dic['v925_orig']/ dic['cntp'])[4::st, 4::st]
    v_orig = (dic['v925_orig']/ dic['cntp'])[4::st, 4::st]


    f = plt.figure(figsize=(15,7))
    ax = f.add_subplot(231)


    plt.contourf((dic['plsta'] / dic['plcnt'])-np.mean((dic['plsta'] / dic['plcnt'])), cmap='RdBu_r', levels=np.linspace(-1.5,1.5,12), extend='both') #-(rkernel2_sum / rcnt_sum)
    plt.plot(extent, extent, 'bo')
    plt.colorbar(label='K', format='%1.2f')
    plt.text(0.02,0.08, 'ITD 0-line', color='turquoise', fontsize=12, transform=ax.transAxes)
    plt.text(0.02, 0.03, 'ITD anomaly 0-line', color='k', fontsize=12, transform=ax.transAxes)

    # plt.annotate('ITD 0-line', xy=(0.04, 0.1), xytext=(0, 4), size=15, color='turquoise', xycoords=('figure fraction', 'figure fraction'))
    #              #             textcoords='offset points')   #transform=ax.transAxes
    #pdb.set_trace()

    contours = plt.contour((dic['t'] / dic['cntp']), extend='both',levels=[ -0.8,-0.7, -0.6,-0.5,-0.4,-0.2,0, 0.2,0.4,0.5,0.6, 0.7, 0.8], cmap='PuOr_r', linewidths=2) # #, levels=np.arange(1,5, 0.5)
    plt.clabel(contours, inline=True, fontsize=11, fmt='%1.1f')
    #contours2 = plt.contour((dic['v925']) / dic['cntp'], extend='both', cmap='RdBu', levels=np.linspace(-1, 1,9))  # , levels=np.linspace(-1,1,10)#(dic['probc']/ dic['cntc'])*100, extend='both', levels=np.arange(15,70,12), cmap='jet') # #, levels=np.arange(1,5, 0.5)
    #plt.clabel(contours2, inline=True, fontsize=11, fmt='%1.0f')

    contours = plt.contour((dic['v925']/ dic['cntc']), extend='both',colors='k', linewidths=4, levels=[-50,0,50])  #np.arange(-15,-10,0.5)
    #plt.clabel(contours, inline=True, fontsize=9, fmt='%1.1f')

    contours = plt.contour((dic['v925_orig']/ dic['cntc']), extend='both',colors='turquoise', linewidths=4, levels=[-50,0,50])  #np.arange(-15,-10,0.5)
    #plt.clabel(contours, inline=True, fontsize=9, fmt='%1.1f')

    plt.plot(extent, extent, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6 , dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.title(str(h).zfill(2)+'00UTC | '+str(np.max(dic['cnt']))+' cores, LSTA day-1, ERA5$_{noon}$ 925hPa T anomaly', fontsize=9)


    ax1 = f.add_subplot(232)
    plt.contourf(((dic['lsta'])/ dic['cnt'])-np.mean((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r', levels=np.linspace(-1.5,1.5,12)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='K', format='%1.2f')
    contours = plt.contour((dic['probc']/ dic['cntc'])*100, extend='both', levels=np.arange(15,70,12), cmap='jet', linewidths=2)  #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.1f')



    #qu = ax1.quiver(xquiv, yquiv, u, v, scale=50)
    plt.plot(extent, extent, 'bo')
    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('LSTA day0, Contours: CMORPH rainP>5mm [6am|day-1 to 10am|day0]', fontsize=9)

    ax1 = f.add_subplot(233)
    plt.contourf(((dic['q'])*1000/ dic['cntp']), extend='both',  cmap='jet', levels=np.arange(11,17,0.5)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='g kg-1')
    contours = plt.contour((dic['shear'] / dic['cntp']), extend='both',levels=np.arange(-17,-12,0.5), cmap='viridis_r', linewidths=2) #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    plt.plot(extent, extent, 'bo')
    #qu = ax1.quiver(xquiv, yquiv, u, v, scale=50)
    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('925hPa q anomaly, contours: 650hPa-925hPa zonal shear', fontsize=9)

    ax1 = f.add_subplot(234)
 #   plt.contourf(((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r', vmin=-1.5, vmax=1.5) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.contourf(((dic['div'])/ dic['cntp'])*100, extend='both',  cmap='PuOr', levels=np.linspace(-0.7,0.7,10)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='10$^{-2}$ s$^{-1}$', format='%1.3f')
    plt.plot(extent, extent, 'bo')
    # contours = plt.contour((dic['v650_orig'] / dic['cntp']), extend='both', cmap='viridis') #np.arange(-15,-10,0.5)
    # plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    # qu = ax1.quiver(xquiv, yquiv, u, v, scale=15)
    #qk = plt.quiverkey(qu, 0.2, 0.02,1, '1 m s$^{-1}$',
    #                   labelpos='E', coordinates='figure')

    ax1.streamplot(xv, yv, (dic['u925']/ dic['cntp']), (dic['v925']/ dic['cntp']), density=[0.5, 1])

    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('925hPa divergence, vectors: 925hPa wind anomaly', fontsize=9)


    ax1 = f.add_subplot(235)
    #   plt.contourf(((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r', vmin=-1.5, vmax=1.5) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
     # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.contourf((dic['tciwmid'] / dic['cntp']), extend='both', cmap='PuOr', levels=np.linspace(-0.02, 0.02, 10))
    plt.colorbar(label=r'Pa s$^{-1}$', format='%1.3f')
    plt.plot(extent, extent, 'bo')
    # contours = plt.contour(((dic['v925_orig']) / dic['cntp']) , extend='both', cmap='RdBu', levels=np.linspace(-2,2,11), linewidths=2) #np.arange(-15,-10,0.5)
    # plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    ax1.streamplot(xv, yv, dic['u925_orig']/ dic['cntp'], dic['v925_orig']/ dic['cntp'], density=[0.5, 1])

    # qu = ax1.quiver(xquiv, yquiv, u_orig, v_orig, scale=30)
    # qk = plt.quiverkey(qu, 0.55, 0.02, 1, '1 m s$^{-1}$',
    #                    labelpos='E', coordinates='figure')

    ax1.set_xticklabels(np.array((np.linspace(0, extent * 2, 9) - 100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent * 2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title(r'500hpa omega, vectors: 925hPa wind', fontsize=9)

    ax1 = f.add_subplot(236)
    #   plt.contourf(((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r', vmin=-1.5, vmax=1.5) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.contourf(((dic['theta_e']) / dic['cntp']), extend='both', cmap='RdBu', levels=np.arange(5,13,0.5)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label=r'K')
    plt.plot(extent, extent, 'bo')
    contours = plt.contour((dic['u650'] / dic['cntp']), extend='both', cmap='PuOr', levels=[-2,-1.5,-1,-0.5,-0.2,-0.1, 0,0.1,0.2,0.5,1,1.5,2], linewidths=2) #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    # qu = ax1.quiver(xquiv, yquiv, u600, v600, scale=20)
    # qk = plt.quiverkey(qu, 0.9, 0.02, 1, '1 m s$^{-1}$',
    #                    labelpos='E', coordinates='figure')

    ax1.set_xticklabels(np.array((np.linspace(0, extent * 2, 9) - 100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent * 2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title(r'$\Delta \theta_{e}$ anomaly, contours: 650hPa u-wind anomaly', fontsize=9)

    plt.tight_layout()
    plt.show()
    plt.savefig(cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/plots/"+name+str(h).zfill(2)+'_'+str(eh).zfill(2)+'_CLIM.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png)
    #plt.close()


def plot_doug_paper(h, eh):

    dic = {}
    name = "ERA5_composite_cores_LSTA_500w04_15k_"


    def coll(dic, h, eh, year):
        print(h)
        core = pkl.load(open(
            cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/"+name+str(eh) + "UTCERA"+str(h).zfill(2)+'_'+str(year)+"_small_cores.p", "rb"))
        for id, k in enumerate(core.keys()):
            try:
                dic[k] = dic[k] + core[k]
            except KeyError:
                dic[k] = core[k]


    for y in range(2006,2011):
        coll(dic, h, eh, y)

    extent = (dic['lsta'].shape[1]-1)/2
    xlen = dic['lsta'].shape[1]
    ylen = dic['lsta'].shape[0]

    xv, yv = np.meshgrid(np.arange(ylen), np.arange(xlen))
    st=30
    xquiv = xv[4::st, 4::st]
    yquiv = yv[4::st, 4::st]

    u = (dic['u925']/ dic['cntp'])[4::st, 4::st]
    v = (dic['v925']/ dic['cntp'])[4::st, 4::st]

    u600 = (dic['u650_orig']/ dic['cntp'])[4::st, 4::st]
    v600 = (dic['v650_orig']/ dic['cntp'])[4::st, 4::st]

    u_orig = (dic['v925_orig']/ dic['cntp'])[4::st, 4::st]
    v_orig = (dic['v925_orig']/ dic['cntp'])[4::st, 4::st]


    f = plt.figure(figsize=(15,7))


    f = plt.figure(figsize=(10,8))
    ax = f.add_subplot(221)

    plt.contourf((dic['plsta'] / dic['plcnt']), cmap='RdBu_r', levels=[-1.5,-1.25,-1 ,-0.75,-0.5,-0.25,0.25,0.5,0.75,1,1.25, 1.5], extend='both') #-(rkernel2_sum / rcnt_sum)
    plt.plot(extent, extent, 'bo')
    plt.colorbar(label='K')
    #pdb.set_trace()

    contours = plt.contour((dic['probc']/ dic['cntc'])*100, extend='both', levels=np.arange(15,70,12), cmap='jet') # #, levels=np.arange(1,5, 0.5)
    plt.clabel(contours, inline=True, fontsize=11, fmt='%1.0f')
    plt.plot(extent, extent, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6 , dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.title('1800UTC | '+str(np.max(dic['cnt']))+' cores, LSTA', fontsize=9)


    ax1 = f.add_subplot(222)
    plt.contourf(((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r',levels=[-1.5,-1.25,-1 ,-0.75,-0.5,-0.25,0.25,0.5,0.75,1,1.25, 1.5]) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='K')
    contours = plt.contour((dic['t'] / dic['cntp']), extend='both',levels=[ -0.8,-0.7, -0.6,-0.5,-0.4,-0.2,-0.1, 0.1,0.2,0.4,0.5,0.6, 0.7, 0.8], cmap='PuOr_r') #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.1f')
    #qu = ax1.quiver(xquiv, yquiv, u, v, scale=50)
    plt.plot(extent, extent, 'bo')
    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('Shading: MSG LSTA, Contours: ERA5 925hPa T', fontsize=9)

    ax1 = f.add_subplot(223)
    plt.contourf(((dic['q'])*1000/ dic['cntp']), extend='both',  cmap='RdBu',levels=np.arange(-1,1.1,0.1)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='g kg-1')
    contours = plt.contour((dic['shear'] / dic['cntp']), extend='both',levels=np.arange(-17,-12,0.25), cmap='viridis') #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    plt.plot(extent, extent, 'bo')
    #qu = ax1.quiver(xquiv, yquiv, u, v, scale=50)
    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('Shading: 925hPa q anomaly, Contours: 650hPa-925hPa wind shear ', fontsize=9)

    ax1 = f.add_subplot(224)
 #   plt.contourf(((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r', vmin=-1.5, vmax=1.5) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.contourf(((dic['div'])/ dic['cntp'])*100, extend='both',  cmap='RdBu', levels=np.linspace(-0.5,0.5,10)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
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
    plt.title('Shading: Divergence, vectors: 925hPa wind anomaly', fontsize=9)


    plt.tight_layout()
    plt.show()
    plt.savefig(cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/plots/ERA5_"+str(h).zfill(2)+'_'+str(eh).zfill(2)+'.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png)
    plt.close()

def plot_doug_small(h, eh):

    dic = {}
    dic2 = {}

    name = "ERA5_composite_cores_LSTA_500w04_15k_"#"ERA5_composite_cores_AMSRE_w1_15k_minusMean"


    def coll(dic, h, eh, year):
        print(h)
        core = pkl.load(open(
            cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/"+name+str(eh) + "UTCERA"+str(h).zfill(2)+'_'+str(year)+"_small_cores.p", "rb"))
        for id, k in enumerate(core.keys()):
            try:
                dic[k] = dic[k] + core[k]
            except KeyError:
                dic[k] = core[k]


    for y in range(2006,2011):
        coll(dic, h, eh, y)

    extent = (dic['lsta'].shape[1]-1)/2
    xlen = dic['lsta'].shape[1]
    ylen = dic['lsta'].shape[0]

    xv, yv = np.meshgrid(np.arange(ylen), np.arange(xlen))
    st=30
    xquiv = xv[4::st, 4::st]
    yquiv = yv[4::st, 4::st]

    u = (dic['u925']/ dic['cntp'])[4::st, 4::st]
    v = (dic['v925']/ dic['cntp'])[4::st, 4::st]

    u600 = (dic['u650_orig']/ dic['cntp'])[4::st, 4::st]
    v600 = (dic['v650_orig']/ dic['cntp'])[4::st, 4::st]

    u_orig = (dic['v925_orig']/ dic['cntp'])[4::st, 4::st]
    v_orig = (dic['v925_orig']/ dic['cntp'])[4::st, 4::st]

    f = plt.figure(figsize=(10, 4))
    ax = f.add_subplot(121)

    plt.contourf((dic['lsta'] / dic['cnt']), cmap='RdBu_r', levels = np.linspace(-1.5, 1.5, 12), extend = 'both')
    # plt.plot(extent, extent, 'bo')
    plt.colorbar(label='K')
    # pdb.set_trace()

    contours = plt.contour(((dic['div'])/ dic['cntp'])*100, extend='both',  cmap='PuOr_r', levels=np.linspace(-0.7,0.7,10)) # #, levels=np.arange(1,5, 0.5)
    qu = ax.quiver(xquiv, yquiv, u, v, scale=15)
    qk = plt.quiverkey(qu, 0.9, 0.02,1, '1 m s$^{-1}$',
                       labelpos='E', coordinates='figure')


    plt.clabel(contours, inline=True, fontsize=11, fmt='%1.0f')
    plt.plot(extent, extent, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, extent * 2, 9) - 100) * 6, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent * 2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.title('23-01UTC | ' + str(np.max(dic['cnt'])) + ' cores, LSTA & 06-06UTC antecedent rain', fontsize=9)


    ax1 = f.add_subplot(122)
    plt.contourf(((dic['q'])*1000/ dic['cntp']), extend='both',  cmap='RdBu',levels=np.arange(-1,1.1,0.1)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='g kg-1')
    contours = plt.contour((dic['shear'] / dic['cntp']), extend='both',levels=np.arange(-15.5,-12,0.5), cmap='viridis_r') #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    plt.plot(extent, extent, 'bo')
    #qu = ax1.quiver(xquiv, yquiv, u, v, scale=50)
    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('Shading: 950hPa q anomaly, Contours: 600hPa-925hPa wind shear ', fontsize=9)

    plt.tight_layout()
    plt.show()
    #plt.savefig('/users/global/cornkle/figs/LSTA-bullshit/corrected_LSTA/system_scale/doug/large_scale/'+str(hour).zfill(2)+'_JUN.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png')
    # plt.close()
