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
    all = afternoon + night

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

    all = np.arange(-41,1,3)


    #hlist = []
    # for hh in all:
    #     if hh >= 14:
    #         hlist.append((hh,12-hh))
    #     else:
    #         hlist.append((hh, 12-(hh+24)))

    for l in all:
        print('Doing '+str(l))
        composite(20, l)


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
        # res = []
        # for m in msgo[0:100]:
        #     out = file_loop(m)
        #     res.append(out)


        pkl.dump(dic, open(cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/ERA5_composite_cores_AMSRE_500w04_15k_p90"+str(eh) + "UTCERA"+str(hour).zfill(2)+'_'+str(year)+"_small_cores.p", "wb"))
        del dic
        print('Dumped file')



def cut_kernel(xpos, ypos, arrlist, dist, era=False, msg=False, cmorph=False):

    smoutlist  = []
    scntlist = []
    wetflag = 0
    for ids, arr in enumerate(arrlist):
        kernel = ua.cut_kernel(arr,xpos, ypos,dist)
        kernel = kernel - np.nanmean(kernel)
        if kernel.shape != (dist*2+1, dist*2+1):
            print('Kernels shape wrong!')
            ipdb.set_trace()
        if ids == 0:
            if (np.sum(np.isfinite(kernel)) < 2):
                return

        if (ids == 0) |(ids == 1):

            #ycirc100e, xcirc100e = ua.draw_circle(dist + 100, dist + 1, 100)  # at - 150km, draw 50km radius circle
            e100 = np.nansum(kernel[dist-30:dist+30, dist:dist+100]>=1)/np.sum(np.isfinite(kernel[dist-30:dist+30, dist:dist+100]))

            # if e100 >= -5:  ### random LSTA p10
            #     return
            #
            if e100 <= 0.5:   ### random LSTA p90
                wetflag +=1

            # if e100 >= -3:  ### core LSTA p10
            #     return

            # if e100 <= 2.88:  ### random LSTA p90
            #     return

            # if e100 >= -1.36:  ### core LSTA p25
            #     return
            #
            # if e100 <= 1.45:  ### random LSTA p75
            #     return

            # if (e100 >= 1.1) | (e100<=1):
            #     return


        if wetflag != 0:
            return
        cnt = np.zeros_like(kernel)
        cnt[np.isfinite(kernel)] = 1
        smoutlist.append(kernel)
        scntlist.append(cnt)

    vdic = {}

    for d in era.data_vars:

        var = ua.cut_kernel(era[d].values,xpos, ypos,dist)
        vdic[d] = var

    cntera = np.zeros_like(kernel)
    cntera[np.isfinite(vdic[list(vdic.keys())[0]])] = 1


    if np.nansum(msg) > 0:
        probm = ua.cut_kernel(msg,xpos, ypos,dist)
        cntmsg = np.zeros_like(kernel)
        cntmsg[np.isfinite(probm)] = 1
    else:
        probm = np.zeros_like(kernel)
        cntmsg = np.zeros_like(kernel)

    if np.nansum(cmorph) > 0:
        probcm = ua.cut_kernel(cmorph,xpos, ypos,dist)
        cntcm = np.zeros_like(kernel)
        cntcm[np.isfinite(probcm)] = 1

    else:
        probcm = np.zeros_like(kernel)
        cntcm = np.zeros_like(kernel)


    return smoutlist,  scntlist, vdic,cntera, probm, cntmsg, probcm, cntcm

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

    pos = np.where((fi.values >= 5) & (fi.values <= 65)) #np.where((fi.values >= 5) & (fi.values <= 65)) #np.where((fi.values == 1)) #np.where((fi.values >= 5) & (fi.values <= 65)) #np.where((fi.values == 1)) #(fi.values >= 5) & (fi.values < 65)

    if (np.sum(pos) == 0):
        print('No blobs found')
        return None

    topo = xr.open_dataset(cnst.LSTA_TOPO)
    ttopo = topo['h']
    #
    grad = np.gradient(ttopo.values)
    gradsum = abs(grad[0])+abs(grad[1])

    smpath = [cnst.AMSRE_ANO_DAY + 'sma_' + fdate + '.nc',
              cnst.AMSRE_ANO_NIGHT + 'sma_' + fdate + '.nc',
              cnst.AMSRE_ANO_DAY + 'sma_' + pfdate + '.nc',
              cnst.AMSRE_ANO_NIGHT + 'sma_' + pfdate + '.nc'
              ]

    smlist = []
    smkernels = []
    smcnt = []
    dist = 200

    for sid , sp in enumerate(smpath):

        try:
            lsta = xr.open_dataset(sp)
        except OSError:
                return None
        print('Doing '+ sp)

        lsta = lsta.sel(lon=slice(-11, 11), lat=slice(9, 21))

        lsta_da = lsta['SM'].squeeze()

        if sid == 0:
            if (np.sum(np.isfinite(lsta_da)) / lsta_da.size) < 0.05:
                print('Not enough valid')
                return None

        try:
            lsta_da = topo.salem.transform(lsta_da)
        except RuntimeError:
            print('lsta_da on LSTA interpolation problem')
            return None


        lsta_da.values[ttopo.values>=450] = np.nan
        lsta_da.values[gradsum>30] = np.nan

        smlist.append(lsta_da)
        smkernels.append(np.zeros((dist*2+1, dist*2+1)))
        smcnt.append(np.zeros((dist*2 + 1, dist * 2 + 1)))
        del lsta

    if len(smlist)!=4:
        return None

    cnte_sum = np.zeros((dist*2+1, dist*2+1))
    cntm_sum = np.zeros((dist*2+1, dist*2+1))
    probm_sum = np.zeros((dist*2+1, dist*2+1))
    probc_sum = np.zeros((dist * 2 + 1, dist * 2 + 1))
    cntc_sum = np.zeros((dist * 2 + 1, dist * 2 + 1))

    edic = {}

    probs = get_previous_hours(storm_date,lsta_date, fi.attrs['eh'], fi.attrs['refhour'])
    print('Era5 collect')
    try:
        probs_on_lsta = topo.salem.transform(probs)
    except RuntimeError:
        print('Era5 on LSTA interpolation problem')
        return None
    del probs

    probs_msg = get_previous_hours_msg(storm_date,lsta_date, fi.attrs['eh'], fi.attrs['refhour'])
    probsm_on_lsta = topo.salem.transform(probs_msg, interp='nearest')
    del probs_msg

    probs_cm = get_previous_hours_CMORPH(storm_date)   # get previous rain to storm
    try:
        probscm_on_lsta = topo.salem.transform(probs_cm)
    except RuntimeError:
        return None
    del probs_cm
    del topo

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

        mhour = mcs_hour.sel(time=pd.datetime(int(fdate[0:4]), int(fdate[4:6]), int(fdate[6:8]), 14), lon=lon,
                             lat=lat).values
        if mhour <= 13:
            mhour += 24

        chour = fi['time.hour'].values

        if (chour >= 0) & (chour <= 13):
            chour += 24
        if (mhour < chour) | (np.isnan(mhour)):
            print('Core overlaps: earliest:', mhour, ' core: ', chour)
            continue

        point = (smlist[0]).sel(lat=lat, lon=lon, method='nearest')
        plat = point['lat'].values
        plon = point['lon'].values

        xpos = np.where((smlist[0])['lon'].values == plon)
        xpos = int(xpos[0])
        ypos = np.where((smlist[0])['lat'].values == plat)
        ypos = int(ypos[0])
        try:
            smoutlist, scntlist, vdic, cntera, probm, cntmsg, probcm, cntcm = cut_kernel(xpos, ypos, smlist, dist, era=probs_on_lsta, msg=probsm_on_lsta, cmorph=probscm_on_lsta)
        except TypeError:
            continue

        if np.nansum(probm[dist-50:dist+50,dist-30:dist+100])>=2:   # filter out cases with MCSs at 12 [dist-50:dist+50,dist-30:dist+100]
            print('Meteosat MCS continue')
            continue
        # #
        # if np.nanmax(vdic['tciw'][:,dist::])>=0.02:   # filter out cases with MCSs at 12
        #     print('ERA MCS continue')
        #     continue

        if np.nanmin((vdic['tciwmid'][dist-50:dist+50,dist-30:dist+100]))<=-0.4:   # 0.03 for tciw, -0.3 for w [dist-50:dist+50,dist-30:dist+100] ,filter out cases with MCSs at 12
            print('ERA MCS continue')
            continue

        # if np.nansum(probcm[100:200,:])>=100:   # filter out cases with rainfall to the south
        #     print('Southward rainfall, continue')
        #     continue

        for ids in np.arange(4):
            smkernels[ids] = np.nansum(np.stack([smkernels[ids], smoutlist[ids]]), axis=0)
            smcnt[ids] = np.nansum(np.stack([smcnt[ids], scntlist[ids]]), axis=0)


        cnte_sum = np.nansum(np.stack([cnte_sum, cntera]), axis=0)

        cntm_sum = np.nansum(np.stack([cntm_sum, cntmsg]), axis=0)
        probm_sum = np.nansum(np.stack([probm_sum, probm]), axis=0)

        probc_sum = np.nansum(np.stack([probc_sum, probcm]), axis=0)
        cntc_sum = np.nansum(np.stack([cntc_sum, cntcm]), axis=0)

        for ks in vdic.keys():
            if ks in edic:
                edic[ks] = np.nansum(np.stack([edic[ks], vdic[ks]]), axis=0)
            else:
                edic[ks] = vdic[ks]
        counter += 1
        print('Saved core ', counter)

    del mcs_hour

    # ipdb.set_trace()
    if np.sum(smcnt) == 0:
        return None

    if not edic:
        return None

    outlist = smkernels + smcnt + [cnte_sum, cntm_sum, probm_sum, cntc_sum,probc_sum]
    outnames = ['lsta0', 'lsta-1', 'lsta-2','lsta-3','cnt0','cnt-1','cnt-2','cnt-3', 'cnte', 'cntm', 'probmsg', 'cntc', 'probc']
    for ek in edic.keys():
        outnames.append(ek)
        outlist.append(edic[ek])

    print('Returning')

    return outlist, outnames


###################################################

###################################################
###################################################


def plot_AMSR_all(h, eh):

    dic = {}
    dic2 = {}

    name = "ERA5_composite_cores_AMSRE_500w04_15k_p90"#"ERA5_composite_cores_AMSRE_w1_15k_minusMean"

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

    extent = (dic['lsta0'].shape[1]-1)/2
    xlen = dic['lsta0'].shape[1]
    ylen = dic['lsta0'].shape[0]

    xv, yv = np.meshgrid(np.arange(ylen), np.arange(xlen))
    st=30
    xquiv = xv[4::st, 4::st]
    yquiv = yv[4::st, 4::st]

    u = (dic['u925']/ dic['cnte'])[4::st, 4::st]
    v = (dic['v925']/ dic['cnte'])[4::st, 4::st]

    u600 = (dic['u650']/ dic['cnte'])[4::st, 4::st]
    v600 = (dic['v650']/ dic['cnte'])[4::st, 4::st]

    u_orig = (dic['v925_orig']/ dic['cnte'])[4::st, 4::st]
    v_orig = (dic['v925_orig']/ dic['cnte'])[4::st, 4::st]


    f = plt.figure(figsize=(15,6), dpi=200)
    ax = f.add_subplot(241)

    plt.contourf((dic['lsta-3']) / (dic['cnt-3'] ), cmap='RdBu', levels=np.linspace(-1,1,10), extend='both') #-(rkernel2_sum / rcnt_sum)
    plt.plot(extent, extent, 'bo')
    plt.colorbar(label='%')
    #pdb.set_trace()

    contours = plt.contour((dic['v925']/ dic['cntc']), extend='both',colors='k', linewidths=4, levels=[-50,0,50])  #np.arange(-15,-10,0.5)
    #plt.clabel(contours, inline=True, fontsize=9, fmt='%1.1f')

    contours = plt.contour((dic['v925_orig']/ dic['cntc']), extend='both',colors='turquoise', linewidths=4, levels=[-50,0,50])  #np.arange(-15,-10,0.5)
    #plt.clabel(contours, inline=True, fontsize=9, fmt='%1.1f')
    plt.plot(extent, extent, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6 , dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.title('1700UTC | '+str(np.max(dic['cnt0']))+' cores, LSTA on day-1')#, CMORPH rainP>10mm [6am|day-1 to 10am|day0]', fontsize=9)


    ax1 = f.add_subplot(242)
    plt.contourf((dic['lsta-2'])/ dic['cnt-2'], extend='both',  cmap='RdBu',levels=np.linspace(-1,1,10)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='%')
    contours = plt.contour((dic['t'] / dic['cnte']), extend='both',levels=np.linspace(-0.8,0.8,9), cmap='PuOr_r') #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.1f')
    #qu = ax1.quiver(xquiv, yquiv, u, v, scale=50)
    plt.plot(extent, extent, 'bo')
    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    #plt.title('Shading: MSG LSTA day0, Contours: ERA5 925hPa T anomaly', fontsize=9)
#
    ax1 = f.add_subplot(243)
    plt.contourf((dic['lsta-1']) / (dic['cnt-1']), extend='both',  cmap='RdBu',levels=np.linspace(-1,1,10)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='g kg-1')
    contours = plt.contour((dic['shear'] / dic['cnte']), extend='both',levels=np.arange(-17,-12,0.25), cmap='viridis') #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    plt.plot(extent, extent, 'bo')
    #qu = ax1.quiver(xquiv, yquiv, u, v, scale=50)
    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    #plt.title('Shading: 925hPa q anomaly, Contours: 650hPa-925hPa zonal shear', fontsize=9)

    ax1 = f.add_subplot(244)
 #   plt.contourf(((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r', vmin=-1.5, vmax=1.5) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.contourf(((dic['lsta0'])/ dic['cnt0']), extend='both',  cmap='RdBu',levels=np.linspace(-1,1,10)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
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
    #plt.title('Shading: Divergence, vectors: 925hPa wind anomaly', fontsize=9)


    ax1 = f.add_subplot(245)
    #   plt.contourf(((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r', vmin=-1.5, vmax=1.5) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.contourf(dic['tciwmid'] / dic['cnte'], extend='both', cmap='RdBu',levels=np.linspace(-0.03,0.03,10))  # , levels=np.linspace(-2,2,10) (dic['theta_e']) / dic['cntp']), #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='K')
    plt.plot(extent, extent, 'bo')
    contours = plt.contour((dic['tciwlow'] / dic['cnte']), extend='both', cmap='viridis') #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    qu = ax1.quiver(xquiv, yquiv, u600, v600, scale=20)
    qk = plt.quiverkey(qu, 0.9, 0.02, 1, '1 m s$^{-1}$',
                       labelpos='E', coordinates='figure')

    ax1.set_xticklabels(np.array((np.linspace(0, extent * 2, 9) - 100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent * 2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    #plt.title(r'Shading: $\Delta \theta_{e}$ anomaly, contours: 925hPa u-wind anomaly, vectors: 650hPa wind anomaly', fontsize=9)

    ax1 = f.add_subplot(246)
    #   plt.contourf(((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r', vmin=-1.5, vmax=1.5) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.contourf(((dic['div'])/ dic['cnte'])*100, extend='both',  cmap='RdBu', levels=np.linspace(-0.5,0.5,10))  # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='K')
    plt.plot(extent, extent, 'bo')
    contours = plt.contour((dic['u650'] / dic['cnte']), extend='both', cmap='viridis') #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    qu = ax1.quiver(xquiv, yquiv, u_orig, v_orig, scale=20)
    qk = plt.quiverkey(qu, 0.9, 0.02, 1, '1 m s$^{-1}$',
                       labelpos='E', coordinates='figure')

    ax1.set_xticklabels(np.array((np.linspace(0, extent * 2, 9) - 100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent * 2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    #plt.title(r'Shading: 925hPa v wind, contours: 650hPa u-wind anomaly, vectors: 925hpa wind', fontsize=9)

    ax1 = f.add_subplot(247)
    #   plt.contourf(((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r', vmin=-1.5, vmax=1.5) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.contourf(((dic['sh']) / dic['cnte']) , extend='both', cmap='RdBu', levels=np.linspace(-20,20,10))  # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='K')
    plt.plot(extent, extent, 'bo')
    contours = plt.contour((dic['skt'] / dic['cnte']), extend='both', cmap='PuOr') #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    qu = ax1.quiver(xquiv, yquiv, u_orig, v_orig, scale=20)
    qk = plt.quiverkey(qu, 0.9, 0.02, 1, '1 m s$^{-1}$',
                       labelpos='E', coordinates='figure')

    ax1.set_xticklabels(np.array((np.linspace(0, extent * 2, 9) - 100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent * 2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    #plt.title(r'Shading: 925hPa v wind, contours: 650hPa u-wind anomaly, vectors: 925hpa wind', fontsize=9)

    ax1 = f.add_subplot(248)
    #   plt.contourf(((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r', vmin=-1.5, vmax=1.5) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.contourf(((dic['q'])*1000/ dic['cnte']), extend='both',  cmap='RdBu',levels=np.arange(-1,1.1,0.1))  # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='K')
    plt.plot(extent, extent, 'bo')
    contours = plt.contour((dic['u650'] / dic['cnte']), extend='both', cmap='viridis') #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    qu = ax1.quiver(xquiv, yquiv, u_orig, v_orig, scale=20)
    qk = plt.quiverkey(qu, 0.9, 0.02, 1, '1 m s$^{-1}$',
                       labelpos='E', coordinates='figure')

    ax1.set_xticklabels(np.array((np.linspace(0, extent * 2, 9) - 100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent * 2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    #plt.title(r'Shading: 925hPa v wind, contours: 650hPa u-wind anomaly, vectors: 925hpa wind', fontsize=9)

    plt.tight_layout()
    plt.show()
    plt.savefig(cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/plots/"+name+str(h).zfill(2)+'_'+str(eh).zfill(2)+'.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png)



def plot_AMSR_small(h, eh):

    dic = {}
    dic2 = {}

    name = "ERA5_composite_cores_AMSRE_500w04_15k_minusMean"#"ERA5_composite_cores_AMSRE_w1_15k_minusMean"

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

    extent = (dic['lsta0'].shape[1]-1)/2
    xlen = dic['lsta0'].shape[1]
    ylen = dic['lsta0'].shape[0]

    xv, yv = np.meshgrid(np.arange(ylen), np.arange(xlen))
    st=30
    xquiv = xv[4::st, 4::st]
    yquiv = yv[4::st, 4::st]

    u = (dic['u925']/ dic['cnte'])[4::st, 4::st]
    v = (dic['v925']/ dic['cnte'])[4::st, 4::st]

    u600 = (dic['u650']/ dic['cnte'])[4::st, 4::st]
    v600 = (dic['v650']/ dic['cnte'])[4::st, 4::st]

    u_orig = (dic['v925_orig']/ dic['cnte'])[4::st, 4::st]
    v_orig = (dic['v925_orig']/ dic['cnte'])[4::st, 4::st]


    f = plt.figure(figsize=(10,7), dpi=200)
    ax = f.add_subplot(221)

    plt.contourf((( dic['lsta-2'] +dic['lsta-3'])  / ( dic['cnt-2']+dic['cnt-3'])), cmap='RdBu', levels=np.linspace(-1.2,1.2,10), extend='both') #-(rkernel2_sum / rcnt_sum)
    plt.plot(extent, extent, 'bo')
    plt.colorbar(label='%')
    #pdb.set_trace()

    contours = plt.contour((dic['v925']) / dic['cnte'] , extend='both', cmap='RdBu', levels=np.linspace(-1,1,9)) # , levels=np.linspace(-1,1,10)#(dic['probc']/ dic['cntc'])*100, extend='both', levels=np.arange(15,70,12), cmap='jet') # #, levels=np.arange(1,5, 0.5)
    plt.clabel(contours, inline=True, fontsize=11, fmt='%1.1f')
    plt.plot(extent, extent, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6 , dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.title('1700UTC | '+str(np.max(dic['cnt0']))+' cores, LSTA on day-1')#, CMORPH rainP>10mm [6am|day-1 to 10am|day0]', fontsize=9)

#
    ax1 = f.add_subplot(222)
    plt.contourf(((dic['lsta0']) / (dic['cnt0'])), extend='both',  cmap='RdBu',levels=np.linspace(-1.2,1.2,10)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='g kg-1')
    contours = plt.contour(((dic['div'])/ dic['cnte'])*100, extend='both',  cmap='PuOr', levels=np.linspace(-0.7,0.7,9)) #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    plt.plot(extent, extent, 'bo')
    qu = ax1.quiver(xquiv, yquiv, u, v, scale=15)
    qk = plt.quiverkey(qu, 0.9, 0.02,1, '1 m s$^{-1}$',
                       labelpos='E', coordinates='figure')
    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    #plt.title('Shading: 925hPa q anomaly, Contours: 650hPa-925hPa zonal shear', fontsize=9)


    ax1 = f.add_subplot(223)
    #   plt.contourf(((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r', vmin=-1.5, vmax=1.5) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.contourf(((dic['div'])/ dic['cnte'])*100, extend='both',  cmap='PuOr', levels=np.linspace(-0.7,0.7,10))  # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='K')
    plt.plot(extent, extent, 'bo')
    contours = plt.contour((dic['u650'] / dic['cnte']), extend='both', cmap='viridis') #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    qu = ax1.quiver(xquiv, yquiv, u_orig, v_orig, scale=20)
    qk = plt.quiverkey(qu, 0.9, 0.02, 1, '1 m s$^{-1}$',
                       labelpos='E', coordinates='figure')

    ax1.set_xticklabels(np.array((np.linspace(0, extent * 2, 9) - 100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent * 2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    #plt.title(r'Shading: 925hPa v wind, contours: 650hPa u-wind anomaly, vectors: 925hpa wind', fontsize=9)


    ax1 = f.add_subplot(224)
    #   plt.contourf(((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r', vmin=-1.5, vmax=1.5) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.contourf(((dic['q'])*1000/ dic['cnte']), extend='both',  cmap='RdBu',levels=np.arange(-1,1.1,0.1))  # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='K')
    plt.plot(extent, extent, 'bo')
    contours = plt.contour(((dic['v925_orig']) / dic['cnte']) , extend='both', cmap='RdBu', levels=np.linspace(-1,1,11))
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    qu = ax1.quiver(xquiv, yquiv, u_orig, v_orig, scale=20)
    qk = plt.quiverkey(qu, 0.9, 0.02, 1, '1 m s$^{-1}$',
                       labelpos='E', coordinates='figure')

    ax1.set_xticklabels(np.array((np.linspace(0, extent * 2, 9) - 100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent * 2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    #plt.title(r'Shading: 925hPa v wind, contours: 650hPa u-wind anomaly, vectors: 925hpa wind', fontsize=9)

    plt.tight_layout()
    plt.show()
    plt.savefig(cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/plots/"+name+str(h).zfill(2)+'_'+str(eh).zfill(2)+'_small.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png)
    plt.close()


def plot_AMSR_itd(h, eh):

    dic = {}
    dic2 = {}

    name = "ERA5_composite_cores_AMSRE_500w04_15k_minusMean" #ERA5_composite_cores_AMSRE_500w04_15k_minusMean"#"ERA5_composite_cores_AMSRE_w1_15k_minusMean"

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

    extent = (dic['lsta0'].shape[1]-1)/2
    xlen = dic['lsta0'].shape[1]
    ylen = dic['lsta0'].shape[0]

    xv, yv = np.meshgrid(np.arange(ylen), np.arange(xlen))
    st=30
    xquiv = xv[4::st, 4::st]
    yquiv = yv[4::st, 4::st]

    u = (dic['u925']/ dic['cnte'])[4::st, 4::st]
    v = (dic['v925']/ dic['cnte'])[4::st, 4::st]

    u600 = (dic['u650']/ dic['cnte'])[4::st, 4::st]
    v600 = (dic['v650']/ dic['cnte'])[4::st, 4::st]

    u_orig = (dic['v925_orig']/ dic['cnte'])[4::st, 4::st]
    v_orig = (dic['v925_orig']/ dic['cnte'])[4::st, 4::st]

    levels = list(np.arange(-1.2, 0, 0.3)) + list(np.arange(0.3, 1.5, 0.3))
    f = plt.figure(figsize=(10, 8), dpi=200)
    ax = f.add_subplot(221)

    plt.title(r'1700UTC | 2707 cores, NIGHT-1 SM anomaly')#, fontsize=9)

    plt.contourf(((dic['lsta-3'])  / (dic['cnt-3'])), cmap='RdBu', levels=levels, extend='both') #-(rkernel2_sum / rcnt_sum)
    plt.plot(extent, extent, 'bo')
    plt.axvline(x=200, linestyle='dashed', color='k', linewidth=1)
    plt.axhline(y=200, linestyle='dashed', color='k', linewidth=1)
    plt.colorbar(label='%')
    #pdb.set_trace()

    #contours = plt.contour((dic['v925']) / dic['cnte'] , extend='both', cmap='RdBu', levels=np.linspace(-1,1,9)) # , levels=np.linspace(-1,1,10)#(dic['probc']/ dic['cntc'])*100, extend='both', levels=np.arange(15,70,12), cmap='jet') # #, levels=np.arange(1,5, 0.5)
    #plt.clabel(contours, inline=True, fontsize=11, fmt='%1.1f')
    plt.plot(extent, extent, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6 , dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    #, CMORPH rainP>10mm [6am|day-1 to 10am|day0]', fontsize=9)

#
    ax1 = f.add_subplot(222)
    plt.contourf(((dic['lsta-2']) / (dic['cnt-2'])), extend='both',  cmap='RdBu',levels=levels) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='%')
    # contours = plt.contour((dic['v925']/ dic['cnte']), extend='both',colors='k', linewidths=4, levels=[-50,0,50])  #np.arange(-15,-10,0.5)
    # #plt.clabel(contours, inline=True, fontsize=9, fmt='%1.1f')
    #
    # contours = plt.contour((dic['v925_orig']/ dic['cnte']), extend='both',colors='turquoise', linewidths=4, levels=[-50,0,50])  #np.arange(-15,-10,0.5)
    # #plt.clabel(contours, inline=True, fontsize=9, fmt='%1.1f')
    plt.plot(extent, extent, 'bo')
    plt.axvline(x=200, linestyle='dashed', color='k', linewidth=1)
    plt.axhline(y=200, linestyle='dashed', color='k', linewidth=1)
    # qu = ax1.quiver(xquiv, yquiv, u, v, scale=15)
    # qk = plt.quiverkey(qu, 0.9, 0.02,1, '1 m s$^{-1}$',
    #                    labelpos='E', coordinates='figure')
    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('DAY-1: SM anomaly')#, vectors: 925hPa wind anomaly (Day0)', fontsize=9)


    ax1 = f.add_subplot(223)
    #   plt.contourf(((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r', vmin=-1.5, vmax=1.5) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.contourf(((dic['lsta-1']) / (dic['cnt-1'])), extend='both', cmap='RdBu', levels=levels) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='%')
    plt.plot(extent, extent, 'bo')
    plt.axvline(x=200, linestyle='dashed', color='k', linewidth=1)
    plt.axhline(y=200, linestyle='dashed', color='k', linewidth=1)
    # contours = plt.contour((dic['t'] / dic['cnte']), extend='both',
    #                        levels=[-0.8, -0.7, -0.6, -0.5, -0.4, -0.2, 0, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8], cmap='RdBu_r',
    #                        linewidths=2) #np.arange(-15,-10,0.5)
    # plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    # qu = ax1.quiver(xquiv, yquiv, u_orig, v_orig, scale=20)
    # qk = plt.quiverkey(qu, 0.9, 0.02, 1, '1 m s$^{-1}$',
    #                    labelpos='E', coordinates='figure')

    ax1.set_xticklabels(np.array((np.linspace(0, extent * 2, 9) - 100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent * 2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title(r'NIGHT0: SM anomaly')#, contours: 925hPa T anomaly (Day0)', fontsize=9)


    ax1 = f.add_subplot(224)
    #   plt.contourf(((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r', vmin=-1.5, vmax=1.5) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.contourf(((dic['lsta0']) / (dic['cnt0'])), extend='both', cmap='RdBu', levels=levels)  # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='%')
    plt.plot(extent, extent, 'bo')
    plt.axvline(x=200, linestyle='dashed', color='k', linewidth=1)
    plt.axhline(y=200, linestyle='dashed', color='k', linewidth=1)
    # contours = plt.contour(((dic['div']) / dic['cnte']) * 100, extend='both', cmap='PuOr',
    #                        levels=np.linspace(-0.7, 0.7, 9))  # np.arange(-15,-10,0.5)
    # plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    # qu = ax1.quiver(xquiv, yquiv, u, v, scale=15)
    # qk = plt.quiverkey(qu, 0.9, 0.02, 1, '1 m s$^{-1}$',
    #                    labelpos='E', coordinates='figure')

    ax1.set_xticklabels(np.array((np.linspace(0, extent * 2, 9) - 100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent * 2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title(r'DAY0: SM anomaly')#, contours: 925hPa divergence', fontsize=9)

    plt.tight_layout()
    plt.show()
    plt.savefig(cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/plots/"+name+str(h).zfill(2)+'_'+str(eh).zfill(2)+'_small_ITD.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png)
    plt.close()


def plot_AMSR_itd_small(h, eh):

    dic = {}
    dic2 = {}

    name = "ERA5_composite_cores_AMSRE_500w04_15k_minusMean" #ERA5_composite_cores_AMSRE_500w04_15k_minusMean"#"ERA5_composite_cores_AMSRE_w1_15k_minusMean"

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

    for k in dic.keys():
        dic[k]=dic[k][100:301,100:301]

    extent = (dic['lsta0'].shape[1]-1)/2
    xlen = dic['lsta0'].shape[1]
    ylen = dic['lsta0'].shape[0]

    xv, yv = np.meshgrid(np.arange(ylen), np.arange(xlen))
    st=30
    xquiv = xv[4::st, 4::st]
    yquiv = yv[4::st, 4::st]

    u = (dic['u925']/ dic['cnte'])[4::st, 4::st]
    v = (dic['v925']/ dic['cnte'])[4::st, 4::st]

    u600 = (dic['u650']/ dic['cnte'])[4::st, 4::st]
    v600 = (dic['v650']/ dic['cnte'])[4::st, 4::st]

    u_orig = (dic['v925_orig']/ dic['cnte'])[4::st, 4::st]
    v_orig = (dic['v925_orig']/ dic['cnte'])[4::st, 4::st]

    levels = list(np.arange(-1.5, 0, 0.25)) + list(np.arange(0.3, 1.75, 0.25))
    f = plt.figure(figsize=(10, 8), dpi=200)
    ax = f.add_subplot(221)

    plt.title(r'1700UTC | 2707 cores, NIGHT-1 SM anomaly')#, fontsize=9)

    plt.contourf(((dic['lsta-3'])  / (dic['cnt-3'])), cmap='RdBu', levels=levels, extend='both') #-(rkernel2_sum / rcnt_sum)
    plt.plot(extent, extent, 'bo')
    plt.axvline(x=100, linestyle='dashed', color='k', linewidth=1)
    plt.axhline(y=100, linestyle='dashed', color='k', linewidth=1)
    plt.colorbar(label='%')
    #pdb.set_trace()

    #contours = plt.contour((dic['v925']) / dic['cnte'] , extend='both', cmap='RdBu', levels=np.linspace(-1,1,9)) # , levels=np.linspace(-1,1,10)#(dic['probc']/ dic['cntc'])*100, extend='both', levels=np.arange(15,70,12), cmap='jet') # #, levels=np.arange(1,5, 0.5)
    #plt.clabel(contours, inline=True, fontsize=11, fmt='%1.1f')
    plt.plot(extent, extent, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6 +300 , dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    #, CMORPH rainP>10mm [6am|day-1 to 10am|day0]', fontsize=9)

#
    ax1 = f.add_subplot(222)
    plt.contourf(((dic['lsta-2']) / (dic['cnt-2'])), extend='both',  cmap='RdBu',levels=levels) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='%')
    # contours = plt.contour((dic['v925']/ dic['cnte']), extend='both',colors='k', linewidths=4, levels=[-50,0,50])  #np.arange(-15,-10,0.5)
    # #plt.clabel(contours, inline=True, fontsize=9, fmt='%1.1f')
    #
    # contours = plt.contour((dic['v925_orig']/ dic['cnte']), extend='both',colors='turquoise', linewidths=4, levels=[-50,0,50])  #np.arange(-15,-10,0.5)
    # #plt.clabel(contours, inline=True, fontsize=9, fmt='%1.1f')
    plt.plot(extent, extent, 'bo')
    plt.axvline(x=100, linestyle='dashed', color='k', linewidth=1)
    plt.axhline(y=100, linestyle='dashed', color='k', linewidth=1)
    # qu = ax1.quiver(xquiv, yquiv, u, v, scale=15)
    # qk = plt.quiverkey(qu, 0.9, 0.02,1, '1 m s$^{-1}$',
    #                    labelpos='E', coordinates='figure')
    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6+300, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('DAY-1: SM anomaly')#, vectors: 925hPa wind anomaly (Day0)', fontsize=9)


    ax1 = f.add_subplot(223)
    #   plt.contourf(((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r', vmin=-1.5, vmax=1.5) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.contourf(((dic['lsta-1']) / (dic['cnt-1'])), extend='both', cmap='RdBu', levels=levels) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='%')
    plt.plot(extent, extent, 'bo')
    plt.axvline(x=100, linestyle='dashed', color='k', linewidth=1)
    plt.axhline(y=100, linestyle='dashed', color='k', linewidth=1)
    # contours = plt.contour((dic['t'] / dic['cnte']), extend='both',
    #                        levels=[-0.8, -0.7, -0.6, -0.5, -0.4, -0.2, 0, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8], cmap='RdBu_r',
    #                        linewidths=2) #np.arange(-15,-10,0.5)
    # plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    # qu = ax1.quiver(xquiv, yquiv, u_orig, v_orig, scale=20)
    # qk = plt.quiverkey(qu, 0.9, 0.02, 1, '1 m s$^{-1}$',
    #                    labelpos='E', coordinates='figure')

    ax1.set_xticklabels(np.array((np.linspace(0, extent * 2, 9) - 100) * 6+300, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent * 2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title(r'NIGHT0: SM anomaly')#, contours: 925hPa T anomaly (Day0)', fontsize=9)


    ax1 = f.add_subplot(224)
    #   plt.contourf(((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r', vmin=-1.5, vmax=1.5) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.contourf((dic['lsta0']/ dic['cnt0']) , extend='both', cmap='RdBu', levels=levels)  # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='%')
    plt.plot(extent, extent, 'bo')
    plt.axvline(x=100, linestyle='dashed', color='k', linewidth=1)
    plt.axhline(y=100, linestyle='dashed', color='k', linewidth=1)
    # contours = plt.contour(((dic['div']) / dic['cnte']) * 100, extend='both', cmap='PuOr',
    #                        levels=np.linspace(-0.7, 0.7, 9))  # np.arange(-15,-10,0.5)
    # plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    # qu = ax1.quiver(xquiv, yquiv, u, v, scale=15)
    # qk = plt.quiverkey(qu, 0.9, 0.02, 1, '1 m s$^{-1}$',
    #                    labelpos='E', coordinates='figure')

    ax1.set_xticklabels(np.array((np.linspace(0, extent * 2, 9) - 100) * 6+300, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent * 2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title(r'DAY0: SM anomaly')#, contours: 925hPa divergence', fontsize=9)

    plt.tight_layout()
    plt.show()
    plt.savefig(cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/plots/"+name+str(h).zfill(2)+'_'+str(eh).zfill(2)+'_300km_ITD.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png)
    plt.close()
