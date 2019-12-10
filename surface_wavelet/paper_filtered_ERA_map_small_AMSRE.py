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
import glob
import os
from collections import OrderedDict
import salem
from utils import u_met, u_parallelise, u_gis, u_arrays as ua, constants as cnst, u_grid, u_darrays
from scipy.interpolate import griddata
import multiprocessing

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
            hlist.append(hh, 12-(hh+24))

    for l in hlist:
        print('Doing '+str(l))
        composite(l[0], l[1])



def eh_loop():

    all = np.arange(-38,4,3)
    #all= np.arange(4,30,3)

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

    msgopen = pd.read_csv(cnst.network_data + 'figs/LSTA/corrected_LSTA/new/ERA5/core_txt/cores_gt15000km2_table_AMSRE_tracking_'+str(h)+'.csv')
    hour = h
    msg = pd.DataFrame.from_dict(msgopen)# &  &
    msg['eh'] = eh
    msg['refhour'] = h

    msg['date'] = pd.to_datetime(msg[['year','month','day']])
    print('Start core number ', len(msg))

    #msgin = msg[(msg['SMmean0']<-3) & (msg['SMmean-1']<-1.5) & (msg['SMmean0']>-900)& (msg['SMmean-1']>-900)] #[(msg['initTime']>=12.5) & (msg['initTime']<h)]#[msg['SMwet']==2]
    #msgin = msg[(msg['SMdry0'] ==1) & (msg['SMdry-1']==1) & (msg['SMmean0']>-900)& (msg['SMmean-1']>-900)& (msg['SMmean-1']<-1.8)]
    msgin = msg[(msg['SMwet0'] == 1) & (msg['SMwet-1'] == 1) & (msg['SMmean0'] > 0.3) & (msg['SMmean-1'] > 0.9)]
    #msgin = msg[(msg['SMmean0'] > 0.01) & (msg['SMmean-1'] > 0.01)]
    print('Number of cores', len(msgin))

    # calculate the chunk size as an integer
    #'chunk_size = int(msg.shape[0] / pnumber)
    msgin.sort_values(by='date')

    for year in np.arange(2006, 2011):

        msgy = msgin[msgin['year']==year]

        chunk, chunk_ind, chunk_count = np.unique(msgy.date, return_index=True, return_counts=True)

        chunks = [msgy.ix[msgy.index[ci:ci + cc]] for ci, cc in zip(chunk_ind, chunk_count)] # daily chunks

        # res = []
        # for m in chunks[0:30]:
        #     out = file_loop(m)
        #     res.append(out)
        #
        # ipdb.set_trace()
        # return
        dic = u_parallelise.era_run_arrays(4, file_loop, chunks)

        pkl.dump(dic, open(cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/core_txt/ERA5_cores_WET_TAG_noMeteosatFilter_AMSRE"+str(eh) + "UTCERA"+str(hour).zfill(2)+'_'+str(year)+".p", "wb"))
        del dic
        print('Dumped file')



def cut_kernel(xpos, ypos, arrlist, dist, era=False, msg=False, cmorph=False):

    smoutlist  = []
    scntlist = []
    wetflag = 0
    for ids, arr in enumerate(arrlist):
        kernel = ua.cut_kernel(arr,xpos, ypos,dist)
        kernel = kernel #- np.nanmean(kernel)
        if kernel.shape != (dist*2+1, dist*2+1):
            print('Kernels shape wrong!')
            ipdb.set_trace()
        if ids == 0:
            if (np.sum(np.isfinite(kernel)) < 2):
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


    cmp = cmp.sel(longitude=slice(-13.5, 13.5), latitude=slice(8, 21))
    #css = css.sel(longitude=slice(-13, 13), latitude=slice(8, 21))
    pl_clim = pl_clim.sel(longitude=slice(-13.5, 13.5), latitude=slice(8, 21))

    #srfc_clim = srfc_clim.sel(longitude=slice(-13, 13), latitude=slice(8, 21))
    esrfc_clim = esrfc_clim.sel(longitude=slice(-13.5, 13.5), latitude=slice(8, 21))

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

    uwind_uup = cmm['u'].sel(level=700).squeeze() #- pl_clim['u'].sel(level=650).squeeze()
    vwind_uup = cmm['v'].sel(level=700).squeeze() #- pl_clim['v'].sel(level=650).squeeze()
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
    theta_mid = u_met.theta_e(650, cmm['t'].sel(level=650).squeeze().values - 273.15,
                                 cmm['q'].sel(level=650).squeeze()) - u_met.theta_e(650, pl_clim['t'].sel(level=650).squeeze().values - 273.15,
                                 pl_clim['q'].sel(level=650).squeeze())

    theta_diff = u_met.theta(925, cmm['t'].sel(level=925).squeeze().values - 273.15) - \
                   u_met.theta(650, cmm['t'].sel(level=650).squeeze().values - 273.15)

    # theta_diff_clim = u_met.theta(925, pl_clim['t'].sel(level=925).squeeze().values - 273.15) - \
    #                u_met.theta(650, pl_clim['t'].sel(level=650).squeeze().values - 273.15)
    #
    #
    # theta = theta_diff - theta_diff_clim

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

    cm['u700_orig'] = uwind_uup
    cm['v700_orig'] = vwind_uup

    cm['u650'] = uwind_up_ano
    cm['v650'] = vwind_up_ano
    cm['sh'] = sh
    cm['ev'] = ev
    cm['skt'] = skt

    cm['div'] = div *1000
    cm['q'] = q
    cm['t'] = t
    cm['theta_e'] = theta_e
    cm['theta_mid'] = theta_mid
    #cm['theta'] = theta

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


def file_loop(df):

    date = df['date'].iloc[0]
    eh = df['eh'].iloc[0]
    hour = df['hour'].iloc[0]
    print('Doing day: ', date)

    storm_date = date

    dayd = pd.Timedelta('1 days')

    if (hour) <= 13:
        print('Nighttime')
        lsta_date = storm_date - dayd
    else:
        print('Daytime')
        lsta_date = storm_date

    prev_lsta_date = lsta_date - dayd

    fdate = str(lsta_date.year) + str(lsta_date.month).zfill(2) + str(lsta_date.day).zfill(2)
    pfdate = str(prev_lsta_date.year) + str(prev_lsta_date.month).zfill(2) + str(prev_lsta_date.day).zfill(2)

    topo = xr.open_dataset(cnst.WA_TOPO_3KM)
    topo = topo.sel(lon=slice(-13.5, 13.5), lat=slice(7.8, 23.2))
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

        lsta = lsta.sel(lon=slice(-13, 13), lat=slice(8, 23))

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
    try:
        probs = get_previous_hours(storm_date,lsta_date, eh, hour)
    except:
        print('ERA failed')
        return
    print('Era5 collect')
    try:
        probs_on_lsta = topo.salem.transform(probs)
    except RuntimeError:
        print('Era5 on LSTA interpolation problem')
        return None
    del probs


    probs_msg = get_previous_hours_msg(storm_date,lsta_date, eh, hour)
    probsm_on_lsta = topo.salem.transform(probs_msg, interp='nearest')
    del probs_msg

    probs_cm = get_previous_hours_CMORPH(storm_date)   # get previous rain to storm
    try:
        probscm_on_lsta = topo.salem.transform(probs_cm)
    except RuntimeError:
        return None
    del probs_cm
    del topo

    counter=0

    #for lat, lon in zip(df.lat, df.lon):
    for dids, dit in df.iterrows():

        lat = dit.lat
        lon = dit.lon

        #initiation filter:
        # initpath = cnst.network_data + 'data/OBS/MSG_WA30/track_back_cores_vn1_'+str(hour)+'Z.txt'
        # if os.path.isfile(initpath):
            # dic = pd.read_table(initpath, delim_whitespace=True, header=None,
            #                     names=['year', 'mon', 'day', 'i_core', 'j_core', 'i_initiation', 'j_initiation',
            #                            'core_time', 'initiation_time'])
            # ddic = dic[
            #     (dic['i_core'] == dit['xloc']) & (dic['j_core'] == dit['yloc']) & (dic['year'] == dit['year']) & (
            #                 dic['mon'] == dit['month']) & (dic['day'] == dit['day'])]
            #
            # #ipdb.set_trace()
            # if len(ddic) == 0:
            #     continue
            # #ipdb.set_trace()
            # if (ddic['initiation_time'].values >3 ) | (ddic['initiation_time'].values < 0  ): #& (ddic['initiation_time'].values >= 12):
            #     continue

            # if (ddic['initiation_time'].values <=12 ) | (ddic['initiation_time'].values >= hour ): #& (ddic['initiation_time'].values >= 12):
            #     continue

            # initiation filter:

            # if (dit['xdiff'] < 100) & (dit['initTime']!=2):
            #     # & (ddic['initiation_time'].values >= 12):
            #     print('Initiation point too close')
            #     continue

                # if np.abs((dit['xdiff']) > 33) | (dit['xinit'] < 0):
                #
                #     print('Initiation point too far', dit['initTime'])
                #     continue



        # msg_latlon = np.load(cnst.network_data + 'data/OBS/MSG_WA30/MSG_1640_580_lat_lon.npz')
        #
        # lon = (msg_latlon['lon'])[ddic['j_initiation'], ddic['i_initiation']]
        # lat = (msg_latlon['lat'])[ddic['j_initiation'], ddic['i_initiation']]
        try:
            point = lsta_da.sel(lat=lat, lon=lon, method='nearest', tolerance=0.05)
        except KeyError:
            continue
        plat = point['lat'].values
        plon = point['lon'].values


        xpos = np.where(lsta_da['lon'].values == plon)
        xpos = int(xpos[0])
        ypos = np.where(lsta_da['lat'].values == plat)
        ypos = int(ypos[0])


        try:
            smoutlist, scntlist, vdic, cntera, probm, cntmsg, probcm, cntcm = cut_kernel(xpos, ypos, smlist, dist, era=probs_on_lsta, msg=probsm_on_lsta, cmorph=probscm_on_lsta)
        except TypeError:
            print('TypeError')
            continue
        #
        # if np.nansum(probm[dist-30:dist+30,dist-30:dist+67])>=2:   # filter out cases with MCSs at 12 [dist-50:dist+50,dist-30:dist+100]
        #     print('Meteosat MCS continue')
        #     continue

        if np.nanmin((vdic['tciwmid'][dist-30:dist+30,dist-30:dist+67]))<=-0.4:   # 0.03 for tciw, -0.3 for w [dist-50:dist+50,dist-30:dist+100] ,filter out cases with MCSs at 12
            print('ERA MCS continue')
            continue



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



    # ipdb.set_trace()
    if np.sum(smcnt) == 0:
        return None

    if not edic:
        return None

    outlist = smkernels + smcnt + [cnte_sum, cntm_sum, probm_sum, cntc_sum, probc_sum]
    outnames = ['lsta0', 'lsta-1', 'lsta-2', 'lsta-3', 'cnt0', 'cnt-1', 'cnt-2', 'cnt-3', 'cnte', 'cntm', 'probmsg',
                'cntc', 'probc']
    for ek in edic.keys():
        outnames.append(ek)
        outlist.append(edic[ek])

    print('Returning')

    return outlist, outnames


def plot_doug_all(h, eh):

    dic = {}
    dic2 = {}

    name = "ERA5_cores_ALL_WET_noMeteosatFilter_AMSRE"#"ERA5_composite_cores_AMSRE_w1_15k_minusMean"


    def coll(dic, h, eh, year):
        print(h)
        core = pkl.load(open(
            cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/core_txt/"+name+str(eh) + "UTCERA"+str(h).zfill(2)+'_'+str(year)+".p", "rb"))
        for id, k in enumerate(core.keys()):
            try:
                dic[k] = dic[k] + core[k]
            except KeyError:
                dic[k] = core[k]


    # for y in range(2006,2011):
    #     for h, eh in zip([16,17,18],[-4,-5,-6]):
    #         coll(dic, h, eh, y)

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

    u600 = (dic['u650_orig']/ dic['cnte'])[4::st, 4::st]
    v600 = (dic['v650_orig']/ dic['cnte'])[4::st, 4::st]

    u_orig = (dic['v925_orig']/ dic['cnte'])[4::st, 4::st]
    v_orig = (dic['v925_orig']/ dic['cnte'])[4::st, 4::st]


    f = plt.figure(figsize=(15,7))
    ax = f.add_subplot(231)


    plt.contourf((dic['lsta-2'] / dic['cnt-2']), cmap='RdBu_r', levels=np.linspace(-2,2,16), extend='both') #-(rkernel2_sum / rcnt_sum)
    plt.plot(extent, extent, 'bo')
    plt.colorbar(label='K', format='%1.2f')
    plt.text(0.02,0.08, 'ITD 0-line', color='turquoise', fontsize=12, transform=ax.transAxes)
    plt.text(0.02, 0.03, 'ITD anomaly 0-line', color='k', fontsize=12, transform=ax.transAxes)

    # plt.annotate('ITD 0-line', xy=(0.04, 0.1), xytext=(0, 4), size=15, color='turquoise', xycoords=('figure fraction', 'figure fraction'))
    #              #             textcoords='offset points')   #transform=ax.transAxes
    #pdb.set_trace()

    contours = plt.contour((dic['t'] / dic['cnte']), extend='both',levels=[ -0.8,-0.7, -0.6,-0.5,-0.4,-0.2,0, 0.2,0.4,0.5,0.6, 0.7, 0.8], cmap='PuOr_r', linewidths=2) # #, levels=np.arange(1,5, 0.5)
    plt.clabel(contours, inline=True, fontsize=11, fmt='%1.1f')
    #contours2 = plt.contour((dic['v925']) / dic['cntp'], extend='both', cmap='RdBu', levels=np.linspace(-1, 1,9))  # , levels=np.linspace(-1,1,10)#(dic['probc']/ dic['cntc'])*100, extend='both', levels=np.arange(15,70,12), cmap='jet') # #, levels=np.arange(1,5, 0.5)
    #plt.clabel(contours2, inline=True, fontsize=11, fmt='%1.0f')

    contours = plt.contour((dic['v925']/ dic['cnte']), extend='both',colors='k', linewidths=4, levels=[-50,0,50])  #np.arange(-15,-10,0.5)
    #plt.clabel(contours, inline=True, fontsize=9, fmt='%1.1f')

    contours = plt.contour((dic['v925_orig']/ dic['cnte']), extend='both',colors='turquoise', linewidths=4, levels=[-50,0,50])  #np.arange(-15,-10,0.5)
    #plt.clabel(contours, inline=True, fontsize=9, fmt='%1.1f')

    plt.plot(extent, extent, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6 , dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.title(str(h).zfill(2)+'00UTC | '+str(np.max(dic['cnt0']))+' cores, LSTA day-1, ERA5$_{noon}$ 925hPa T anomaly', fontsize=9)


    ax1 = f.add_subplot(232)
    plt.contourf(((dic['lsta0'])/ dic['cnt0']), extend='both',  cmap='RdBu_r', levels=np.linspace(-2,2,16)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
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
    plt.contourf(((dic['q'])*1000/ dic['cnte']), extend='both',  cmap='RdBu',levels=np.linspace(-0.9,0.9,16)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='g kg-1')
    # contours = plt.contour((dic['shear'] / dic['cnte']), extend='both',levels=np.arange(-17,-12,0.5), cmap='viridis_r', linewidths=2) #np.arange(-15,-10,0.5)
    # plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    contours = plt.contour((dic['u650'] / dic['cnte']), extend='both', levels=np.arange(-12, -7, 0.5),
                           cmap='viridis_r', linewidths=2)  # np.arange(-15,-10,0.5)
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
    plt.contourf(((dic['div'])/ dic['cnte'])*100, extend='both',  cmap='PuOr', levels=np.linspace(-0.7,0.7,10)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='10$^{-2}$ s$^{-1}$', format='%1.3f')
    plt.plot(extent, extent, 'bo')
    contours = plt.contour((dic['skt'] / dic['cnte']), extend='both', cmap='RdBu_r') #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    # qu = ax1.quiver(xquiv, yquiv, u, v, scale=15)
    #qk = plt.quiverkey(qu, 0.2, 0.02,1, '1 m s$^{-1}$',
    #                   labelpos='E', coordinates='figure')

    ax1.streamplot(xv, yv, (dic['u925']/ dic['cnte']), (dic['v925']/ dic['cnte']), density=[0.5, 1])

    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('925hPa divergence, vectors: 925hPa wind anomaly', fontsize=9)


    ax1 = f.add_subplot(235)
    #   plt.contourf(((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r', vmin=-1.5, vmax=1.5) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
     # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.contourf((dic['probmsg'] / dic['cntm']), extend='both', cmap='PuOr', levels=np.linspace(-0.02, 0.02, 10))
    plt.colorbar(label=r'Pa s$^{-1}$', format='%1.3f')
    plt.plot(extent, extent, 'bo')
    # contours = plt.contour(((dic['v925_orig']) / dic['cntp']) , extend='both', cmap='RdBu', levels=np.linspace(-2,2,11), linewidths=2) #np.arange(-15,-10,0.5)
    # plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    ax1.streamplot(xv, yv, dic['u925_orig']/ dic['cnte'], dic['v925_orig']/ dic['cnte'], density=[0.5, 1])

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
    plt.contourf(((dic['theta_e']) / dic['cnte']), extend='both', cmap='RdBu', levels=np.linspace(-3,3,18)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label=r'K')
    plt.plot(extent, extent, 'bo')
    contours = plt.contour((dic['u650'] / dic['cnte']), extend='both', cmap='PuOr', levels=[-2,-1.5,-1,-0.5,-0.2,-0.1, 0,0.1,0.2,0.5,1,1.5,2], linewidths=2) #np.arange(-15,-10,0.5)
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
    plt.savefig(cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/plots/"+name+str(h).zfill(2)+'_'+str(eh).zfill(2)+'_AMSR_big.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png)
    #plt.close()


def plot_doug_small(h, eh):

    dic = {}
    name = "ERA5_cores_propagation_WET_noERAfilter_AMSRE"  # "ERA5_composite_cores_AMSRE_w1_15k_minusMean"

    def coll(dic, h, eh, year):
        print(h)
        core = pkl.load(open(
            cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/core_txt/" + name + str(eh) + "UTCERA" + str(
                h).zfill(2) + '_' + str(year) + ".p", "rb"))
        for id, k in enumerate(core.keys()):
            try:
                dic[k] = dic[k] + core[k]
            except KeyError:
                dic[k] = core[k]

    for y in range(2006,2011):
            coll(dic, h, eh, y)

    extent = (dic['lsta0'].shape[1] - 1) / 2
    xlen = dic['lsta0'].shape[1]
    ylen = dic['lsta0'].shape[0]

    xv, yv = np.meshgrid(np.arange(ylen), np.arange(xlen))
    st = 30
    xquiv = xv[4::st, 4::st]
    yquiv = yv[4::st, 4::st]

    u = (dic['u925'] / dic['cnte'])[4::st, 4::st]
    v = (dic['v925'] / dic['cnte'])[4::st, 4::st]

    f = plt.figure(figsize=(10, 8), dpi=200)
    ax = f.add_subplot(221)

    plt.contourf((dic['lsta-1'] / dic['cnt-1']), cmap='RdBu',
                 levels=[-1.5,-1,-0.8, -0.6, -0.4, -0.2, -0.1, 0.1, 0.2, 0.4, 0.6, 0.8,1,1.5],
                 extend='both')  # -(rkernel2_sum / rcnt_sum)
    # plt.plot(extent, extent, 'bo')
    plt.colorbar(label='%')
    # pdb.set_trace()

    contours = plt.contour((dic['t'] / dic['cnte']), extend='both',levels=[ -0.8,-0.7, -0.6,-0.5,-0.4,-0.2,0, 0.2,0.4,0.5,0.6, 0.7, 0.8], cmap='PuOr_r', linewidths=1) # #, levels=np.arange(1,5, 0.5)
    plt.clabel(contours, inline=True, fontsize=11, fmt='%1.1f')

    contours = plt.contour((dic['v925']/ dic['cnte']), extend='both',colors='k', linewidths=4, levels=[-50,0,50])  #np.arange(-15,-10,0.5)
    #plt.clabel(contours, inline=True, fontsize=9, fmt='%1.1f')

    contours = plt.contour((dic['v925_orig']/ dic['cnte']), extend='both',colors='turquoise', linewidths=4, levels=[-50,0,50])  #np.arange(-15,-10,0.5)
    #plt.clabel(contours, inline=True, fontsize=9, fmt='%1.1f')

    plt.text(0.02,0.08, 'ITD 0-line', color='turquoise', fontsize=12, transform=ax.transAxes)
    plt.text(0.02, 0.03, 'ITD anomaly 0-line', color='k', fontsize=12, transform=ax.transAxes)

    # qu = ax.quiver(xquiv, yquiv, u, v, scale=15)
    # qk = plt.quiverkey(qu, 0.9, 0.02,1, '1 m s$^{-1}$',
    #                    labelpos='E', coordinates='figure')


    plt.clabel(contours, inline=True, fontsize=11, fmt='%1.0f')
    plt.plot(extent, extent, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, extent * 2, 9) - 100) * 6, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent * 2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.title(str(h).zfill(2)+'00UTC | '+str(np.max(dic['cnt0']))+' cores, SM night0, ERA5$_{noon}$ 925hPa T anomaly', fontsize=9)


    ax = f.add_subplot(222)

    plt.contourf((dic['lsta0'] / dic['cnt0']), cmap='RdBu',
                 levels=[-1.5,-1,-0.8, -0.6, -0.4, -0.2, -0.1, 0.1, 0.2, 0.4, 0.6, 0.8,1,1.5],
                 extend='both')  # -(rkernel2_sum / rcnt_sum)
    # plt.plot(extent, extent, 'bo')
    plt.colorbar(label='%')
    # pdb.set_trace()

    contours = plt.contour(((dic['div'])/ dic['cnte'])*100, extend='both',  cmap='PuOr', levels=[-0.8,-0.4,-0.2, 0.2,0.4,0.8]) # #, levels=np.arange(1,5, 0.5)
    qu = ax.quiver(xquiv, yquiv, u, v, scale=15)
    qk = plt.quiverkey(qu, 0.9, 0.02,1, '1 m s$^{-1}$',
                       labelpos='E', coordinates='figure')


    plt.clabel(contours, inline=True, fontsize=11, fmt='%1.1f')
    plt.plot(extent, extent, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, extent * 2, 9) - 100) * 6, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent * 2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.title('SM day0, Contours: 925hPa divergence', fontsize=9)




    ax1 = f.add_subplot(223)
    plt.contourf(((dic['q'])*1000/ dic['cnte']), extend='both',  cmap='RdBu',levels=np.arange(-1,1.1,0.1)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='g kg-1')
    contours = plt.contour((dic['shear'] / dic['cnte']), extend='both',levels=np.arange(-17,-12,0.5), cmap='jet') #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    # contours = plt.contour((dic['u650'] / dic['cnte']), extend='both', cmap='jet') #np.arange(-15,-10,0.5)
    # plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    plt.plot(extent, extent, 'bo')
    #qu = ax1.quiver(xquiv, yquiv, u, v, scale=50)
    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('Shading: 950hPa q anomaly, Contours: 600hPa-925hPa wind shear ', fontsize=9)


    ax1 = f.add_subplot(224)
    #   plt.contourf(((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r', vmin=-1.5, vmax=1.5) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.contourf(((dic['theta_e']) / dic['cnte']), extend='both', cmap='RdBu', levels=np.linspace(-2.3,2.3,16)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label=r'K')

    plt.plot(extent, extent, 'bo')
    #contours = plt.contour((dic['probmsg'] / dic['cntm']), extend='both', cmap='PuOr_r', levels=[-0.4,-0.3,-0.2,-0.1, 0,0.1,0.2,0.3,0.4]) #np.arange(-15,-10,0.5)
    contours = plt.contour((dic['probmsg'] / dic['cntm'])*100, extend='both', cmap='PuOr_r', levels=[-30,-20,-10,0,10,20,30])  # np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.0f')
    qu = ax.quiver(xquiv, yquiv, u, v, scale=15)
    qk = plt.quiverkey(qu, 0.9, 0.02,1, '1 m s$^{-1}$',
                       labelpos='E', coordinates='figure')

    ax1.set_xticklabels(np.array((np.linspace(0, extent * 2, 9) - 100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent * 2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title(r'$\Delta \theta_{e}$ anomaly, contours: MSG -40C cloud probability', fontsize=9)


    plt.tight_layout()
    plt.show()
    plt.savefig(
        cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/plots/" + name + str(h).zfill(2) + '_' + str(eh).zfill(
            2) + '_AMSR_small.png')




def plot_timeseries_small():

    x = len(list(range(-38, 4, 3)))
    y = 401

    # outticks = list(range(-30, 1, 5))
    # ranges = np.arange(-30,1,3)
    #
    # outticks = [12,17,22,3,8,13,18]
    #outticks = [6, 11, 16, 21, 2, 7, 12, 17]
    ranges = np.arange(-38, 4, 3)

    h = 17

    outdic = {}
    dummyfile = glob.glob(
        cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/core_txt/ERA5_cores_WET_SM0GT0.01-1GT0.01_noMeteosatFilter_AMSRE*.p")
    dummy = pkl.load(open(dummyfile[0], "rb"))

    #file = "figs/LSTA/corrected_LSTA/new/ERA5/core_txt/ERA5_cores_DRY_SM0LT3-1LT1.5_noMeteosatFilter_AMSRE"
    #file = "figs/LSTA/corrected_LSTA/new/ERA5/core_txt/ERA5_cores_WET_SM0GT0.01-1GT0.01_noMeteosatFilter_AMSRE"
    file = "figs/LSTA/corrected_LSTA/new/ERA5/core_txt/ERA5_cores_WET_TAG_noMeteosatFilter_AMSRE"
    for k in dummy.keys():
        outdic[k] = np.zeros((y, x))

    for ids, eh in enumerate(range(-38, 4, 3)):

        dic = {}

        def coll(dic, h, eh, year, file):
            print(h)


            core = pkl.load(open(
                cnst.network_data + file + str(eh) + "UTCERA" + str(h).zfill(2) + '_' + str(year) + ".p", "rb"))

            for id, k in enumerate(core.keys()):
                try:
                    dic[k] = dic[k] + core[k]
                except KeyError:
                    dic[k] = core[k]

        for y in range(2006, 2011):
            coll(dic, h, eh, y, file)

        for k in dic.keys():
            outdic[k][:, ids] = dic[k][:, 190:211].mean(axis=1)
    #
    # pos = np.isclose(outdic['v925'], np.zeros_like(outdic['v925']), atol=0.7)
    # #ipdb.set_trace()
    # ppos = np.where(pos)
    #
    # for pp in ppos:
    #     if pp[0] < 150:
    #         outdic['v925'][pp] = np.nan


    # plt.figure()
    # plt.pcolormesh(outdic['v925'], vmin=0)
    # plt.colorbar()

    if 'WET' in file:
        tag = 'WET'
    elif 'DRY' in file:
        tag = 'DRY'
    else:
        tag = 'ALL'


    #outdic['v925'][50:90,5:9] = np.nan


    print(outdic.keys())
    f = plt.figure(figsize=(6, 8), dpi=200)


    ax = f.add_subplot(211)
    plt.contourf(ranges, np.arange(401), (outdic['t']) / outdic['cnte'], extend='both', cmap='RdBu_r',
                levels=[-0.8, -0.7, -0.6, -0.5, -0.4, -0.2, -0.1, 0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8])
    plt.colorbar(label=r'K')
    plt.contour(ranges, np.arange(401), (outdic['t']) / outdic['cnte'], extend='both', colors='k', linewidths=0.1,
                 levels=[-0.8, -0.7, -0.6, -0.5, -0.4, -0.2, -0.1, 0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8])

    contours = plt.contour(ranges, np.arange(401), (outdic['u650'] / outdic['cnte']), extend='both', colors='k',
                           levels=np.arange(-1,0,0.2), linewidths=1)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')


    contours = plt.contour(ranges, np.arange(401),(outdic['v925']/ outdic['cnte']), extend='both',colors='k', linewidths=5, levels=[-50,0,50])  #np.arange(-15,-10,0.5)
    contours = plt.contour(ranges, np.arange(401),(outdic['v925_orig']/ outdic['cnte']), extend='both',colors='k', linewidths=5, levels=[-50,0,50])

    contours = plt.contour(ranges, np.arange(401),(outdic['v925']/ outdic['cnte']), extend='both',colors='k', linewidths=3, levels=[-50,0,50])  #np.arange(-15,-10,0.5)
    contours = plt.contour(ranges, np.arange(401),(outdic['v925_orig']/ outdic['cnte']), extend='both',colors='r', linewidths=3, levels=[-50,0,50])


    contours = plt.contour(ranges, np.arange(401),(outdic['v650_orig']/ outdic['cnte']), extend='both',colors=['r','r','white','b', 'b'], linewidths=0.5, levels=[-2,-1,0,0.25,0.5])  #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
   # contours = plt.contour(ranges, np.arange(401),(outdic['v650_orig']/ outdic['cnte']), extend='both',colors='b', linewidths=0.5, levels=[-50,0,50])

    plt.axvline(x=-5, color='slategrey')
    plt.axvline(x=-29, color='slategrey')
    plt.axhline(y=200, color='slategrey')
    plt.plot(-5, 200, 'ko')
    plt.plot(0, 200, 'ro')

    plt.text(0.02,0.1, 'ITD 0-line', color='red', fontsize=14, transform=ax.transAxes)
    plt.text(0.02, 0.03, 'v-wind anomaly 0-line', color='k', fontsize=14, transform=ax.transAxes)


    ax.set_yticklabels(np.array((np.linspace(0, 200 * 2, 9) - 200) * 3, dtype=int))
    # ax.set_xticklabels(outticks)

    plt.title('Shading:t-anomaly, contours: 650hpa u-wind anomaly')
    #plt.hlines(200, xmin=ranges[0], xmax=ranges[-1], linewidth=1)
    plt.xlabel('Hour relative to t0 [1700UTC]')
    plt.ylabel('North-South distance from core (km)')

    ax = f.add_subplot(212)


    plt.contourf(ranges, np.arange(401), (outdic['q']) * 1000 / outdic['cnte'], levels=[-0.8, -0.7, -0.6, -0.5, -0.4, -0.2, -0.1, 0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8],
                 cmap='RdBu', extend='both')
    plt.colorbar(label=r'g kg$^{-1}$')

    plt.contour(ranges, np.arange(401), (outdic['q']) * 1000 / outdic['cnte'], levels=[-0.8, -0.7, -0.6, -0.5, -0.4, -0.2, -0.1, 0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8],
     colors = 'k', linewidths = 0.1)

    contours = plt.contour(ranges, np.arange(401), outdic['theta_e'] / outdic['cnte'],colors='white', levels=np.arange(1,3,0.5), linewidths=1.5)

    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')


    contours = plt.contour(ranges, np.arange(401), (outdic['v925'] / outdic['cnte']), extend='both', colors='k',
                           linewidths=5, levels=[-50, 0, 50])  # np.arange(-15,-10,0.5)
    contours = plt.contour(ranges, np.arange(401), (outdic['v925_orig'] / outdic['cnte']), extend='both', colors='k',
                           linewidths=5, levels=[-50, 0, 50])

    contours = plt.contour(ranges, np.arange(401), (outdic['v925'] / outdic['cnte']), extend='both', colors='k',
                           linewidths=3, levels=[-50, 0, 50])  # np.arange(-15,-10,0.5)
    contours = plt.contour(ranges, np.arange(401), (outdic['v925_orig'] / outdic['cnte']), extend='both', colors='r',
                           linewidths=3, levels=[-50, 0, 50])

    plt.axvline(x=-5, color='slategrey')
    plt.axvline(x=-29, color='slategrey')
    plt.axhline(y=200, color='slategrey')
    plt.plot(-5, 200, 'ko')
    plt.plot(0, 200, 'ro')

    ax.set_yticklabels(np.array((np.linspace(0, 200 * 2, 9) - 200) * 3, dtype=int))
    # ax.set_xticklabels(outticks)

    plt.title(r'Shading:q-anomaly, contours: $\Delta \theta_{e}$ anomaly')
    # plt.hlines(200, xmin=ranges[0], xmax=ranges[-1], linewidth=1)
    plt.xlabel('Hour relative to t0 [1700UTC]')
    plt.ylabel('North-South distance from core (km)')

    plt.tight_layout()
    #plt.show()
    plt.savefig(cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/plots/ERA5_amsre_"+str(h).zfill(2)+'_timeseries_SMALL_'+tag+'.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png)
    plt.close()


def plot_timeseries_diff():

    x = len(list(range(-38, 4, 3)))
    y = 401

    # outticks = list(range(-30, 1, 5))
    # ranges = np.arange(-30,1,3)
    #
    # outticks = [12,17,22,3,8,13,18]
    #outticks = [6, 11, 16, 21, 2, 7, 12, 17]
    ranges = np.arange(-38, 4, 3)

    h = 17

    outdic1 = {}
    outdic2 = {}
    dummyfile = glob.glob(
        cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/core_txt/ERA5_cores_WET_SM0GT0.01-1GT0.01_noMeteosatFilter_AMSRE*.p")
    dummy = pkl.load(open(dummyfile[0], "rb"))

    file1 = "figs/LSTA/corrected_LSTA/new/ERA5/core_txt/ERA5_cores_DRY_SM0LT3-1LT1.5_noMeteosatFilter_AMSRE"
    file2 = "figs/LSTA/corrected_LSTA/new/ERA5/core_txt/ERA5_cores_WET_SM0GT0.01-1GT0.01_noMeteosatFilter_AMSRE"
    #file2 = "figs/LSTA/corrected_LSTA/new/ERA5/core_txt/ERA5_cores_WET_TAG_noMeteosatFilter_AMSRE"
    files = [file1, file2]
    outdics = [outdic1, outdic2]

    for file, outdic in zip(files,outdics):

        for k in dummy.keys():
            outdic[k] = np.zeros((y, x))

        for ids, eh in enumerate(range(-38, 4, 3)):

            dic = {}

            def coll(dic, h, eh, year, file):
                print(h)

                core = pkl.load(open(
                    cnst.network_data + file + str(eh) + "UTCERA" + str(h).zfill(2) + '_' + str(year) + ".p", "rb"))

                for id, k in enumerate(core.keys()):
                    try:
                        dic[k] = dic[k] + core[k]
                    except KeyError:
                        dic[k] = core[k]

            for yy in range(2006, 2011):
                coll(dic, h, eh, yy, file)

            for k in dic.keys():
                try:
                    outdic[k][:, ids] = dic[k][:, 190:211].mean(axis=1)
                except ValueError:
                    ipdb.set_trace()

    diff = {}
    for k in outdics[0].keys():

        if 'cnt' in k:
            continue

        if 'lsta0' in k:
            diff[k] = (outdics[0][k] / outdics[0]['cnt0']) - (outdics[1][k] / outdics[1]['cnt0'])
        elif 'lsta-1' in k:
            diff[k] = (outdics[0][k] / outdics[0]['cnt-1']) - (outdics[1][k] / outdics[1]['cnt-1'])
        elif 'lsta-2' in k:
            diff[k] = (outdics[0][k] / outdics[0]['cnt-2']) - (outdics[1][k] / outdics[1]['cnt-2'])
        elif 'lsta-3' in k:
            diff[k] = (outdics[0][k] / outdics[0]['cnt-3']) - (outdics[1][k] / outdics[1]['cnt-3'])
        elif 'msg' in k:
            diff[k] = (outdics[0][k] / outdics[0]['cntm']) - (outdics[1][k] / outdics[1]['cntm'])
        elif 'probc' in k:
            diff[k] = (outdics[0][k] / outdics[0]['cntc']) - (outdics[1][k] / outdics[1]['cntc'])
        else:
            diff[k] = (outdics[0][k] / outdics[0]['cnte']) - (outdics[1][k] / outdics[1]['cnte'])

    def groupedAvg(myArray, N=2):
        result = np.cumsum(myArray, 0)[N - 1::N, :] / float(N)
        result[1:,:] = result[1:, :] - result[:-1, :]
        return result

    tcoarse = groupedAvg(diff['t'], N=4)
    plt.figure()
    plt.pcolormesh(tcoarse)
    print(tcoarse.shape)


    f = plt.figure(figsize=(6, 8), dpi=200)




    ax = f.add_subplot(211)
    # plt.contourf(ranges, np.arange(401), (diff['t']), extend='both', cmap='RdBu_r',
    #             levels=[-0.8, -0.7, -0.6, -0.5, -0.4, -0.2, -0.1, 0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8])
    # plt.colorbar(label=r'K')
    # # plt.contour(ranges, np.arange(401), (diff['t']), extend='both', colors='k', linewidths=0.1,
    #              levels=[-0.8, -0.7, -0.6, -0.5, -0.4, -0.2, -0.1, 0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8])
    #ipdb.set_trace()
    plt.contourf(ranges, np.arange(100), tcoarse, extend='both', cmap='RdBu_r',
                levels=[-0.8, -0.7, -0.6, -0.5, -0.4, -0.2, -0.1, 0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8])
    plt.contour(ranges, np.arange(100), tcoarse, extend='both', colors='k', linewidths=0.1,
                 levels=[-0.8, -0.7, -0.6, -0.5, -0.4, -0.2, -0.1, 0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8])

    contours = plt.contour(ranges, np.arange(102), (diff['u650']), extend='both', colors='k',
                           levels=np.arange(-2.5,-1.6,0.5), linewidths=1)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')


    #contours = plt.contour(ranges, np.arange(401),(diff['v925']), extend='both',colors='k', linewidths=5, levels=[-50,0,50])  #np.arange(-15,-10,0.5)
    contours = plt.contour(ranges, np.arange(401),(outdics[0]['v925_orig']/outdics[0]['cnte']), extend='both',colors='r', linewidths=5, levels=[-50,0.06,50])

    #contours = plt.contour(ranges, np.arange(401),(diff['v925']), extend='both',colors='k', linewidths=3, levels=[-50,0,50])  #np.arange(-15,-10,0.5)
    contours = plt.contour(ranges, np.arange(401),(outdics[1]['v925_orig']/outdics[1]['cnte']), extend='both',colors='b', linewidths=5, levels=[-50,0.06,50])



   # contours = plt.contour(ranges, np.arange(401),(outdic['v650_orig']/ outdic['cnte']), extend='both',colors='b', linewidths=0.5, levels=[-50,0,50])

    plt.axvline(x=-5, color='slategrey')
    plt.axvline(x=-29, color='slategrey')
    plt.axhline(y=200, color='slategrey')
    plt.plot(-5, 200, 'ko')
    plt.plot(0, 200, 'ro')

    plt.text(0.02,0.1, 'ITD 0-line', color='red', fontsize=14, transform=ax.transAxes)
    plt.text(0.02, 0.03, 'v-wind anomaly 0-line', color='k', fontsize=14, transform=ax.transAxes)


    ax.set_yticklabels(np.array((np.linspace(0, 200 * 2, 9) - 200) * 3, dtype=int))
    # ax.set_xticklabels(outticks)

    plt.title('Shading:t-anomaly, contours: 650hpa u-wind anomaly')
    #plt.hlines(200, xmin=ranges[0], xmax=ranges[-1], linewidth=1)
    plt.xlabel('Hour relative to t0 [1700UTC]')
    plt.ylabel('North-South distance from core (km)')

    ax = f.add_subplot(212)


    plt.contourf(ranges, np.arange(401), (diff['q']) * 1000, levels=[-0.8, -0.7, -0.6, -0.5, -0.4, -0.2, -0.1, 0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8],
                 cmap='RdBu', extend='both')
    plt.colorbar(label=r'g kg$^{-1}$')

    plt.contour(ranges, np.arange(401), (diff['q']) * 1000, levels=[-0.8, -0.7, -0.6, -0.5, -0.4, -0.2, -0.1, 0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8],
     colors = 'k', linewidths = 0.1)

    contours = plt.contour(ranges, np.arange(401), diff['theta_e'],colors='white', levels=np.arange(1,3,0.5), linewidths=1.5)

    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')


    contours = plt.contour(ranges, np.arange(401),(outdics[0]['v650_orig']/outdics[0]['cnte']), extend='both',colors='r', linewidths=0.5, levels=[-2,-1,0,0.25,0.5])  #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')


    contours = plt.contour(ranges, np.arange(401),(outdics[1]['v650_orig']/outdics[0]['cnte']), extend='both',colors='b', linewidths=0.5, levels=[-2,-1,0,0.25,0.5])  #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')

    plt.axvline(x=-5, color='slategrey')
    plt.axvline(x=-29, color='slategrey')
    plt.axhline(y=200, color='slategrey')
    plt.plot(-5, 200, 'ko')
    plt.plot(0, 200, 'ro')

    ax.set_yticklabels(np.array((np.linspace(0, 200 * 2, 9) - 200) * 3, dtype=int))
    # ax.set_xticklabels(outticks)

    plt.title(r'Shading:q-anomaly, contours: $\Delta \theta_{e}$ anomaly')
    # plt.hlines(200, xmin=ranges[0], xmax=ranges[-1], linewidth=1)
    plt.xlabel('Hour relative to t0 [1700UTC]')
    plt.ylabel('North-South distance from core (km)')

    plt.tight_layout()
    #plt.show()
    plt.savefig(cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/plots/ERA5_amsre_"+str(h).zfill(2)+'_timeseries_SMALL_DIFF_STRONG.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png)
    plt.close()


def plot_timeseries_diff_coarse():

    x = len(list(range(-38, 4, 3)))
    y = 401

    # outticks = list(range(-30, 1, 5))
    # ranges = np.arange(-30,1,3)
    #
    # outticks = [12,17,22,3,8,13,18]
    #outticks = [6, 11, 16, 21, 2, 7, 12, 17]
    ranges = np.arange(-38, 4, 3)

    h = 17

    outdic1 = {}
    outdic2 = {}
    dummyfile = glob.glob(
        cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/core_txt/ERA5_cores_WET_SM0GT0.01-1GT0.01_noMeteosatFilter_AMSRE*.p")
    dummy = pkl.load(open(dummyfile[0], "rb"))

    file1 = "figs/LSTA/corrected_LSTA/new/ERA5/core_txt/ERA5_cores_DRY_SM0LT3-1LT1.5_noMeteosatFilter_AMSRE"
    file2 = "figs/LSTA/corrected_LSTA/new/ERA5/core_txt/ERA5_cores_WET_SM0GT0.01-1GT0.01_noMeteosatFilter_AMSRE"
    #file2 = "figs/LSTA/corrected_LSTA/new/ERA5/core_txt/ERA5_cores_WET_TAG_noMeteosatFilter_AMSRE"
    files = [file1, file2]
    outdics = [outdic1, outdic2]

    for file, outdic in zip(files,outdics):

        for k in dummy.keys():
            outdic[k] = np.zeros((y, x))

        for ids, eh in enumerate(range(-38, 4, 3)):

            dic = {}

            def coll(dic, h, eh, year, file):
                print(h)

                core = pkl.load(open(
                    cnst.network_data + file + str(eh) + "UTCERA" + str(h).zfill(2) + '_' + str(year) + ".p", "rb"))

                for id, k in enumerate(core.keys()):
                    try:
                        dic[k] = dic[k] + core[k]
                    except KeyError:
                        dic[k] = core[k]

            for yy in range(2006, 2011):
                coll(dic, h, eh, yy, file)

            for k in dic.keys():
                try:
                    outdic[k][:, ids] = dic[k][:, 190:211].mean(axis=1)
                except ValueError:
                    ipdb.set_trace()

    def groupedAvg(myArray, N=2):
        result = np.cumsum(myArray, 0)[N - 1::N, :] / float(N)
        result[1:,:] = result[1:, :] - result[:-1, :]
        return result

    diff = {}
    for k in outdics[0].keys():

        if 'cnt' in k:
            continue

        if 'lsta0' in k:
            diff[k] = groupedAvg((outdics[0][k] / outdics[0]['cnt0']) - (outdics[1][k] / outdics[1]['cnt0']), N=4)
        elif 'lsta-1' in k:
            diff[k] = groupedAvg((outdics[0][k] / outdics[0]['cnt-1']) - (outdics[1][k] / outdics[1]['cnt-1']), N=4)
        elif 'lsta-2' in k:
            diff[k] = groupedAvg((outdics[0][k] / outdics[0]['cnt-2']) - (outdics[1][k] / outdics[1]['cnt-2']), N=4)
        elif 'lsta-3' in k:
            diff[k] = groupedAvg((outdics[0][k] / outdics[0]['cnt-3']) - (outdics[1][k] / outdics[1]['cnt-3']), N=4)
        elif 'msg' in k:
            diff[k] = groupedAvg((outdics[0][k] / outdics[0]['cntm']) - (outdics[1][k] / outdics[1]['cntm']), N=4)
        elif 'probc' in k:
            diff[k] = groupedAvg((outdics[0][k] / outdics[0]['cntc']) - (outdics[1][k] / outdics[1]['cntc']), N=4)
        elif 'v925' in k:
            diff[k+'_0'] = groupedAvg((outdics[0][k] / outdics[0]['cnte']) , N=4)
            diff[k + '_1'] = groupedAvg((outdics[1][k] / outdics[1]['cnte']), N=4)
        # elif 'v650' in k:
        #     diff[k+'_0'] = groupedAvg((outdics[0][k] / outdics[0]['cnte']) , N=4)
        #     diff[k + '_1'] = groupedAvg((outdics[1][k] / outdics[1]['cnte']), N=4)
        else:
            diff[k] = groupedAvg((outdics[0][k] / outdics[0]['cnte']) - (outdics[1][k] / outdics[1]['cnte']), N=4)


    print(diff.keys())

    f = plt.figure(figsize=(6, 8), dpi=200)


    yax = np.arange(100)
    ax = f.add_subplot(211)

    #ipdb.set_trace()
    plt.contourf(ranges, yax, (diff['t']), extend='both', cmap='RdBu_r',
                levels=[-0.8, -0.7, -0.6, -0.5, -0.4, -0.2, -0.1, 0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8])
    plt.colorbar(label=r'K')
    plt.contour(ranges, yax, (diff['t']), extend='both', colors='k', linewidths=0.1,
                 levels=[-0.8, -0.7, -0.6, -0.5, -0.4, -0.2, -0.1, 0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8])
    #ipdb.set_trace()


    contours = plt.contour(ranges, yax, (diff['u650']), extend='both', colors='k',
                           levels=np.arange(-2.5,-1,0.5), linewidths=1)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')


    #contours = plt.contour(ranges, yax,(diff['v925']), extend='both',colors='k', linewidths=5, levels=[-50,0,50])  #np.arange(-15,-10,0.5)
    contours = plt.contour(ranges, yax,(diff['v925_orig_0']), extend='both',colors='r', linewidths=5, levels=[-50,0.06,50])

    #contours = plt.contour(ranges, yax,(diff['v925']), extend='both',colors='k', linewidths=3, levels=[-50,0,50])  #np.arange(-15,-10,0.5)
    contours = plt.contour(ranges, yax,(diff['v925_orig_1']), extend='both',colors='b', linewidths=5, levels=[-50,0.06,50])



   # contours = plt.contour(ranges, yax,(outdic['v650_orig']/ outdic['cnte']), extend='both',colors='b', linewidths=0.5, levels=[-50,0,50])

    plt.axvline(x=-5, color='slategrey')
    plt.axvline(x=-29, color='slategrey')
    plt.axhline(y=50, color='slategrey')
    plt.plot(-5, 50, 'ko')
    plt.plot(0, 50, 'ro')

    plt.text(0.02,0.1, 'ITD 0-line DRY', color='r', fontsize=14, transform=ax.transAxes)
    plt.text(0.02, 0.03, 'ITD 0-line WET', color='b', fontsize=14, transform=ax.transAxes)

    ax.set_yticklabels(np.array((np.linspace(0, 100, 6) - 50) * 12, dtype=int))
    # ax.set_xticklabels(outticks)

    plt.title('Shading:t-anomaly, contours: 650hpa u-wind anomaly')
    #plt.hlines(200, xmin=ranges[0], xmax=ranges[-1], linewidth=1)
    plt.xlabel('Hour relative to t0 [1700UTC]')
    plt.ylabel('North-South distance from core (km)')

    ax = f.add_subplot(212)


    plt.contourf(ranges, yax, (diff['q']) * 1000, levels=[-0.8, -0.7, -0.6, -0.5, -0.4, -0.2, -0.1, 0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8],
                 cmap='RdBu', extend='both')
    plt.colorbar(label=r'g kg$^{-1}$')

    plt.contour(ranges, yax, (diff['q']) * 1000, levels=[-0.8, -0.7, -0.6, -0.5, -0.4, -0.2, -0.1, 0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8],
     colors = 'k', linewidths = 0.1)

    contours = plt.contour(ranges, yax, diff['theta_e'],colors='r', levels=np.arange(1,3,0.5), linewidths=1.5)

    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')

    ax.streamplot(ranges, yax, diff['u650_orig'] * 0.1, diff['v650_orig'], density=[0.5, 1], linewidth=0.5, color='k')


    # contours = plt.contour(ranges, yax,(diff['v650_orig_1']), extend='both',colors='b', linewidths=0.5, levels=[-2,-1,0,0.25,0.5])  #np.arange(-15,-10,0.5)
    # plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')

    plt.axvline(x=-5, color='slategrey')
    plt.axvline(x=-29, color='slategrey')
    plt.axhline(y=50, color='slategrey')
    plt.plot(-5, 50, 'ko')
    plt.plot(0, 50, 'ro')

    ax.set_yticklabels(np.array((np.linspace(0, 100, 6) - 50) * 12, dtype=int))
    # ax.set_xticklabels(outticks)

    plt.title(r'Shading:q-anomaly, contours: $\Delta \theta_{e}$ anomaly')
    # plt.hlines(200, xmin=ranges[0], xmax=ranges[-1], linewidth=1)
    plt.xlabel('Hour relative to t0 [1700UTC]')
    plt.ylabel('North-South distance from core (km)')

    plt.tight_layout()
    #plt.show()
    plt.savefig(cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/plots/ERA5_amsre_"+str(h).zfill(2)+'_timeseries_SMALL_DIFF_coarse.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png)
    plt.close()


def plot_timeseries_coarse():

    x = len(list(range(-38, 4, 3)))
    y = 401

    # outticks = list(range(-30, 1, 5))
    # ranges = np.arange(-30,1,3)
    #
    # outticks = [12,17,22,3,8,13,18]
    #outticks = [6, 11, 16, 21, 2, 7, 12, 17]
    ranges = np.arange(-38, 4, 3)

    h = 17

    outdic = {}
    dummyfile = glob.glob(
        cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/core_txt/ERA5_cores_WET_SM0GT0.01-1GT0.01_noMeteosatFilter_AMSRE*.p")
    dummy = pkl.load(open(dummyfile[0], "rb"))

    #file = "figs/LSTA/corrected_LSTA/new/ERA5/core_txt/ERA5_cores_DRY_SM0LT3-1LT1.5_noMeteosatFilter_AMSRE"
    file = "figs/LSTA/corrected_LSTA/new/ERA5/core_txt/ERA5_cores_WET_SM0GT0.01-1GT0.01_noMeteosatFilter_AMSRE"
    #file = "figs/LSTA/corrected_LSTA/new/ERA5/core_txt/ERA5_cores_WET_TAG_noMeteosatFilter_AMSRE"
    for k in dummy.keys():
        outdic[k] = np.zeros((y, x))

    for ids, eh in enumerate(range(-38, 4, 3)):

        dic = {}

        def coll(dic, h, eh, year, file):
            print(h)


            core = pkl.load(open(
                cnst.network_data + file + str(eh) + "UTCERA" + str(h).zfill(2) + '_' + str(year) + ".p", "rb"))

            for id, k in enumerate(core.keys()):
                try:
                    dic[k] = dic[k] + core[k]
                except KeyError:
                    dic[k] = core[k]

        for y in range(2006, 2011):
            coll(dic, h, eh, y, file)

        for k in dic.keys():
            outdic[k][:, ids] = dic[k][:, 190:211].mean(axis=1)
    #
    if 'WET' in file:
        tag = 'WET'
    elif 'DRY' in file:
        tag = 'DRY'
    else:
        tag = 'ALL'

    def groupedAvg(myArray, N=2):
        result = np.cumsum(myArray, 0)[N - 1::N, :] / float(N)
        result[1:,:] = result[1:, :] - result[:-1, :]
        return result

    diff = {}
    for k in outdic.keys():

        if 'cnt' in k:
            continue

        if 'lsta0' in k:
            diff[k] = groupedAvg((outdic[k] / outdic['cnt0']) , N=4)
        elif 'lsta-1' in k:
            diff[k] = groupedAvg((outdic[k] / outdic['cnt-1']) , N=4)
        elif 'lsta-2' in k:
            diff[k] = groupedAvg((outdic[k] / outdic['cnt-2']), N=4)
        elif 'lsta-3' in k:
            diff[k] = groupedAvg((outdic[k] / outdic['cnt-3']) , N=4)
        elif 'msg' in k:
            diff[k] = groupedAvg((outdic[k] / outdic['cntm']) , N=4)
        elif 'probc' in k:
            diff[k] = groupedAvg((outdic[k] / outdic['cntc']) , N=4)
        # elif 'v925' in k:
        #     diff[k+'_0'] = groupedAvg((outdic[k] / outdic['cnte']) , N=4)
        # elif 'v650' in k:
        #     diff[k+'_0'] = groupedAvg((outdic[k] / outdic['cnte']) , N=4)
        else:
            diff[k] = groupedAvg((outdic[k] / outdic['cnte']) , N=4)


    print(diff.keys())

    f = plt.figure(figsize=(6, 8), dpi=200)


    yax = np.arange(100)
    ax = f.add_subplot(211)

    #ipdb.set_trace()
    plt.contourf(ranges, yax, (diff['t']), extend='both', cmap='RdBu_r',
                levels=[-0.8, -0.7, -0.6, -0.5, -0.4, -0.2, -0.1, 0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8])
    plt.colorbar(label=r'K')
    plt.contour(ranges, yax, (diff['t']), extend='both', colors='k', linewidths=0.1,
                 levels=[-0.8, -0.7, -0.6, -0.5, -0.4, -0.2, -0.1, 0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8])
    #ipdb.set_trace()


    contours = plt.contour(ranges, yax, (diff['u650']), extend='both', colors='k',
                           levels=np.arange(-2.5,0,0.5), linewidths=1)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')


    #contours = plt.contour(ranges, yax,(diff['v925']), extend='both',colors='k', linewidths=5, levels=[-50,0,50])  #np.arange(-15,-10,0.5)
    contours = plt.contour(ranges, yax,(diff['v925_orig']), extend='both',colors='r', linewidths=5, levels=[-50,0.06,50])


   # contours = plt.contour(ranges, yax,(outdic['v650_orig']/ outdic['cnte']), extend='both',colors='b', linewidths=0.5, levels=[-50,0,50])

    plt.axvline(x=-5, color='slategrey')
    plt.axvline(x=-29, color='slategrey')
    plt.axhline(y=50, color='slategrey')
    plt.plot(-5, 50, 'ko')
    plt.plot(0, 50, 'ro')

    plt.text(0.02,0.1, 'ITD 0-line', color='red', fontsize=14, transform=ax.transAxes)
    #plt.text(0.02, 0.03, 'v-wind anomaly 0-line', color='k', fontsize=14, transform=ax.transAxes)

    ax.set_yticklabels(np.array((np.linspace(0, 100, 6) - 50) * 12, dtype=int))
    # ax.set_xticklabels(outticks)

    plt.title('Shading:t-anomaly, contours: 650hpa u-wind anomaly')
    #plt.hlines(200, xmin=ranges[0], xmax=ranges[-1], linewidth=1)
    plt.xlabel('Hour relative to t0 [1700UTC]')
    plt.ylabel('North-South distance from core (km)')

    ax = f.add_subplot(212)


    plt.contourf(ranges, yax, (diff['q']) * 1000, levels=[-0.8, -0.7, -0.6, -0.5, -0.4, -0.2, -0.1, 0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8],
                 cmap='RdBu', extend='both')
    plt.colorbar(label=r'g kg$^{-1}$')

    plt.contour(ranges, yax, (diff['q']) * 1000, levels=[-0.8, -0.7, -0.6, -0.5, -0.4, -0.2, -0.1, 0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8],
     colors = 'k', linewidths = 0.1)

    contours = plt.contour(ranges, yax, diff['theta_e'],colors='r', levels=np.arange(1,3,0.5), linewidths=1.5)

    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')


    #contours = plt.contour(ranges, yax,(diff['v650_orig_0']), extend='both',colors='r', linewidths=0.5, levels=[-2,-1,0,0.25,0.5])  #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    #xv, yv = np.meshgrid(np.arange(ylen), np.arange(xlen))
    ax.streamplot(ranges, yax,diff['u650_orig']*0.01, diff['v650_orig'], density=[0.5, 1], linewidth=0.5, color='k')  #*0.01

    plt.axvline(x=-5, color='slategrey')
    plt.axvline(x=-29, color='slategrey')
    plt.axhline(y=50, color='slategrey')
    plt.plot(-5, 50, 'ko')
    plt.plot(0, 50, 'ro')

    ax.set_yticklabels(np.array((np.linspace(0, 100, 6) - 50) *12 , dtype=int))
    # ax.set_xticklabels(outticks)

    plt.title(r'Shading:q-anomaly, contours: $\Delta \theta_{e}$ anomaly')
    # plt.hlines(200, xmin=ranges[0], xmax=ranges[-1], linewidth=1)
    plt.xlabel('Hour relative to t0 [1700UTC]')
    plt.ylabel('North-South distance from core (km)')

    plt.tight_layout()
    #plt.show()
    plt.savefig(cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/plots/ERA5_amsre_"+str(h).zfill(2)+'_timeseries_SMALL_'+tag+'_coarse.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png)
    plt.close()




def plot_all():

    hours = [16,17,18,19,20,21,22,23,0,1,2,3,4,5,6,7]

    for h in hours:
        plot_doug_all(h)
