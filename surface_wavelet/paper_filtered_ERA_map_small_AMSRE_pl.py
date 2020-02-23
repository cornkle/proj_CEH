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
import metpy
from metpy import calc
from metpy.units import units

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
    all = np.arange(-38, 4, 3)
    #all = np.arange(-50,31,3)
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

    key = '2hOverlap'
    msgopen = pd.read_csv(cnst.network_data + 'figs/LSTA/corrected_LSTA/new/ERA5/core_txt/cores_gt15000km2_table_AMSRE_LSTA_tracking_new_2hOverlap_'+str(h)+'.csv', na_values=-999)
    hour = h
    msg = pd.DataFrame.from_dict(msgopen)
    msg['eh'] = eh
    msg['refhour'] = h

    msg['date'] = pd.to_datetime(msg[['year','month','day']])
    print('Start core number ', len(msg))

    msg = msg[ (msg['dtime']<=2) & (msg['LSTAslotfrac'] >= 0.03) ]

    #msg = msg[(msg['dtime'] <= 2) & (msg['SMmean0']>0.01) & (msg['SMmean-1']>0.8) & (msg['LSTAslotfrac'] >= 0.03) ] #& (msg['LSTAslotfrac'] >= 0.1)
    #msg = msg[(msg['dtime'] <= 2) & (msg['SMmean0'] < -6) & (msg['SMmean-1'] < -3.3) & (msg['LSTAslotfrac'] >= 0.03)] # slotfrac is an okay filter to make sure it doesnt rain on itself

    msgin = msg[(msg['lat']>9) & (msg['lat']<19) & (msg['topo']<450) ] #& (msgopen['ERAqmean']>14.8) & (msgopen['ERAqmean']<16.5)

    print('Number of cores', len(msgin))
    #ipdb.set_trace()

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

        pkl.dump(dic, open(cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/core_txt/ERA5_pressureLevels_WE_"+key+"_AMSRE_SMALL_"+str(eh) + "UTCERA"+str(hour).zfill(2)+'_'+str(year)+".p", "wb"))
        del dic
        print('Dumped file')



def cut_kernel(xpos, ypos, era, dist):

    vdic = {}

    for d in era.data_vars:

        var = ua.cut_kernel_3d(era[d].values,xpos, ypos,dist)
        #var = np.mean(var[:,:, dist-1:dist+2], axis=2)
        var = np.mean(var[:, dist - 1:dist + 2, :], axis=1)
        vdic[d] = var

    cntera = np.zeros_like(vdic[list(vdic.keys())[0]])
    cntera[np.isfinite(vdic[list(vdic.keys())[0]])] = 1

    return vdic,cntera




def get_previous_hours(storm_date, ehour, refhour):


    date = storm_date.replace(hour=refhour)
    cm = xr.Dataset()

    if ehour > 0:
        edate = date + pd.Timedelta(str(ehour) + ' hours')
    else:
        edate = date - pd.Timedelta(str(np.abs(ehour)) + ' hours')
    #edate = edate.replace(hour=ehour)


    t1 = edate
    file = cnst.ERA5


    try:
        cmp = xr.open_dataset(file + 'hourly/pressure_levels/ERA5_'+str(edate.year)+'_' + str(edate.month).zfill(2) + '_pl.nc')
        cmp = u_darrays.flip_lat(cmp)
    except:
        return None


    pl_clim = xr.open_dataset(file + 'monthly/synop_selfmade/CLIM_2006-2010/ERA5_2006-2010_CLIM_'+str(edate.month).zfill(2)+'-'+str(edate.day).zfill(2)+'-'+str(edate.hour).zfill(2)+'_pl.nc').load()
    pl_clim = u_darrays.flip_lat(pl_clim)
    try:
        pl_clim = pl_clim.rename({'lat':'latitude', 'lon':'longitude'})
    except:
        pass

    ## latitude in surface is already flipped, not for pressure levels though... ?!

    cmp = cmp.sel(longitude=slice(-12.5, 12.5), latitude=slice(7, 23))
    pl_clim = pl_clim.sel(longitude=slice(-12.5, 12.5), latitude=slice(7, 23))

    cmm = cmp.sel(time=t1)
    pl_clim = pl_clim.squeeze()

    uwind = cmm['u']
    vwind= cmm['v']
    vwind_clim = pl_clim['v']
    uwind_clim = pl_clim['u']

    div = cmm['d']
    divclim = pl_clim['d']

    thetae = u_met.theta_e(pl_clim.level.values, cmm['t'].values - 273.15,
                            cmm['q'])

    thetae_clim = u_met.theta_e(pl_clim.level.values, pl_clim['t'].values - 273.15,
                                pl_clim['q'])

    theta = u_met.theta(pl_clim.level.values, cmm['t'] - 273.15)
    theta_clim = u_met.theta(pl_clim.level.values, pl_clim['t']- 273.15)

    z = cmm['z']/9.81
    zclim = pl_clim['z']/9.81

    q = cmm['q']
    qclim = pl_clim['q']


    cm['u'] = uwind
    cm['v'] = vwind

    cm['uclim'] = uwind_clim
    cm['vclim'] = vwind_clim


    cm['div'] = div *1000
    cm['divclim'] = divclim * 1000

    cm['q'] = q
    cm['qclim'] = qclim

    cm['theta'] = theta
    cm['theta_clim'] = theta_clim
    cm['thetae'] = thetae
    cm['thetae_clim'] = thetae_clim

    cm['z'] = z
    cm['zclim'] = zclim

    del pl_clim
    del cmm
    return cm


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
    next_lsta_date = lsta_date + dayd

    fdate = str(lsta_date.year) + str(lsta_date.month).zfill(2) + str(lsta_date.day).zfill(2)
    pfdate = str(prev_lsta_date.year) + str(prev_lsta_date.month).zfill(2) + str(prev_lsta_date.day).zfill(2)

    topo = xr.open_dataset(cnst.WA_TOPO_3KM)
    topo = topo.sel(lon=slice(-12, 12), lat=slice(7, 23))
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

    smcnt = []
    dist = 22 # (22*27 = 594)

    for sid , sp in enumerate(smpath):

        try:
            lsta = xr.open_dataset(sp)
        except OSError:
                return None
        print('Doing '+ sp)

        lsta = lsta.sel(lon=slice(-12.5, 12.5), lat=slice(7, 23))

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
        del lsta

    if len(smlist)!=4:
        return None

    cnte_sum = np.zeros((17, dist * 2 + 1))
    edic = {}
    try:
        probs = get_previous_hours(storm_date, eh, hour)
    except:
        print('ERA failed')
        return
    print('Era5 collect')

    del topo

    counter=0

    for dids, dit in df.iterrows():

        lat = dit.lat
        lon = dit.lon

        try:
            point = probs.sel(latitude=lat, longitude=lon, method='nearest', tolerance=0.26)
        except KeyError:
            print('Tolerance error')
            continue
        plat = point['latitude'].values
        plon = point['longitude'].values


        xpos = np.where(probs['longitude'].values == plon)
        xpos = int(xpos[0])
        ypos = np.where(probs['latitude'].values == plat)
        ypos = int(ypos[0])

        try:
            vdic, cntera = cut_kernel(xpos, ypos, probs, dist)
        except TypeError:
            print('TypeError')
            continue

        cnte_sum = np.nansum(np.stack([cnte_sum, cntera]), axis=0)

        for ks in vdic.keys():
            if ks in edic:
                edic[ks] = np.nansum(np.stack([edic[ks], vdic[ks]]), axis=0)
            else:
                edic[ks] = vdic[ks]
        counter += 1
        print('Saved core ', counter)

    if not edic:
        return None

    outlist = [cnte_sum]
    outnames = ['cnte']
    for ek in edic.keys():
        outnames.append(ek)
        outlist.append(edic[ek])

    print('Returning')

    return outlist, outnames


def plot_big_SN(h, eh):
    hour = h
    dic = {}

    name = "ERA5_pressureLevels_SN_2hOverlap_AMSRE_SMALL_"

    def coll(dic2, h, eh, year):
        print(h)
        core = pkl.load(open(
            cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/core_txt/" + name + str(eh) + "UTCERA" + str(h).zfill(
                2) + '_' + str(year) + ".p", "rb"))



        for id, k in enumerate(core.keys()):
            try:
                dic2[k] = dic2[k] + core[k]
            except KeyError:
                dic2[k] = core[k]

    for y in range(2006, 2011):
        coll(dic, h, eh, y)


    levels = [300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 825, 850, 875, 900, 925, 950, 975]
    #ipdb.set_trace()
    f = plt.figure(figsize=(15, 8))
    ax = f.add_subplot(231)
    plt.contourf(np.arange(-22, 23) * 27, levels, (dic['u']) / dic['cnte'], cmap='RdBu_r',
                 levels=np.arange(-12, 2.1, 1), extend='both')  # -(rkernel2_sum / rcnt_sum)
    plt.colorbar(label='m s-1')
    plt.gca().invert_yaxis()
    # plt.gca().invert_xaxis()
    plt.title('u wind', fontsize=9)
    plt.xlabel('West-East extent')
    plt.ylabel('Pressure level (hPa)')

    plt.ylim(975,400)

    ax1 = f.add_subplot(232)
    plt.contourf(np.arange(-22, 23) * 27, levels, ((dic['v']-dic['vclim']) / dic['cnte']), extend='both', cmap='RdBu_r',
                 levels=np.arange(-1.6, 1.7, 0.2))  # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='m s-1')
    plt.gca().invert_yaxis()
    # plt.gca().invert_xaxis()
    plt.title('v wind', fontsize=9)
    plt.xlabel('West-East extent')
    plt.ylabel('Pressure level (hPa)')
    plt.ylim(975, 400)

    ax1 = f.add_subplot(233)
    plt.contourf(np.arange(-22, 23) * 27, levels, ((dic['q']-dic['qclim']) * 1000 / dic['cnte']), extend='both', cmap='RdBu',
                 levels=np.arange(-1, 1.1, 0.1))  # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='g kg-1')
    plt.gca().invert_yaxis()
    # plt.gca().invert_xaxis()
    plt.title('specific humidity', fontsize=9)
    plt.xlabel('West-East extent')
    plt.ylabel('Pressure level (hPa)')
    plt.ylim(975, 400)

    ax1 = f.add_subplot(234)
    plt.contourf(np.arange(-22, 23) * 27, levels, (dic['div']/ dic['cnte']), extend='both', cmap='RdBu_r', levels=np.arange(-0.005, 0.0051, 0.001))  # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='K')
    plt.gca().invert_yaxis()
    # plt.gca().invert_xaxis()
    plt.title('Divergence', fontsize=9)
    plt.xlabel('West-East extent')
    plt.ylabel('Pressure level (hPa)')
    plt.ylim(975, 400)

    ax1 = f.add_subplot(235)
    theta = dic['theta']-dic['theta_clim']
    plt.contourf(np.arange(-22, 23) * 27, levels, (theta / dic['cnte']), extend='both', cmap='RdBu_r',
                 vmin=-0.4,
                 vmax=0.4)  # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='K')
    plt.gca().invert_yaxis()
    # plt.gca().invert_xaxis()
    plt.title(r'$\theta$', fontsize=9)
    plt.xlabel('West-East extent')
    plt.ylabel('Pressure level (hPa)')
    plt.ylim(975, 400)

    ax1 = f.add_subplot(236)
    theta = dic['z']-dic['zclim']
    plt.contourf(np.arange(-22, 23) * 27, levels, (theta / dic['cnte']), extend='both', cmap='RdBu_r', levels=np.linspace(-1,1,10))  # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='K')
    plt.gca().invert_yaxis()
    # plt.gca().invert_xaxis()
    plt.title(r'$\theta$', fontsize=9)
    plt.xlabel('West-East extent')
    plt.ylabel('Pressure level (hPa)')
    plt.ylim(975, 700)



    plt.tight_layout()
    # plt.savefig('/users/global/cornkle/figs/LSTA-bullshit/corrected_LSTA/system_scale/doug/' + str(hour).zfill(
    #     2) + '_WE_pl_big.png')  # str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png')
    #plt.close()


    #np.arange(-20,21)*30, levels,


def plot_big_WE(h, eh):
    hour = h
    dic = {}

    name = "ERA5_pressureLevels_WE_2hOverlap_AMSRE_SMALL_"

    def coll(dic2, h, eh, year):
        print(h)
        core = pkl.load(open(
            cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/core_txt/" + name + str(eh) + "UTCERA" + str(h).zfill(
                2) + '_' + str(year) + ".p", "rb"))



        for id, k in enumerate(core.keys()):
            try:
                dic2[k] = dic2[k] + core[k]
            except KeyError:
                dic2[k] = core[k]

    for y in range(2006, 2011):
        coll(dic, h, eh, y)


    levels = [300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 825, 850, 875, 900, 925, 950, 975]
    #ipdb.set_trace()
    f = plt.figure(figsize=(15, 8))
    ax = f.add_subplot(231)
    plt.contourf(np.arange(-22, 23) * 27, levels, (dic['u']) / dic['cnte'], cmap='RdBu_r',
                 levels=np.arange(-12, 2.1, 1), extend='both')  # -(rkernel2_sum / rcnt_sum)
    plt.colorbar(label='m s-1')
    plt.gca().invert_yaxis()
    # plt.gca().invert_xaxis()
    plt.title('u wind', fontsize=9)
    plt.xlabel('West-East extent')
    plt.ylabel('Pressure level (hPa)')

    ax1 = f.add_subplot(232)
    plt.contourf(np.arange(-22, 23) * 27, levels, ((dic['v']-dic['vclim']) / dic['cnte']), extend='both', cmap='RdBu_r',
                 levels=np.arange(-1.6, 1.7, 0.2))  # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='m s-1')
    plt.gca().invert_yaxis()
    # plt.gca().invert_xaxis()
    plt.title('v wind', fontsize=9)
    plt.xlabel('West-East extent')
    plt.ylabel('Pressure level (hPa)')

    ax1 = f.add_subplot(233)
    plt.contourf(np.arange(-22, 23) * 27, levels, ((dic['q']-dic['qclim']) * 1000 / dic['cnte']), extend='both', cmap='RdBu',
                 levels=np.arange(-1, 1.1, 0.1))  # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='g kg-1')
    plt.gca().invert_yaxis()
    # plt.gca().invert_xaxis()
    plt.title('specific humidity', fontsize=9)
    plt.xlabel('West-East extent')
    plt.ylabel('Pressure level (hPa)')

    ax1 = f.add_subplot(234)
    plt.contourf(np.arange(-22, 23) * 27, levels, (dic['div']/ dic['cnte']), extend='both', cmap='RdBu_r', levels=np.arange(-0.005, 0.0051, 0.001))  # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='K')
    plt.gca().invert_yaxis()
    # plt.gca().invert_xaxis()
    plt.title('Divergence', fontsize=9)
    plt.xlabel('West-East extent')
    plt.ylabel('Pressure level (hPa)')

    ax1 = f.add_subplot(235)
    theta = dic['theta']-dic['theta_clim']
    plt.contourf(np.arange(-22, 23) * 27, levels, (theta / dic['cnte']), extend='both', cmap='RdBu_r',
                 vmin=-0.4,
                 vmax=0.4)  # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='K')
    plt.gca().invert_yaxis()
    # plt.gca().invert_xaxis()
    plt.title(r'$\theta$', fontsize=9)
    plt.xlabel('West-East extent')
    plt.ylabel('Pressure level (hPa)')

    ax1 = f.add_subplot(236)
    theta = dic['z']-dic['zclim']
    plt.contourf(np.arange(-22, 23) * 27, levels, (theta / dic['cnte']), extend='both', cmap='RdBu_r', levels=np.linspace(-5,5,10))  # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='K')
    plt.gca().invert_yaxis()
    # plt.gca().invert_xaxis()
    plt.title(r'$\theta$', fontsize=9)
    plt.xlabel('West-East extent')
    plt.ylabel('Pressure level (hPa)')
    plt.ylim(975, 700)


    plt.tight_layout()
    # plt.savefig('/users/global/cornkle/figs/LSTA-bullshit/corrected_LSTA/system_scale/doug/' + str(hour).zfill(
    #     2) + '_WE_pl_big.png')  # str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png')
    #plt.close()