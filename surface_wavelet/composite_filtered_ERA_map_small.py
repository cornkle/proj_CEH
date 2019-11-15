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


def composite(h, eh):

    msgopen = pd.read_pickle(
        '/home/ck/DIR/cornkle/figs/LSTA/corrected_LSTA/new/cores_-60_filteredPoints_gt25000km2_table_new.p')
    hour = h
    msg = pd.DataFrame.from_dict(msgopen)
    msg = msg[(msg.hour==h)]# &  &
    msg['eh'] = eh
    print('Start core number ', len(msg))

    # calculate the chunk size as an integer
    #'chunk_size = int(msg.shape[0] / pnumber)
    msg.sort_values(by='date')
    chunk, chunk_ind, chunk_count = np.unique(msg.date, return_index=True, return_counts=True)


    chunks = [msg.ix[msg.index[ci:ci + cc]] for ci, cc in zip(chunk_ind, chunk_count)] # daily chunks

    # res = []
    # for m in chunks[0:30]:
    #     out = file_loop(m)
    #     res.append(out)
    #
    # ipdb.set_trace()
    # return
    dic = u_parallelise.era_run_arrays(4, file_loop, chunks)

    pkl.dump(dic, open(cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/ERA5_mapSmallfunc_"+str(eh) + "UTCERA"+str(hour).zfill(2)+"_cores.p", "wb"))
    del dic
    print('Dumped file')



def cut_kernel(xpos, ypos, arr, dist, probs=False):

    kernel = u_arrays.cut_kernel(arr,xpos, ypos,dist)

    vdic = {}

    for d in probs.data_vars:

        var = u_arrays.cut_kernel(probs[d].values,xpos, ypos,dist)
        vdic[d] = var

    cnt2 = np.zeros_like(kernel)
    cnt2[np.isfinite(vdic[list(vdic.keys())[0]])] = 1


    cnt = np.zeros_like(kernel)
    cnt[np.isfinite(kernel)] = 1


    if kernel.shape != (dist*2+1, dist*2+1):
        print('Kernels shape wrong!')
        ipdb.set_trace()

    #kernel = kernel - np.nanmean(kernel)


    return kernel, cnt, vdic, cnt2

def get_previous_hours(date, ehour, refhour):


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
    #css = css.sel(longitude=slice(-13, 13), latitude=slice(8, 21))
    pl_clim = pl_clim.sel(longitude=slice(-13, 13), latitude=slice(8, 21))
    #srfc_clim = srfc_clim.sel(longitude=slice(-13, 13), latitude=slice(8, 21))


    cmm = cmm.sel(time=t1)
    pl_clim = pl_clim.squeeze()

    css = css.sel(time=t1)
    #srfc_clim = srfc_clim.squeeze()

    t = cmm['t'].sel(level=925).squeeze() - pl_clim['t'].sel(level=925).squeeze() #* 1000

    shear =  (cmm['u'].sel(level=650).squeeze() - cmm['u'].sel(level=925).squeeze() ) #- (pl_clim['u'].sel(level=600).squeeze() - pl_clim['u'].sel(level=925).squeeze() ) #

    vwind_srfc = cmm['v'].sel(level=925).squeeze() - pl_clim['v'].sel(level=925).squeeze()
    uwind_srfc = cmm['u'].sel(level=925).squeeze() - pl_clim['u'].sel(level=925).squeeze()


    uwind_up = cmm['u'].sel(level=650).squeeze() - pl_clim['u'].sel(level=650).squeeze()

    #wwind_srfc = cmm['w'].sel(level=400).squeeze() - pl_clim['w'].sel(level=400).squeeze()
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

    cm['u650'] = uwind_up

    cm['u925raw'] = cmm['u'].sel(level=925).squeeze()
    cm['v925raw'] = cmm['v'].sel(level=925).squeeze()

    cm['div'] = div *1000
    cm['q'] = q
    cm['t'] = t
    cm['theta_e'] = theta_e

    #del srfc_clim
    del pl_clim
    del css
    del cmm
    return cm


def file_loop(df):

    date = df['date'].iloc[0]
    eh = df['eh'].iloc[0]
    hour = df['hour'].iloc[0]
    print('Doing day: ', date)

    fdate = str(date.year) + str(date.month).zfill(2) + str(date.day).zfill(2)  # date of LSTA

    try:
        lsta = xr.open_dataset(cnst.LSTA_NEW + 'lsta_daily_' + fdate + '.nc')
    except OSError:
        print('No lsta, return')
        return
    print('Doing ' + 'lsta_daily_' + fdate + '.nc')

    dist = 200
    lsta_da = lsta['LSTA'].squeeze()

    edic = {}

    probs = get_previous_hours(date, eh, hour)
    print('Era5 collect')

    try:
        probs_on_lsta = lsta.salem.transform(probs)
    except RuntimeError:
        print('Era5 on LSTA interpolation problem')
        return
    del probs
    del lsta
    counter = 0

    kernel2_sum = np.zeros((dist * 2 + 1, dist * 2 + 1))
    cnt_sum = np.zeros((dist * 2 + 1, dist * 2 + 1))
    cntm_sum = np.zeros((dist * 2 + 1, dist * 2 + 1))

    for index,fi in df.iterrows():


        xpos = fi.x_5k
        ypos = fi.y_5k

        try:
            kernel2, cnt, vdic, ecnt = cut_kernel(xpos, ypos, lsta_da, dist, probs=probs_on_lsta)
        except TypeError:
            print('No kernel for core, continue')
            continue

        kernel2_sum = np.nansum(np.stack([kernel2_sum, kernel2]), axis=0)
        cnt_sum = np.nansum(np.stack([cnt_sum, cnt]), axis=0)
        cntm_sum = np.nansum(np.stack([cntm_sum, ecnt]), axis=0)

        for ks in vdic.keys():
            if ks in edic:
                edic[ks] = np.nansum(np.stack([edic[ks], vdic[ks]]), axis=0)
            else:
                edic[ks] = vdic[ks]
        counter += 1
        print('Saved core ', counter)

    if np.sum(cnt_sum)==0:
        print('Lsta valid count is zero, continue')
        return

    outlist = [kernel2_sum, cnt_sum, cntm_sum]
    outnames = ['lsta',  'cnt', 'ecnt']
    for ek in edic.keys():
        outnames.append(ek)
        outlist.append(edic[ek])

    print('Returning')
    return outlist, outnames



def plot_doug(h, eh):

    # dic = {}
    # dic2 = {}
    #
    # def coll(dic, h, eh, year):
    #     print(h)
    #     core = pkl.load(open(
    #         cnst.network_data + "figs/LSTA/corrected_LSTA/new/wavelet_drywet_power/composite_backtrack"+str(eh) + "UTCERA"+str(h).zfill(2)+'_'+str(year)+"_small_cores.p", "rb"))
    #     for id, k in enumerate(core.keys()):
    #         try:
    #             dic[k] = dic[k] + core[k]
    #         except KeyError:
    #             dic[k] = core[k]

        # cm = pkl.load(open(
        #     cnst.network_data + "figs/LSTA/corrected_LSTA/new/composite_backtrack_CMORPH_" + str(h).zfill(
        #         2) + "_cores.p", "rb"))

        # for id, k in enumerate(cm.keys()):
        #     try:
        #         dic2[k] = dic2[k] + cm[k]
        #     except KeyError:
        #         dic2[k] = cm[k]

    dic = pkl.load(open(
        cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/ERA5_mapSmallfunc_" + str(eh) + "UTCERA" + str(
            h).zfill(2) + "_cores.p", "rb"))


    # for y in range(2006,2011):
    #     coll(dic, h, eh, y)

    extent = (dic['lsta'].shape[1]-1)/2
    xlen = dic['lsta'].shape[1]
    ylen = dic['lsta'].shape[0]

    xv, yv = np.meshgrid(np.arange(ylen), np.arange(xlen))
    st=30
    xquiv = xv[4::st, 4::st]
    yquiv = yv[4::st, 4::st]

    u = (dic['u925']/ dic['ecnt'])[4::st, 4::st]
    v = (dic['v925']/ dic['ecnt'])[4::st, 4::st]


    f = plt.figure(figsize=(10,8))
    ax = f.add_subplot(221)

    plt.contourf((dic['lsta'] / dic['cnt']), cmap='RdBu_r', levels=[-0.8, -0.6,-0.4,-0.2,-0.1, 0.1,0.2,0.4,0.6,  0.8], extend='both') #-(rkernel2_sum / rcnt_sum)
    #plt.plot(extent, extent, 'bo')
    plt.colorbar(label='K')
    #pdb.set_trace()

    #contours = plt.contour((dic2['prob']/ dic2['cntp'])*100, extend='both', levels=np.arange(10,70,10), cmap='jet') # #, levels=np.arange(1,5, 0.5)
    #plt.clabel(contours, inline=True, fontsize=11, fmt='%1.0f')
    plt.plot(extent, extent, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6 , dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.title('1800UTC | '+str(np.max(dic['ecnt']))+' cores, LSTA', fontsize=9)


    ax1 = f.add_subplot(222)
    plt.contourf(((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r',levels=[ -0.8, -0.6,-0.4,-0.2,-0.1, 0.1,0.2,0.4,0.6,  0.8]) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='K')
    contours = plt.contour((dic['t'] / dic['ecnt']), extend='both',levels=[ -0.8,-0.7, -0.6,-0.5,-0.4,-0.2,-0.1, 0.1,0.2,0.4,0.5,0.6, 0.7, 0.8], cmap='PuOr_r') #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.1f')
    #qu = ax1.quiver(xquiv, yquiv, u, v, scale=50)
    plt.plot(extent, extent, 'bo')
    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('Shading: MSG LSTA, Contours: ERA5 925hPa T', fontsize=9)

    ax1 = f.add_subplot(223)
    plt.contourf(((dic['q'])*1000/ dic['ecnt']), extend='both',  cmap='RdBu',levels=np.arange(-1,1.1,0.1)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='g kg-1')
    contours = plt.contour((dic['shear'] / dic['ecnt']), extend='both',levels=np.arange(-17,-12,0.25), cmap='viridis') #np.arange(-15,-10,0.5)
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
    plt.contourf(((dic['div'])/ dic['ecnt'])*100, extend='both',  cmap='RdBu', levels=np.linspace(-0.5,0.5,10)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
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
    #plt.savefig('/users/global/cornkle/figs/LSTA-bullshit/corrected_LSTA/system_scale/doug/large_scale/'+str(hour).zfill(2)+'_JUN.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png')
    #plt.close()

def plot_doug_small(h, eh):

    dic = {}

    # def coll(dic, h, eh, year):
    #     print(h)
    #     core = pkl.load(open(
    #         cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/ERA5_mapSmallfunc_"+str(eh) + "UTCERA"+str(hour).zfill(2)+"_cores.p", "rb"))
    #     for id, k in enumerate(core.keys()):
    #         try:
    #             dic[k] = dic[k] + core[k]
    #         except KeyError:
    #             dic[k] = core[k]
    #
    #
    # for y in range(2006, 2011):
    #     coll(dic, h, eh, y)

    dic = pkl.load(open(
        cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/ERA5_mapSmallfunc_" + str(eh) + "UTCERA" + str(
            h).zfill(2) + "_cores.p", "rb"))

    extent = (dic['lsta'].shape[1] - 1) / 2
    xlen = dic['lsta'].shape[1]
    ylen = dic['lsta'].shape[0]

    xv, yv = np.meshgrid(np.arange(ylen), np.arange(xlen))
    st = 30
    xquiv = xv[4::st, 4::st]
    yquiv = yv[4::st, 4::st]

    u = (dic['u925'] / dic['ecnt'])[4::st, 4::st]
    v = (dic['v925'] / dic['ecnt'])[4::st, 4::st]

    f = plt.figure(figsize=(10, 4))
    ax = f.add_subplot(121)

    plt.contourf((dic['lsta'] / dic['cnt']), cmap='RdBu_r',
                 levels=[-0.8, -0.6, -0.4, -0.2, -0.1, 0.1, 0.2, 0.4, 0.6, 0.8],
                 extend='both')  # -(rkernel2_sum / rcnt_sum)
    # plt.plot(extent, extent, 'bo')
    plt.colorbar(label='K')
    # pdb.set_trace()

    contours = plt.contour(((dic['div'])/ dic['ecnt'])*100, extend='both',  cmap='RdBu', levels=np.linspace(-0.8,0.8,8)) # #, levels=np.arange(1,5, 0.5)
    qu = ax.quiver(xquiv, yquiv, u, v, scale=15)
    qk = plt.quiverkey(qu, 0.9, 0.02,1, '1 m s$^{-1}$',
                       labelpos='E', coordinates='figure')


    plt.clabel(contours, inline=True, fontsize=11, fmt='%1.0f')
    plt.plot(extent, extent, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, extent * 2, 9) - 100) * 6, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent * 2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.title('23-01UTC | ' + str(np.max(dic['ecnt'])) + ' cores, LSTA & 06-06UTC antecedent rain', fontsize=9)


    ax1 = f.add_subplot(122)
    plt.contourf(((dic['q'])*1000/ dic['ecnt']), extend='both',  cmap='RdBu',levels=np.arange(-1,1.1,0.1)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='g kg-1')
    contours = plt.contour((dic['shear'] / dic['ecnt']), extend='both',levels=np.arange(-18,-12,0.5), cmap='viridis') #np.arange(-15,-10,0.5)
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
    # plt.savefig('/users/global/cornkle/figs/LSTA-bullshit/corrected_LSTA/system_scale/doug/large_scale/'+str(hour).zfill(2)+'_JUN.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png')
    # plt.close()


def plot_doug_small_allhours():

    dic = {}

    def coll(dic, h, eh, year):
        print(h)
        core = pkl.load(open(
            cnst.network_data + "figs/LSTA/corrected_LSTA/new/august/composite_backtrack" + str(
                eh) + "UTCERA" + str(h).zfill(2) + '_' + str(year) + "_small_cores.p", "rb"))
        for id, k in enumerate(core.keys()):
            try:
                dic[k] = dic[k] + core[k]
            except KeyError:
                dic[k] = core[k]


    for y in range(2006, 2011):
        for h,eh in zip([17,18,19,20], [-5,-6,-7,-8]):
            coll(dic, h, eh, y)

    extent = (dic['lsta'].shape[1] - 1) / 2
    xlen = dic['lsta'].shape[1]
    ylen = dic['lsta'].shape[0]

    xv, yv = np.meshgrid(np.arange(ylen), np.arange(xlen))
    st = 30
    xquiv = xv[4::st, 4::st]
    yquiv = yv[4::st, 4::st]

    u = (dic['u925'] / dic['cntp'])[4::st, 4::st]
    v = (dic['v925'] / dic['cntp'])[4::st, 4::st]

    f = plt.figure(figsize=(10, 4))
    ax = f.add_subplot(121)

    plt.contourf((dic['lsta'] / dic['cnt']), cmap='RdBu_r',
                 levels=[-1, -0.8, -0.6, -0.4, -0.2, -0.1, 0.1, 0.2, 0.4, 0.6, 0.8,1],
                 extend='both', alpha=0.9)  # -(rkernel2_sum / rcnt_sum)
    # plt.plot(extent, extent, 'bo')
    plt.colorbar(label='K')
    # pdb.set_trace()

    contours = plt.contour(((dic['div'])/ dic['cntp'])*100, extend='both',  cmap='RdBu', levels=np.linspace(-0.49,0.49,6)) # #, levels=np.arange(1,5, 0.5)
    qu = ax.quiver(xquiv, yquiv, u, v, scale=13, width=0.006)
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
    plt.contourf(((dic['q'])*1000/ dic['cntp']), extend='both',  cmap='RdBu',levels=np.arange(-0.5,0.6,0.05)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='g kg-1')
    contours = plt.contour((dic['shear'] / dic['cntp']), extend='both',levels=np.arange(-16,-12,0.5), cmap='viridis') #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    plt.plot(extent, extent, 'bo')
    #qu = ax1.quiver(xquiv, yquiv, u, v, scale=50)
    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) -100) * 6, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('Shading: 925hPa q anomaly, Contours: 650hPa-925hPa wind shear ', fontsize=9)

    plt.tight_layout()
    plt.show()
    plt.savefig(cnst.network_data + "figs/LSTA/corrected_LSTA/new/august/jun-sep_17-20UTC_mean.png")#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png')
    # plt.close()


def plot_all():

    hours = [16,17,18,19,20,21,22,23,0,1,2,3,4,5,6,7]

    for h in hours:
        plot_doug(h)
