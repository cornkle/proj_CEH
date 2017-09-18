# -*- coding: utf-8 -*-


import glob
import numpy as np
from wavelet import util
import xarray as xr
import os
import ipdb
import matplotlib.pyplot as plt
from scipy import ndimage
from utils import u_arrays as ua
import multiprocessing
from scipy.ndimage.measurements import label
import pdb

def run():
    #  (1174, 378)
    #gridsat_folder = '/users/global/cornkle/data/OBS/gridsat/gridsat_netcdf/z18/yearly/'
    gridsat_folder = '/users/global/cornkle/data/OBS/gridsat/gridsat_netcdf/yearly_files/'
    pool = multiprocessing.Pool(processes=7)

    files = glob.glob(gridsat_folder+'gridsat_*.nc')

    res = pool.map(file_loop_circle, files)

    pool.close()

def file_loop_circle(files):

    out = '/users/global/cornkle/data/OBS/gridsat/gridsat_netcdf/circles/single/Mar-Oct-69_18UTC/'

    da = xr.open_dataset(files)
    da = da['t']
    lat = da['lat']
    lon = da['lon']
    year = np.unique(da['time.year'].values)[0]

    da_list = []

    da = da[(da['time.month']>=3) & (da['time.month']<=10)]

    for time in da['time']: # [730:1030]

        print(time)

        outt = np.array(da.sel(time=time).values, dtype=float)
        figure = outt.copy() * 0.
        outt[outt>-40] = np.nan
        image = da.sel(time=time).values

        ## remove small
        image[image >= -40] = 0
        labels, numL = label(image)
        mcs = []
        for gi in np.unique(labels):
            if gi == 0:  # index 0 is always background, ignore!
                continue
            pos = np.where(labels == gi)
            #### only MCSs > 15000km2 with pixel = 64km2
            if np.sum(labels == gi) < 234:
                image[pos] = 0
                continue
            mcs.append(pos)

        # ltgrad = []
        # lpower = []
        for p in mcs:
            yb, xb = p
            smcs = image[np.min(yb):np.max(yb), np.min(xb):np.max(xb)]
            bbla = outt[np.min(yb):np.max(yb), np.min(xb):np.max(xb)]
            fig_box = outt[np.min(yb):np.max(yb), np.min(xb):np.max(xb)]*0.

            ##mark MCS edges with gradient
            smcs[smcs == 0] = 150
            smcs[np.isnan(smcs)] = 150
            grad = np.gradient(smcs)
            smcs[smcs == 150] = -55

            # f = plt.figure()
            # plt.imshow(smcs)

            ###smooth MCS edges in a 3x3 pixel kernel
            nok = np.where(abs(grad[0]) > 80)
            d = 2
            i = nok[0]
            j = nok[1]
            for ii, jj in zip(i, j):
                kern = smcs[ii - d:ii + d + 1, jj - d:jj + d + 1]
                smcs[ii - d:ii + d + 1, jj - d:jj + d + 1] = ndimage.gaussian_filter(kern, 2, mode='nearest')

            ##apply wavelet transform
            wav = util.waveletT8(smcs, 8)
            arr = np.array(wav['scales'], dtype=float)

            arrint = np.array(wav['scales'], dtype=int)

            scale_ind = range(arr.size)
            wll = wav['t']

            # yyy = []
            # xxx = []
            # scal = []
            for nb in scale_ind[::-1]:

                orig = float(arr[nb])

                if orig > 35:  # > 30:
                    continue

                scale = int(np.round(orig))

                #print(np.round(orig))

                wl = wll[nb, :, :]
                # maxout = maxoutt[nb, :, :]

                maxout = (
                    wl == ndimage.maximum_filter(wl, (4, 4), mode='constant',
                                                 cval=np.amax(wl) + 1))  # (np.round(orig / 5))

                try:
                    yy, xx = np.where((maxout == 1) & (bbla <= -69) & ((wl >= np.percentile(wl[wl >= 1], 90)) & (
                    wl > orig ** .5)))  # )& (wl > orig**.5) (wl >= np.percentile(wl[wl >= 0.1], 90)) )#(wl > orig**.5))#  & (wlperc > orig**.5))# & (wlperc > np.percentile(wlperc[wlperc>=0.1], 80)))# & (wlperc > np.percentile(wlperc[wlperc>=0.1], 80) ))  # & (wl100 > 5)
                except IndexError:
                    continue

                for y, x in zip(yy, xx):
                    ss = orig
                    iscale = (np.ceil(ss / 2. / 8.)).astype(int)

                    ycirc, xcirc = ua.draw_cut_circle(x, y, iscale, bbla)

                    fig_box[ycirc, xcirc] = 1
                    tgrad = np.nanmax(bbla[ycirc, xcirc]) - np.nanmin(bbla[ycirc, xcirc])
                    wmean = np.nanmax(wl[ycirc, xcirc])
                    # xxx.append(x)
                    # yyy.append(y)
                    # ltgrad.append(tgrad)
                    # lpower.append(wmean)
                    # scal.append(orig)

            figure[np.min(yb):np.max(yb), np.min(xb):np.max(xb)] = fig_box

        da_new = xr.DataArray(figure, coords={'time': time, 'lat': lat, 'lon': lon},
                              dims=['lat', 'lon'])
        da_list.append(da_new)
    cube = xr.concat(da_list, dim='time')

    cube.name = 'power'

    cube.to_netcdf(out + 'gs_scale_circle_Mar-Oct-69_18UTC_'+str(year)+'.nc')
    #
    # f=plt.figure()
    # plt.scatter(ltgrad, lpower)

def concat_files():
    y1 = 1982
    y2 = 2017
    years = np.arange(y1 + 1, y2)  # 2017)

    msg_folder = '/users/global/cornkle/data/OBS/gridsat/gridsat_netcdf/circles/single/MAM-69/'
    out = '/users/global/cornkle/data/OBS/gridsat/gridsat_netcdf/circles/'
    fname = 'gs_scale_circle_MAM-69_1983-2016.nc'

    arlist = []

    for y in years:
        y = str(y)
        da = xr.open_dataarray(msg_folder + 'gs_scale_circle_MAM-69_' + str(y) + '.nc')
        da[np.isnan(da)]=0
        print('Doing ' + y)

        arlist.append(da)

    all = xr.concat(arlist, 'time', compat='identical')

    #enc = {'power': {'complevel': 5, 'zlib': True}}
    all.to_netcdf(out + fname)#, encoding=enc)


def file_loop_blob_slice(files):

    spath = '/users/global/cornkle/data/OBS/gridsat/gridsat_netcdf/'

    da = xr.open_dataset(files)
    lat = da['lat']
    lon = da['lon']
    year = np.unique(da['time.year'].values)[0]

    da_list = []

    #tcube = da.isel(time=list(np.arange(730,733)))

    #tcube.to_netcdf(spath+'gs_temp_testslice4.nc')
    lt = []
    ltgrad = []
    lpower = []

    for time in da['time'][730:740]: # [730:1030]

        bla = np.array(da['t'].sel(time=time).values, dtype=float)
        bla[bla>-40] = np.nan
        image = da['t'].sel(time=time).values

        ## remove small
        image[image >= -40] = 0
        labels, numL = label(image)
        mcs = []
        for gi in np.unique(labels):
            if gi == 0:  # index 0 is always background, ignore!
                continue
            pos = np.where(labels == gi)
            #### only MCSs > 15000km2 with pixel = 64km2
            if np.sum(labels == gi) < 234:
                image[pos] = 0
                continue

            #print(np.min(bla[pos]), np.max(bla[pos]))
            mcs.append(pos)

        dummy = util.waveletT(image, 8)
        dummy = dummy['t'].copy()*0

        for p in mcs:
            y, x = p
            smcs = image[np.min(y):np.max(y), np.min(x):np.max(x)]
            bbla = bla[np.min(y):np.max(y), np.min(x):np.max(x)]

            ##mark MCS edges with gradient
            smcs[smcs == 0] = 150
            smcs[np.isnan(smcs)] = 150
            grad = np.gradient(smcs)
            smcs[smcs == 150] = -65

            # f = plt.figure()
            # plt.imshow(smcs)

            ###smooth MCS edges in a 3x3 pixel kernel
            nok = np.where(abs(grad[0]) > 80)
            d = 2
            i = nok[0]
            j = nok[1]
            for ii, jj in zip(i, j):
                kern = smcs[ii - d:ii + d + 1, jj - d:jj + d + 1]
                smcs[ii - d:ii + d + 1, jj - d:jj + d + 1] = ndimage.gaussian_filter(kern, 3, mode='nearest')

            ##apply wavelet transform
            wav = util.waveletT8(smcs, 8)
            arr = np.array(wav['scales'], dtype=float)

            scale_ind = range(arr.size)
            wl = wav['t']
            dummy2 = wav['t'].copy()*0
            for nb in scale_ind[::-1]:

                orig = float(arr[nb])

                wl = wav['t'][nb,:,:]
                try:
                    wl[(wl < orig**.5) | (wl < 5) | (wl < np.percentile(wl[wl >= 0.5], 50)) ] = 0
                except IndexError:
                    wl[(wl < orig ** .5) | (wl < 5)] = 0

                wl[np.isnan(bbla)] = 0

                wlabels, num = label(wl)

                for gi in np.unique(wlabels):

                    if gi == 0:  # index 0 is always background, ignore!
                        continue

                    pos = np.where(wlabels == gi)
                    tmean = np.median(bbla[pos])

                    tgrad = np.percentile(bbla[pos], 95) - np.percentile(bbla[pos], 5)
                    if not np.isfinite(tgrad):
                        print('Not finite')
                        return

                    wmean = np.median(wl[pos])

                    if (tmean > -55) or (wmean < 2):
                        wl[pos] = 0
                        continue

                    lt.append(tmean)
                    ltgrad.append(tgrad)
                    lpower.append(wmean)

                dummy2[nb, :, :] = wl

            dummy[: ,np.min(y):np.max(y),np.min(x):np.max(x)] = dummy2

    #     da_new = xr.DataArray(dummy, coords={'scales': wav['scales'],'time': time, 'lat': lat, 'lon': lon}, dims=['scales', 'lat', 'lon'])
    #     da_list.append(da_new)
    # cube = xr.concat(da_list, dim='time')
    #
    # cube.name = 'power'
    # enc = {'power': {'complevel': 5, 'zlib': True}}
    #cube.to_netcdf(spath + 'gs_power_testslice4.nc', encoding=enc)
    f = plt.figure()
    plt.scatter(ltgrad, lpower, c=lt, cmap='viridis')
    plt.xlabel('SCF gradient')
    plt.ylabel('SCF median power')
    plt.colorbar(label='SCF T')


def file_loop_circle_slice(files):

    spath = '/users/global/cornkle/data/OBS/gridsat/gridsat_netcdf/'

    da = xr.open_dataset(files)
    lat = da['lat']
    lon = da['lon']
    year = da['time.year'][0]

    da_list = []

    tcube = da.isel(time=list(np.arange(730,733)))

    tcube.to_netcdf(spath+'gs_temp_circle_slice.nc')

    for time in da['time'][730:740]: # [730:1030]
        print(time)
        outt = np.array(da['t'].sel(time=time).values, dtype=float)
        outt[outt>-40] = np.nan
        image = da['t'].sel(time=time).values

        ## remove small
        image[image >= -40] = 0
        labels, numL = label(image)
        mcs = []
        for gi in np.unique(labels):
            if gi == 0:  # index 0 is always background, ignore!
                continue
            pos = np.where(labels == gi)
            #### only MCSs > 15000km2 with pixel = 64km2
            if np.sum(labels == gi) < 234:
                image[pos] = 0
                continue
            mcs.append(pos)

        dummy = util.waveletT(image, 8)
        figure = outt.copy()*0.
        ltgrad = []
        lpower = []
        ltt = []
        for p in mcs:
            yb, xb = p
            smcs = image[np.min(yb):np.max(yb), np.min(xb):np.max(xb)]
            bbla = outt[np.min(yb):np.max(yb), np.min(xb):np.max(xb)]
            fig_box = outt[np.min(yb):np.max(yb), np.min(xb):np.max(xb)]*0.

            ##mark MCS edges with gradient
            smcs[smcs == 0] = 150
            smcs[np.isnan(smcs)] = 150
            grad = np.gradient(smcs)
            smcs[smcs == 150] = -65

            # f = plt.figure()
            # plt.imshow(smcs)

            ###smooth MCS edges in a 3x3 pixel kernel
            nok = np.where(abs(grad[0]) > 80)
            d = 2
            i = nok[0]
            j = nok[1]
            for ii, jj in zip(i, j):
                kern = smcs[ii - d:ii + d + 1, jj - d:jj + d + 1]
                smcs[ii - d:ii + d + 1, jj - d:jj + d + 1] = ndimage.gaussian_filter(kern, 3, mode='nearest')

            ##apply wavelet transform
            wav = util.waveletT8(smcs, 8)
            arr = np.array(wav['scales'], dtype=float)

            arrint = np.array(wav['scales'], dtype=int)

            scale_ind = range(arr.size)
            wll = wav['t']

            yyy = []
            xxx = []
            scal = []
            for nb in scale_ind[::-1]:

                orig = float(arr[nb])

                if orig > 35:  # > 30:
                    continue

                scale = int(np.round(orig))

                print(np.round(orig))

                wl = wll[nb, :, :]
                # maxout = maxoutt[nb, :, :]

                maxout = (
                    wl == ndimage.maximum_filter(wl, (5, 5), mode='constant',
                                                 cval=np.amax(wl) + 1))  # (np.round(orig / 5))

                try:
                    yy, xx = np.where((maxout == 1) & (bbla <= -60) & ((wl >= np.percentile(wl[wl >= 1], 90)) & (
                    wl > orig ** .5)))  # )& (wl > orig**.5) (wl >= np.percentile(wl[wl >= 0.1], 90)) )#(wl > orig**.5))#  & (wlperc > orig**.5))# & (wlperc > np.percentile(wlperc[wlperc>=0.1], 80)))# & (wlperc > np.percentile(wlperc[wlperc>=0.1], 80) ))  # & (wl100 > 5)
                except IndexError:
                    continue

                for y, x in zip(yy, xx):
                    ss = orig
                    iscale = (np.ceil(ss / 2. / 8.)).astype(int)

                    ycirc, xcirc = ua.draw_cut_circle(x, y, iscale, bbla)

                    fig_box[ycirc, xcirc] = 1
                    tgrad = np.nanmax(bbla[ycirc, xcirc]) - np.nanmin(bbla[ycirc, xcirc])
                    wmean = np.nanmax(wl[ycirc, xcirc])
                    tt = np.nanmin(bbla[ycirc, xcirc])
                    xxx.append(x)
                    yyy.append(y)
                    ltgrad.append(tgrad)
                    lpower.append(wmean)
                    ltt.append(tt)
                    scal.append(orig)

            figure[np.min(yb):np.max(yb), np.min(xb):np.max(xb)] = fig_box

        da_new = xr.DataArray(figure, coords={'time': time, 'lat': lat, 'lon': lon},
                              dims=['lat', 'lon'])
        da_list.append(da_new)
    cube = xr.concat(da_list, dim='time')

    cube.name = 'power'
    enc = {'power': {'complevel': 5, 'zlib': True}}
    cube.to_netcdf(spath + 'gs_scale_circle_slice.nc', encoding=enc)

    f=plt.figure()
    plt.scatter(ltt, lpower, c=ltgrad, cmap='viridis')




def check_wavelets():

    tpath ='/users/global/cornkle/data/OBS/gridsat/gridsat_netcdf/gs_temp_testslice4.nc'
    wpath = '/users/global/cornkle/data/OBS/gridsat/gridsat_netcdf/gs_power_testslice4.nc'

    tar = xr.open_dataarray(tpath)
    war = xr.open_dataarray(wpath)

    time = tar['time'][1] # [730:1030]

    bla = np.array(tar.sel(time=time).values, dtype=float)
    power = np.array(war.sel(time=time).values, dtype=float)
    bla[bla>=-40] = np.nan
    image = tar.sel(time=time).values

    ## remove small
    image[image >= -40] = 0
    labels, numL = label(image)
    mcs = []
    for gi in np.unique(labels):
        if gi == 0:  # index 0 is always background, ignore!
            continue
        pos = np.where(labels == gi)
        #### only MCSs > 15000km2 with pixel = 64km2
        if np.sum(labels == gi) < 234:
            image[pos] = 0
            continue

        #print(np.min(bla[pos]), np.max(bla[pos]))
        mcs.append(pos)

    for p in mcs:
        f = plt.figure(figsize=(12,6))
        f.add_subplot(121)
        bla[np.isnan(bla)] = -40
        gplot = bla
        dat = gplot[np.min(p[0]):np.max(p[0]), np.min(p[1]):np.max(p[1])]
        scale = power[0, :, :]

        sc = scale[np.min(p[0]):np.max(p[0]), np.min(p[1]):np.max(p[1])]
        if np.sum(sc) == 0:
            continue
        plt.imshow(dat) # levels=np.linspace(np.min(dat), np.max(dat), 40)
        plt.colorbar()

        plt.contour(sc, levels=np.linspace(np.min(sc), np.max(sc), 10), cmap='jet')
        plt.colorbar()
        f.add_subplot(122)
        plt.imshow(dat)#, levels=np.linspace(np.min(dat), np.max(dat), 40), cmap='jet')

def check_wavelets_circle():

    tpath ='/users/global/cornkle/data/OBS/gridsat/gridsat_netcdf/gs_temp_circle_slice.nc'
    wpath = '/users/global/cornkle/data/OBS/gridsat/gridsat_netcdf/gs_scale_circle_slice.nc'

    tar = xr.open_dataarray(tpath)
    war = xr.open_dataarray(wpath)

    time = tar['time'][1] # [730:1030]

    bla = np.array(tar.sel(time=time).values, dtype=float)
    power = np.array(war.sel(time=time).values, dtype=float)
    bla[bla>=-40] = np.nan
    image = tar.sel(time=time).values

    ## remove small
    image[image >= -40] = 0
    labels, numL = label(image)
    mcs = []
    for gi in np.unique(labels):
        if gi == 0:  # index 0 is always background, ignore!
            continue
        pos = np.where(labels == gi)
        #### only MCSs > 15000km2 with pixel = 64km2
        if np.sum(labels == gi) < 234:
            image[pos] = 0
            continue

        #print(np.min(bla[pos]), np.max(bla[pos]))
        mcs.append(pos)

    for p in mcs:
        f = plt.figure(figsize=(12,6))
        f.add_subplot(121)
        bla[np.isnan(bla)] = -40
        gplot = bla
        dat = gplot[np.min(p[0]):np.max(p[0]), np.min(p[1]):np.max(p[1])]
        scale = power
        scale[np.isnan(scale)]=0

        sc = scale[np.min(p[0]):np.max(p[0]), np.min(p[1]):np.max(p[1])]
        if np.nansum(sc) == 0:
            continue
        plt.imshow(dat) # levels=np.linspace(np.min(dat), np.max(dat), 40)
        plt.colorbar()

        plt.contour(sc, levels=np.linspace(np.min(sc), np.max(sc), 10), cmap='jet')
        plt.colorbar()
        f.add_subplot(122)
        plt.imshow(dat)#, levels=np.linspace(np.min(dat), np.max(dat), 40), cmap='jet')