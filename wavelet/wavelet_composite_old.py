
import numpy as np
import xarray as xr
from wavelet import util
from utils import u_arrays as ua
import pickle as pkl
from scipy import ndimage
import matplotlib.pyplot as plt
from eod import tm_utils
from collections import defaultdict
import multiprocessing


def composite():
    pool = multiprocessing.Pool(processes=6)
    files = ua.locate(".nc", '/users/global/cornkle/MCSfiles/')
    print('Nb files', len(files))
    dic = xr.open_dataset(files[0])
    dummy = util.waveletTP(dic['t_lag0'].values, dic['p'].values, 5)

    comp_collect = defaultdict(list)

    res = pool.map(file_loop, files)
    pool.close()

    res = [item for sublist in res for item in sublist] # flatten list of lists

    for v in res:
        comp_collect[v[0]].append(v[1])

    keys = comp_collect.keys()
    print(keys)
    xp = (np.arange(1, 42) - 20)*5
    f = plt.figure()
    siz = 3

    ######### 2d plots
    for ind, k in enumerate(keys):
        num = len(keys)
        arr = np.array(comp_collect[k])

        bla = np.nanmean(arr, 0)
        isbigger = arr >= 30
        isfinite = np.isfinite(arr)
        blab = np.sum(isbigger, 0) / np.sum(isfinite,0)*100

        ax = f.add_subplot(3, num, 1 + ind)
        plt.imshow(bla, vmin=0, vmax=3, cmap='viridis')
        plt.title(str(k) + ' km', fontsize=9)
        plt.plot(20, 20, 'ro', markersize=siz)
        cbar = plt.colorbar()
        cbar.set_label('mm h-1')

        ax = f.add_subplot(3, num, 6 + ind)
        plt.imshow(bla, cmap='viridis')
        plt.title(str(k) + ' km', fontsize=9)
        plt.plot(20, 20, 'ro', markersize=siz)
        cbar = plt.colorbar()
        cbar.set_label('mm h-1')

        ax = f.add_subplot(3, num, 11 + ind)
        plt.imshow(blab, vmin=0, vmax=3, cmap='viridis')
        plt.title(str(k) + ' km', fontsize=9)
        plt.plot(20, 20, 'ro', markersize=siz)
        cbar = plt.colorbar()
        cbar.set_label('P(>30mm) %')

    ############# start cross mean
    pplot = []
    for ind, k in enumerate(keys):
        arr = np.array(comp_collect[k])

        bla = np.mean(np.nanmean(arr, 0), 0)
        blab = np.mean(np.nanmean(arr, 0), 1)
        pplot.append((k, bla, blab))

    f = plt.figure()
    col = ['r', 'b', 'g', 'y', 'black']
    ax = f.add_subplot(1, 2, 1)
    for p, c in zip(pplot, col):
        plt.plot(xp, p[1], c, markersize=siz, label=str(p[0]))
    plt.title('West-East')
    ax.vlines(0, 0, 4, linestyles='dashed')

    ax = f.add_subplot(1, 2, 2)
    for p, c in zip(pplot, col):
        plt.plot((np.arange(1, 42) - 20) * 5, p[2], c, markersize=siz, label=str(p[0]))
    ax.vlines(0, 0, 4, linestyles='dashed')
    plt.title('North-South')
    plt.legend()

    ############   start cross probability 5-5 km surrounding point
    pplot = []
    for ind, k in enumerate(keys):
        arr = np.array(comp_collect[k])

        isbigger = arr >= 30
        isfinite = np.isfinite(arr)
        ib = np.sum(isbigger, 0) / np.sum(isfinite, 0)*100

        bla = np.mean(ib[19:21,:], 0)
        blab = np.mean(ib[:,19:21], 1)
        pplot.append((k, bla, blab))

    f = plt.figure()
    col = ['r', 'b', 'g', 'y', 'black']
    ax = f.add_subplot(1, 2, 1)
    for p, c in zip(pplot, col):
        plt.plot((np.arange(1, 42) - 20) * 5, p[1], c, markersize=siz, label=str(p[0]))
    ax.vlines(0, 0, 5, linestyles='dashed')
    plt.title('West-East: >30 5-5')
    # ax.vlines(0, 3, 6, linestyles='dashed')

    ax = f.add_subplot(1, 2, 2)
    for p, c in zip(pplot, col):
        plt.plot((np.arange(1, 42) - 20) * 5, p[2], c, markersize=siz, label=str(p[0]))
    ax.vlines(0, 0, 5, linestyles='dashed')
    plt.title('North-South: >30 5-5')
    plt.legend()

    ############################ start cross probability mean
    pplot = []
    for ind, k in enumerate(keys):
        num = len(keys)
        arr = np.array(comp_collect[k])

        isbigger = arr >= 35
        isfinite = np.isfinite(arr)
        ib = np.sum(isbigger, 0) / np.sum(isfinite, 0) * 100

        bla = np.mean(ib, 0)
        blab = np.mean(ib, 1)
        pplot.append((k, bla, blab))

    f = plt.figure()
    col = ['r', 'b', 'g', 'y', 'black']
    ax = f.add_subplot(1, 2, 1)
    for p, c in zip(pplot, col):
        plt.plot((np.arange(1, 42) - 20) * 5, p[1], c, markersize=siz, label=str(p[0]))
    plt.title('West-East: >30')
    ax.vlines(0, 0,2, linestyles='dashed')

    ax = f.add_subplot(1, 2, 2)
    for p, c in zip(pplot, col):
        plt.plot((np.arange(1, 42) - 20) * 5, p[2], c, markersize=siz, label=str(p[0]))
    ax.vlines(0, 0, 2, linestyles='dashed')
    plt.title('North-South: >30')
    plt.legend()



def file_loop(f):

    scale_ind = [0,12,24,-1,-13]

    ret = []

    print('Doing file: ' + f)
    dic = xr.open_dataset(f)

    outt = dic['tc_lag0'].values
    p = dic['p'].values

    outt[np.isnan(outt)] = 150
    outt[outt >= -40] = 150
    grad = np.gradient(outt)
    outt[outt >= -40] =  np.percentile(outt, 30)
    o2 = outt.copy()
    nok = np.where(abs(grad[0]) > 80)
    d = 2
    i = nok[0]
    j = nok[1]

    for ii, jj in zip(i, j):
        kern = o2[ii - d:ii + d + 1, jj - d:jj + d + 1]
        o2[ii - d:ii + d + 1, jj - d:jj + d + 1] = ndimage.gaussian_filter(kern, 3, mode='nearest')

    wav = util.waveletT(o2, 5)
    wl100 = wav['t'][24, :, :]  # 24 is 60km

    arr = np.array(wav['scales'], dtype=str)

    for nb in scale_ind:
        orig = float(arr[nb])
        scale = np.round(orig)

        wl = wav['t'][nb, :, :]
        maxoutt = (wl == ndimage.maximum_filter(wl,np.round(orig / 5) , mode='constant', cval=np.amax(wl) + 1))  # scale/5

     #   maxoutt = (outt == ndimage.minimum_filter(outt, np.round(orig / 5) , mode='constant', cval=np.amax(outt) + 1))
        maxoutt = maxoutt.astype(int)
        yy, xx = np.where((maxoutt == 1) & (outt < -60) & (wl >= np.percentile(wl, 90)))


        for y, x in zip(yy, xx):

            r = 20
            kernel = tm_utils.cut_kernel(p, x, y, r)

            if kernel.shape != (r * 2 + 1, r * 2 + 1):
                continue

            ret.append((scale,kernel))

    return ret




if __name__ == "__main__":
    #composite()
    llists = composite([0, 12, 24, -1, -13])