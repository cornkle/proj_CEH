import seaborn as sns

pal = sns.color_palette('Blues')
sns.set_context("paper", font_scale=1.5)
sns.set_style("ticks")
import pdb

from utils import u_statistics as ug

import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt


def run_scales():
    #  df = pd.read_pickle('/users/global/cornkle/C_paper/wavelet/saves/pandas/3dmax_gt15000.p')

    path = '/users/global/cornkle/C_paper/wavelet/saves/pandas/'
  #  path = 'D://data/wavelet/saves/pandas/'

    df = pkl.load(open(path + '3dmax_gt15000.p', 'rb'))

    hours = np.arange(0, 24, 1)
    center = (np.arange(24)) + 0.5

    scales = np.arange(15, 185, 20)

    print(scales)

    sc = np.array(df['scale'])
    lat = np.array(df['clat'])
    hour = np.array(df['hour'])
    tmin = np.array(df['circle_Tcentre'])
    area = np.array(df['area'])

    arr40 = np.zeros((scales.size - 1, hours.size))
    arr70 = np.zeros((scales.size - 1, hours.size))
    arr2 = np.zeros((scales.size - 1, hours.size))
    print(center)
    print(hours)

    for iind, s in enumerate(scales):
        if iind == 0:
            continue
        diurn40 = []
        diurn70 = []
        for ind, siz in enumerate(hours):
            sums = sc[
                (sc >= scales[iind - 1]) & (sc < s) & (hour == siz) & (lat > 5) & (
                    tmin < -50)].size
            diurn40.append(sums)

            sums = sc[
                (sc >= scales[iind - 1]) & (sc < s) & (hour == siz) & (lat > 5) & (
                    tmin < -70)].size
            diurn70.append(sums)
        print(s, diurn70)

        arr40[iind - 1, :] = diurn40
        arr70[iind - 1, :] = diurn70
    nb = []
    narea = []
    for ind, siz in enumerate(hours):
        sums = np.unique(area[(sc >= scales[iind - 1]) & (sc < s) & (hour == siz) & (
            lat > 4) & (tmin < -50)]).size
        ar = np.sum(area[(sc >= scales[iind - 1]) & (sc < s) & (hour == siz) & (
            lat > 4) & (tmin < -50)])
        nb.append(sums)
        narea.append(ar)

    arr40 = np.transpose(arr40.T / np.sum(arr40, axis=1)) * 100
    arr70 = np.transpose(arr70.T / np.sum(arr70, axis=1)) * 100
    mean40 = np.mean(arr40, axis=0)
    mean70 = np.mean(arr70, axis=0)
    stddev40 = np.std(arr40, axis=1)
    stddev70 = np.std(arr70, axis=1)
    nb = nb / np.sum(nb)
    narea = narea / np.sum(narea)

    f = plt.figure(figsize=(15, 10), dpi=300)

    ax = f.add_subplot(221)
    plt.contourf(center, scales[0:-1] + 25, (arr40), cmap='viridis')
    plt.title('Normalised diurnal cycle of number of power maxima  < -40degC')
    plt.xlabel('Hours of day')
    plt.ylabel('Scales (km)')
    plt.xticks(np.arange(int(min(center) - 0.5), int(max(center) - 0.5 + 1), 5))
    plt.minorticks_on()
    plt.colorbar(label='%')

    ax = f.add_subplot(222)
    plt.contourf(center, scales[0:-1] + 25, (arr70), cmap='viridis')
    plt.title('Normalised diurnal cycle of power maxima  < -70degC')
    plt.xlabel('Hours of day')
    plt.ylabel('Scales (km)')
    plt.xticks(np.arange(int(min(center) - 0.5), int(max(center) - 0.5 + 1), 5))
    plt.minorticks_on()
    plt.colorbar(label='%')

    ax = f.add_subplot(223)
    plt.contourf(center, scales[0:-1] + 25, (arr40 - mean40), cmap='RdBu',
                 levels=np.arange(-3, 3.5, 0.5), extend='both')
    plt.title('Deviation from mean diurnal cycle < -40degC')
    plt.xlabel('Hours of day')
    plt.ylabel('Scales (km)')
    plt.xticks(np.arange(int(min(center) - 0.5), int(max(center) - 0.5 + 1), 5))
    plt.minorticks_on()
    plt.colorbar(label='%')
    #  plt.contourf(center, scales[0:-1] + 25, std40, hatches=['..'], colors='none', levels=[0.5, 1.5])

    ax = f.add_subplot(224)
    plt.contourf(center, scales[0:-1] + 25, (arr70 - mean70), cmap='RdBu', levels=np.arange(-3, 3.5, 0.5),
                 extend='both')
    plt.title('Deviation from mean diurnal cycle < -70degC')
    plt.xlabel('Hours of day')
    plt.ylabel('Scales (km)')
    plt.xticks(np.arange(int(min(center) - 0.5), int(max(center) - 0.5 + 1), 5))
    plt.minorticks_on()
    plt.colorbar(label='%')
    #  plt.contourf(center, scales[0:-1] + 25, std70, hatches=['..'], colors='none', levels=[0.5, 1.5])
    plt.tight_layout()
    # plt.savefig('/users/global/cornkle/C_paper/wavelet/figs/diurnal/scale_contour.png')

    f = plt.figure()
    ax = f.add_subplot(111)
    plt.contourf(center, scales[0:-1] + 25, (arr70) * 100, cmap='viridis')  # , levels=np.arange(-80,81,0.5))
    plt.colorbar()

    # f = plt.figure()
    # ax = f.add_subplot(111)
    # plt.contourf(center, scales[0:-1] + 25, (arr2) * 100, cmap='viridis')  # , levels=np.arange(-80,81,0.5))
    # plt.colorbar()

    #
    f = plt.figure()
    norm = ug.MidPointNorm(midpoint=0)
    ax = f.add_subplot(111)
    plt.imshow(arr70 - mean70, cmap='RdBu', norm=norm)
    plt.colorbar()

    f = plt.figure()
    ax = f.add_subplot(111)
    plt.plot(hours, arr70[6, :], color='b', label='100-125km')
    plt.plot(hours, arr70[0, :], color='r', label='25-30km')
    # plt.plot(hour, stddev70, color='y', label='standarddev')
    plt.plot(hours, mean70, color='orange', label='all scales mean')
    plt.plot(hours, nb, color='g', label='storm number')
    plt.plot(hours, narea, color='c', label='area')
    plt.legend()

    f = plt.figure()
    ax = f.add_subplot(111)
    plt.plot(scales[0:-1] + 25, stddev40, color='b', label='-40')
    plt.plot(scales[0:-1] + 25, stddev70, color='r', label='-70')
    plt.legend()
    plt.show()


def run_tmin():
    path = '/users/global/cornkle/C_paper/wavelet/figs/paper/'
    path = 'D://data/wavelet/saves/pandas/'

    df = pkl.load(open(path + '3dmax_gt15000_no0.5.p', 'rb'))
    hour = np.arange(0, 23, 1)
    center = (np.arange(23)) + 0.5

    scales = np.arange(15, 185, 20)

    arr40 = np.zeros((scales.size - 1, hour.size))

    for iind, s in enumerate(scales):
        if iind == 0:
            continue
        diurn40 = []
        for ind, siz in enumerate(hour):
            sums = np.nanmean(df['tmin'][
                                  (df['scale'] >= scales[iind - 1]) & (df['scale'] < s) & (df['hour'] == siz) & (
                                  df['clat'] > 4) & (
                                      df['tmin'] < -40)])
            if np.isnan(sums):
                sums = -70
            diurn40.append(sums)

        arr40[iind - 1, :] = diurn40
        print(s, diurn40)
    arr40 = np.transpose(arr40.T - np.mean(arr40, axis=1))

    f = plt.figure()  # figsize=(15, 10), dpi=300)

    mean40 = np.mean(arr40, axis=0)

    stddev40 = np.std(arr40, axis=1)

    std40 = np.transpose(np.abs(arr40 - mean40).T > stddev40 * 2)

    f = plt.figure()  # figsize=(15, 10), dpi=300)

    ax = f.add_subplot(121)
    plt.contourf(center, scales[0:-1] + 25, (arr40), cmap='RdBu_r', levels=np.arange(-5, 6, 0.5), extend='both')
    plt.title('Normalised diurnal cycle of Tmin at power maxima')
    plt.xlabel('Hours of day')
    plt.ylabel('Scales (km)')
    plt.colorbar(label='K')

    ax = f.add_subplot(122)
    plt.contourf(center, scales[0:-1] + 25, (arr40 - mean40), cmap='RdBu_r', extend='both',
                 levels=np.arange(-5, 6, 0.5))
    plt.title('Deviation from mean diurnal cycle')
    plt.xlabel('Hours of day')
    plt.ylabel('Scales (km)')
    plt.colorbar(label='K')
    #  plt.contourf(center, scales[0:-1] + 25, std40, hatches=['..'], colors='none', levels=[0.5, 1.5])

    f = plt.figure()
    ax = f.add_subplot(111)
    plt.plot(hour, arr40[6, :], color='b', label='100-125km')
    plt.plot(hour, arr40[0, :], color='r', label='25-30km')
    # plt.plot(hour, stddev70, color='y', label='standarddev')
    plt.plot(hour, mean40, color='orange', label='all scales mean')
    plt.legend()


def run_pcp():
    path = '/users/global/cornkle/C_paper/wavelet/figs/paper/'
    path = 'D://data/wavelet/saves/pandas/'

    df = pkl.load(open(path + '3dmax_gt15000_no0.5.p', 'rb'))
    hours = np.arange(0, 23, 1)
    center = (np.arange(0, 23, 1)) + 0.5

    scales = np.arange(10, 190, 20)

    print(scales)

    sc = np.array(df['scale'])
    lat = np.array(df['clat'])
    hour = np.array(df['hour'])
    tmin = np.array(df['circle_Tcentre'])
    area = np.array(df['area'])
    p = np.array(df['circle_sum'])

    arr40 = np.zeros((scales.size - 1, hours.size))
    arr70 = np.zeros((scales.size - 1, hours.size))
    arr2 = np.zeros((scales.size - 1, hours.size))
    print(center)
    print(hours)

    for iind, s in enumerate(scales):
        if iind == 0:
            continue
        diurn40 = []

        diurn70 = []

        for ind, siz in enumerate(hours):
            sums = np.nansum(p[
                                 (sc >= scales[iind - 1]) & (sc < s) & (hour == siz) & (lat > 5) & (
                                     tmin < -50)])

            diurn40.append(sums)

            sums = np.nansum(p[
                                 (sc >= scales[iind - 1]) & (sc < s) & (hour == siz) & (lat > 5) & (
                                     tmin < -70)])

            diurn70.append(sums)
        print(s, diurn70)

        arr40[iind - 1, :] = diurn40
        arr70[iind - 1, :] = diurn70
    nb = []
    narea = []
    for ind, siz in enumerate(hours):
        sums = np.unique(area[(sc >= scales[iind - 1]) & (sc < s) & (hour == siz) & (
            lat > 4) & (tmin < -50)]).size
        ar = np.sum(area[(sc >= scales[iind - 1]) & (sc < s) & (hour == siz) & (
            lat > 4) & (tmin < -50)])
        nb.append(sums)
        narea.append(ar)

        # pos = np.isnan(arr40)

        # arr40[:, pos[1]]=np.nan

    arr40 = np.transpose(arr40.T / np.sum(arr40, axis=1))
    arr70 = np.transpose(arr70.T / np.sum(arr70, axis=1))
    mean40 = np.mean(arr40, axis=0)
    mean70 = np.mean(arr70, axis=0)

    nb = nb / np.sum(nb)
    narea = narea / np.sum(narea)

    f = plt.figure(figsize=(15, 10), dpi=300)

    ax = f.add_subplot(221)
    plt.contourf(center, scales[0:-1] + 25, (arr40) * 100, cmap='viridis')
    plt.title('Normalised diurnal cycle of number of power maxima  < -40degC')
    plt.xlabel('Hours of day')
    plt.ylabel('Scales (km)')
    plt.xticks(np.arange(int(min(center) - 0.5), int(max(center) - 0.5 + 1), 5))
    plt.minorticks_on()
    plt.colorbar(label='%')

    ax = f.add_subplot(222)
    plt.contourf(center, scales[0:-1] + 25, (arr70) * 100, cmap='viridis')
    plt.title('Normalised diurnal cycle of power maxima  < -70degC')
    plt.xlabel('Hours of day')
    plt.ylabel('Scales (km)')
    plt.xticks(np.arange(int(min(center) - 0.5), int(max(center) - 0.5 + 1), 5))
    plt.minorticks_on()
    plt.colorbar(label='%')

    ax = f.add_subplot(223)
    plt.contourf(center, scales[0:-1] + 25, (arr40 - mean40) * 100, cmap='RdBu',
                 levels=np.arange(-3, 3.5, 0.5), extend='both')
    plt.title('Deviation from mean diurnal cycle < -40degC')
    plt.xlabel('Hours of day')
    plt.ylabel('Scales (km)')
    plt.xticks(np.arange(int(min(center) - 0.5), int(max(center) - 0.5 + 1), 5))
    plt.minorticks_on()
    plt.colorbar(label='%')
    #  plt.contourf(center, scales[0:-1] + 25, std40, hatches=['..'], colors='none', levels=[0.5, 1.5])

    ax = f.add_subplot(224)
    plt.contourf(center, scales[0:-1] + 25, (arr70 - mean70) * 100, cmap='RdBu', levels=np.arange(-3, 3.5, 0.5),
                 extend='both')
    plt.title('Deviation from mean diurnal cycle < -70degC')
    plt.xlabel('Hours of day')
    plt.ylabel('Scales (km)')
    plt.xticks(np.arange(int(min(center) - 0.5), int(max(center) - 0.5 + 1), 5))
    plt.minorticks_on()
    plt.colorbar(label='%')
    #  plt.contourf(center, scales[0:-1] + 25, std70, hatches=['..'], colors='none', levels=[0.5, 1.5])
    plt.tight_layout()
    # plt.savefig('/users/global/cornkle/C_paper/wavelet/figs/diurnal/scale_contour.png')

    f = plt.figure()
    ax = f.add_subplot(111)
    plt.contourf(center, scales[0:-1] + 25, (arr70) * 100, cmap='viridis')  # , levels=np.arange(-80,81,0.5))
    plt.colorbar()

    # f = plt.figure()
    # ax = f.add_subplot(111)
    # plt.contourf(center, scales[0:-1] + 25, (arr2) * 100, cmap='viridis')  # , levels=np.arange(-80,81,0.5))
    # plt.colorbar()

    #
    f = plt.figure()
    norm = ug.MidPointNorm(midpoint=0)
    ax = f.add_subplot(111)
    plt.imshow(arr70 - mean70, cmap='RdBu', norm=norm)
    plt.colorbar()

    f = plt.figure()
    ax = f.add_subplot(111)
    plt.plot(hours, arr40[4, :], color='b', label='100-125km')
    plt.plot(hours, arr40[0, :], color='r', label='25-30km')
    # plt.plot(hour, stddev70, color='y', label='standarddev')
    plt.plot(hours, mean40, color='orange', label='all scales mean')
    plt.plot(hours, nb, color='g', label='storm number')
    plt.plot(hours, narea, color='c', label='area')
    plt.legend()

    plt.show()

def run_pcp_T():
    fpath = '/users/global/cornkle/C_paper/wavelet/figs/paper/'
    path = '/users/global/cornkle/C_paper/wavelet/saves/pandas/'
   # path = 'D://data/wavelet/saves/pandas/'

    df = pkl.load(open(path + '3dmax_gt15000.p', 'rb'))
    hours = np.arange(0, 23, 1)
    center = (np.arange(0, 23, 1)) + 0.5

    scales = np.arange(10, 180, 10)

    print(scales)

    sc = np.array(df['scale'])
    lat = np.array(df['clat'])
    hour = np.array(df['hour'])
    tmin = np.array(df['circle_Tcentre'])
    area = np.array(df['area'])
    p = np.array(df['circle_sum'])

    arr40 = np.zeros((scales.size - 1, hours.size))
    arr70 = np.zeros((scales.size - 1, hours.size))
    arr2 = np.zeros((scales.size - 1, hours.size))
    print(center)
    print(hours)

    for iind, s in enumerate(scales):
        if iind == 0:
            continue
        diurn40 = []

        diurn70 = []

        for ind, siz in enumerate(hours):
            sums = np.nanmean(tmin[
                                  (sc >= scales[iind - 1]) & (sc < s) & (hour == siz) & (
                                      lat > 5) & (
                                      tmin <= -40)])
            if np.isnan(sums):
                sums = -40
            diurn40.append(sums)

            sums = np.nansum(p[
                                 (sc >= scales[iind - 1]) & (sc < s) & (hour == siz) & (lat > 5) & (
                                     tmin < -70)])

            diurn70.append(sums)
        print(s, diurn70)

        arr40[iind - 1, :] = diurn40
        arr70[iind - 1, :] = diurn70
    nb = []
    narea = []
    for ind, siz in enumerate(hours):
        sums = np.unique(area[(sc >= scales[iind - 1]) & (sc < s) & (hour == siz) & (
            lat > 4) & (tmin < -50)]).size
        ar = np.sum(area[(sc >= scales[iind - 1]) & (sc < s) & (hour == siz) & (
            lat > 4) & (tmin < -50)])
        nb.append(sums)
        narea.append(ar)

        # pos = np.isnan(arr40)

        # arr40[:, pos[1]]=np.nan

    arr40 = np.transpose(arr40.T - np.mean(arr40, axis=1))
    mean40 = np.mean(arr40, axis=0)
    arr70 = np.transpose(arr70.T / np.sum(arr70, axis=1))
    mean70 = np.mean(arr70, axis=0)

    nb = nb / np.sum(nb)
    narea = narea / np.sum(narea)

    f = plt.figure(figsize=(15, 10), dpi=300)

    ax = f.add_subplot(221)
    plt.contourf(center, scales[0:-1] + 25, (arr40), cmap='RdBu_r', vmin=-5, vmax=5, levels=np.arange(-5,5.1,0.5))
    plt.title('Normalised diurnal cycle of number of power maxima  < -40degC')
    plt.xlabel('Hours of day')
    plt.ylabel('Scales (km)')
    plt.xticks(np.arange(int(min(center) - 0.5), int(max(center) - 0.5 + 1), 5))
    plt.minorticks_on()
    plt.colorbar(label='%')

    ax = f.add_subplot(222)
    plt.contourf(center, scales[0:-1] + 25, (arr70) * 100, cmap='viridis', levels=np.arange(0,12.1,0.5))
    plt.title('Normalised diurnal cycle of power maxima  < -70degC')
    plt.xlabel('Hours of day')
    plt.ylabel('Scales (km)')
    plt.xticks(np.arange(int(min(center) - 0.5), int(max(center) - 0.5 + 1), 5))
    plt.minorticks_on()
    plt.colorbar(label='%')

    ax = f.add_subplot(223)
    plt.contourf(center, scales[0:-1] + 25, (arr40 - mean40) , cmap='RdBu_r',
                 levels=np.arange(-3, 3.5, 0.5), extend='both')
    plt.title('Deviation from mean diurnal cycle < -40degC')
    plt.xlabel('Hours of day')
    plt.ylabel('Scales (km)')
    plt.xticks(np.arange(int(min(center) - 0.5), int(max(center) - 0.5 + 1), 5))
    plt.minorticks_on()
    plt.colorbar(label='%')
    #  plt.contourf(center, scales[0:-1] + 25, std40, hatches=['..'], colors='none', levels=[0.5, 1.5])

    ax = f.add_subplot(224)
    plt.contourf(center, scales[0:-1] + 25, (arr70 - mean70) * 100, cmap='RdBu', levels=np.arange(-3, 3.5, 0.5),
                 extend='both')
    plt.title('Deviation from mean diurnal cycle < -70degC')
    plt.xlabel('Hours of day')
    plt.ylabel('Scales (km)')
    plt.xticks(np.arange(int(min(center) - 0.5), int(max(center) - 0.5 + 1), 5))
    plt.minorticks_on()
    plt.colorbar(label='%')
    #  plt.contourf(center, scales[0:-1] + 25, std70, hatches=['..'], colors='none', levels=[0.5, 1.5])
    plt.tight_layout()
    # plt.savefig(fpath+'/scale_contour.png')



    plt.show()