import numpy as np
import xarray as xr
from utils import u_arrays as ua
from scipy import ndimage
import matplotlib.pyplot as plt
import multiprocessing


def minmax():

    pool = multiprocessing.Pool(processes=5)
    files = ua.locate(".nc", '/users/global/cornkle/MCSfiles/WA30/')
    print('Nb files', len(files))

    res = pool.map(file_loop, files)
    pool.close()
    #
    res = [item for sublist in res for item in sublist]  # flatten list of lists
    print('test')
    #
    temp = []
    ptemp = []
    grad = []
    pgrad = []

    for v in res:
        if v[0] == 't':
            temp.append(v[1])
            ptemp.append(v[2])
        if v[0] == 'g':
            grad.append(v[1])
            pgrad.append(v[2])

    f = plt.figure()


    siz = 3

    ax = f.add_subplot(1, 2, 1)
    plt.scatter(temp, ptemp)
    plt.title('temp', fontsize=9)

    ax = f.add_subplot(1, 2, 2)
    plt.scatter(grad, pgrad)
    plt.title('grad', fontsize=9)


def file_loop(f):
    print('Doing file: ' + f)
    dic = xr.open_dataset(f)
    res = []
    outt = dic['t_lag0'].values
    outp = dic['p'].values

    for nb in range(5):
        boole = np.isnan(outp)
        outp[boole] = -1000
        gg = np.gradient(outp)
        outp[boole] = np.nan
        outp[abs(gg[1]) > 300] = np.nan
        outp[abs(gg[0]) > 300] = np.nan

    gradi = np.gradient(outt)
    grad = gradi[1]

    maxoutt = (outt == ndimage.minimum_filter(outt, 20, mode='constant', cval=np.amin(outt) - 1))
    maxoutt = maxoutt.astype(int)
    yt, xt = np.where((maxoutt == 1) & (outt < -50))
    tstr = ['t'] * len(yt)

    maxoutp = (outp == ndimage.maximum_filter(outp, 20, mode='constant', cval=np.amax(outp) + 1))
    maxoutp = maxoutp.astype(int)
    yp, xp = np.where((maxoutp == 1) & (outp > 10))

    maxoutg = (grad == ndimage.minimum_filter(grad, 20, mode='constant', cval=np.amin(grad) - 1))
    maxoutg = maxoutg.astype(int)
    yg, xg = np.where((maxoutg == 1) & (grad < -10))
    gstr = ['g'] * len(yg)

    tstr.extend(gstr)
    tglist = tstr
    tgx = [xt, xg]
    tgx = [item for sublist in tgx for item in sublist]
    tgy = [yt, yg]
    tgy = [item for sublist in tgy for item in sublist]
    points = np.array(list(zip(tgy, tgx)))

    for point in zip(yp, xp):
        try:
            pos = ua.closest_point(point, points)
        except ValueError:
            continue

        if tglist[pos] == 't':
            res.append((tglist[pos], outt[tuple(points[pos])], outp[point]))
        if tglist[pos] == 'g':
            res.append((tglist[pos], grad[tuple(points[pos])], outp[point]))
    dic.close()
    return res
