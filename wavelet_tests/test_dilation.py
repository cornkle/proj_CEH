import numpy as np
import xarray as xr
from wavelet import util
from scipy import ndimage
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from utils import u_arrays as ua


matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)

def run(fi):
    ret = []

    print('Doing file: ' + fi)

    dic = xr.open_dataset(fi)

    outt = dic['tc_lag0'].values
    outp = dic['p'].values
    outpc = dic['pconv'].values

    outplot = outp.copy()

    outt[np.isnan(outt)] = 150
    outt[outt >= -40] = 150
    grad = np.gradient(outt)
    outt[outt == 150] = np.nan
    outp[np.isnan(outt)] = np.nan
    outpc[np.isnan(outt)] = np.nan


    area = np.sum(outt <= -40)
    try:
        bulk_pmax = np.max(outp[(np.isfinite(outp)) & (np.isfinite(outt))])
    except ValueError:
        return ret
    try:
        bulk_pmin = np.min(outp[(np.isfinite(outp)) & (np.isfinite(outt))])
    except ValueError:
        return ret

    if (area * 25 < 15000) or (area * 25 > 800000)  or (bulk_pmax > 200) or (bulk_pmin < 0):
        print(area*25)
        print('throw out')
        return

    perc = np.percentile(outt[np.isfinite(outt)], 60)  # 60
    print('perc:',perc)
   # perc=-60

    clat = np.min(dic.lat.values) + ((np.max(dic.lat.values) - np.min(dic.lat.values)) * 0.5)
    clon = np.min(dic.lon.values) + ((np.max(dic.lon.values) - np.min(dic.lon.values)) * 0.5)

    if (clon > 28) or (clon < -17.2) or (clat < 4.1):
        return

    figure = np.zeros_like(outt)
    outt[np.isnan(outt)] = perc
    o1 = outt.copy()
    o2 = o1.copy()
    nok = np.where(abs(grad[0]) > 80)
    d = 2
    i = nok[0]
    j = nok[1]

    for ii, jj in zip(i, j):
        kern = outt[ii - d:ii + d + 1, jj - d:jj + d + 1]
        o1[ii - d:ii + d + 1, jj - d:jj + d + 1] = ndimage.gaussian_filter(kern, 3, mode='nearest')
        o2[ii - d:ii + d + 1, jj - d:jj + d + 1] = ndimage.grey_dilation(kern, size=(3,3))


    f = plt.figure(figsize = (7,6), dpi=300)
    ax2 = f.add_subplot(221)
    mp = ax2.imshow(outt)
    ax2.set_xticklabels([])
    plt.title('Original MCS', fontsize=8)
    plt.colorbar(mp)
    ax = f.add_subplot(222)
    mp=ax.imshow(o1)
    plt.title('Gaussian smoothing, sigma=3', fontsize=8)
    plt.colorbar(mp)
    ax3 = f.add_subplot(223)
    plt.title('Greyscale dilation, 3x3 kernel', fontsize=8)
    mp1 = ax3.imshow(o2)
    clb = plt.colorbar(mp1)
    clb.set_label('Cloud top temperature (degC)', fontsize=8)
    #plt.tight_layout()


if __name__ == "__main__":
    files = ua.locate(".nc", '/users/global/cornkle/MCSfiles/WA15_big_-40_15W-20E_size_zR/')
    run(files[238])