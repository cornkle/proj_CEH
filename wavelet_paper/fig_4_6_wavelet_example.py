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

    clat = np.min(dic.lat.values) + ((np.max(dic.lat.values) - np.min(dic.lat.values)) * 0.5)
    clon = np.min(dic.lon.values) + ((np.max(dic.lon.values) - np.min(dic.lon.values)) * 0.5)

    if (clon > 28) or (clon < -17.2) or (clat < 4.1):
        return

    figure = np.zeros_like(outt)

    o2 = outt.copy()
    o2[np.isnan(o2)] = perc
    nok = np.where(abs(grad[0]) > 80)
    d = 2
    i = nok[0]
    j = nok[1]

    for ii, jj in zip(i, j):
        kern = o2[ii - d:ii + d + 1, jj - d:jj + d + 1]
        o2[ii - d:ii + d + 1, jj - d:jj + d + 1] = ndimage.gaussian_filter(kern, 3, mode='nearest')

    wav = util.waveletT(o2, 5)

    arr = np.array(wav['scales'], dtype=str)
    arrf = np.array(wav['scales'], dtype=float)

    scale_ind = range(arr.size)

    yp, xp = np.where(outp > 30)

    figure = np.zeros_like(outt)

    wll = wav['t']
    maxs = np.zeros_like(wll)


    yyy=[]
    xxx=[]
    scal=[]
    for nb in scale_ind[::-1]:

        orig = float(arr[nb])
        print(np.round(orig))

        wl = wll[nb, :, :]
        maxout = (
            wl == ndimage.maximum_filter(wl, (5,5), mode='constant', cval=np.amax(wl) + 1))  # (np.round(orig / 5))

        try:
            yy, xx = np.where((maxout == 1) & (outt <= -40)  & (wl >= np.percentile(wl[wl >= 0.5], 90)) & (wl > orig**.5))# & (wl > orig**.5) )) #(wl >= np.percentile(wl[wl >= 0.5], 90)))# & (wl > orig**.5))#& (wl >= np.percentile(wl[wl >= 0.5], 90))) #)& (wl > orig**.5) (wl >= np.percentile(wl[wl >= 0.1], 90)) )#(wl > orig**.5))#  & (wlperc > orig**.5))# & (wlperc > np.percentile(wlperc[wlperc>=0.1], 80)))# & (wlperc > np.percentile(wlperc[wlperc>=0.1], 80) ))  # & (wl100 > 5)
        except IndexError:
            continue

        for y, x in zip(yy, xx):

            ss = orig
            iscale = (np.ceil(ss / 2. / 5.)).astype(int)
            if ss <= 20:
                iscale = iscale+1

            ycirc, xcirc = ua.draw_cut_circle(x, y, iscale, outp)

            figure[ycirc, xcirc] = np.round(orig)
            xxx.append(x)
            yyy.append(y)
            scal.append(orig)

    figure[np.isnan(outt)]=0

    ##file 130!!! nR
    spos = np.where(np.array(scal, dtype=int) == 15)
    figure[figure == 0] = np.nan


    f = plt.figure(figsize = (7,6), dpi=300)
    ax2 = f.add_subplot(111)
    # ax2.autoscale = False
    ttest = outplot.copy()
    ttest = ttest + 1
    ttest[np.isnan(ttest)] = 0
    ax2.contourf(np.arange(wll.shape[2]) ,  np.arange(wll.shape[1]) , outplot, cmap='Blues', zorder=1)

    ax2.imshow(outt, cmap='Greys', vmax=-40, zorder=2)
    mt = ax2.imshow(figure, cmap='OrRd_r', vmax=180, zorder=3)
    ax2.contour(np.arange(wll.shape[2]), np.arange(wll.shape[1]), ttest, cmap='Blues', levels=[-0.5, 0.5], zorder=4)
    plt.plot(np.array(xxx)[spos], np.array(yyy)[spos], 'wo', markersize=3, label='Wavelet power maximum', zorder=5)

    ax2.invert_yaxis()
    ax2.set_xlim(20, 140)
    ax2.set_ylim(20, 140)
    ax2.set_xticklabels(np.arange(100, 800, 100)-100)
    ax2.set_yticklabels(np.arange(100, 800, 100)-100)
    ax2.plot(xp , yp , 'o', markersize=3, label='Rain > 30mm h$^{-1}$', zorder=10)
    ax2.set_xlabel('Location (km)')
    ax2.set_ylabel('Location (km)')
    plt.colorbar(mt, label = 'Scale (km)')

    plt.tight_layout()
    plt.show()
    spath = '/users/global/cornkle/C_paper/wavelet/figs/paper/'
    plt.savefig(spath + '/method_circles2.png', dpi=300)


    #bla = wcno.file_loop(files[130])
    f = plt.figure(figsize = (6.5,11), dpi=300)

    gridspec.GridSpec(4,1)
    posi = 116
    ax1 = plt.subplot2grid((4,1),(0,0),rowspan=2)
    ax2 = plt.subplot2grid((4,1),(2,0))
    ax3 = plt.subplot2grid((4,1),(3,0))

    lev = np.arange(-90,-39,4)

    ax1.contourf(np.arange(wll.shape[2]) * 5, np.arange(wll.shape[1]) * 5, outplot, cmap='Blues')
    mt = ax1.contourf(np.arange(wll.shape[2])*5,np.arange(wll.shape[1])*5,outt, cmap='Greys', vmax=-40, levels = lev)
    ax1.plot(np.arange(wll.shape[2])*5, [posi*5]*len(np.arange(wll.shape[2])*5), linestyle = '--', linewidth=2, color = 'black')

    ttest = outplot.copy()
    ttest = ttest+1
    ttest[np.isnan(ttest)] = 0
    ax1.contour(np.arange(wll.shape[2]) * 5, np.arange(wll.shape[1]) * 5, ttest, cmap='Blues' , levels=[-0.5,0.5])

    ax1.invert_yaxis()
    ax1.set_xlim(100,700)
    ax1.set_ylim(100, 700)

    ax1.set_xticklabels(np.arange(100, 800, 100) - 100)
    ax1.set_yticklabels(np.arange(100, 800, 100) - 100)

    ax1.plot(xp*5, yp*5, 'o', markersize=3, label='Rain > 30mm h$^{-1}$')
    ax1.legend(loc=4)
    ax1.set_ylabel('Location (km)')
    ax1.set_title(str(dic['time.year'].values)+'-'+str(dic['time.month'].values)+'-'+str(dic['time.day'].values)+' '+str(dic['time.hour'].values)+':'+str(dic['time.minute'].values)+'UTC')

    colors = cm.viridis(np.linspace(0, 1, len([0,1, 2,5,10,20,40,60,80])))

    ax2.plot(np.arange(wll.shape[2])*5, outt[posi,:], color='r')  #118
    ax2.set_xlim(100, 700)
    ax2.set_xticklabels(np.arange(100, 800, 100) - 100)


    ax2.set_ylabel('Cloud-top temperature ($^{\circ}$C)')
    ax22 = ax2.twinx()
    ax22.set_xlim(100, 700)
    ax22.plot(np.arange(wll.shape[2])*5, outp[posi,:])
    ax22.set_ylabel('Rain (mm h$^{-1}$)')
    print(np.nanmax(wll[:,posi,:]))

    mp = ax3.contourf(np.arange(wll.shape[2])*5, arr,wll[:,posi,:], levels=[0,1, 2,5,10,20,40,80,100], colors=colors)
    maxs = np.mean(maxs[:,posi-1:posi+2, :], 1) # -1, +2

    ppos = np.where(maxs)

    for p1, p2 in zip(ppos[1], ppos[0]):
        ax3.errorbar((np.arange(wll.shape[2])*5)[p1], arrf[p2], xerr=arrf[p2]/2, fmt='o', ecolor='white', color='white', capthick=3, ms=3, elinewidth=0.7)
    ax3.set_xlim(100,700)
    ax3.set_xticklabels(np.arange(100, 800, 100) - 100)

    ax3.set_ylim(15, 180)
    ax3.set_xlabel('Location (km)')
    ax3.set_ylabel('Length scale (km)')

    plt.tight_layout()

    f.subplots_adjust(right=0.86)

    cax = f.add_axes([0.87, 0.545, 0.025, 0.415])
    cb = plt.colorbar(mt, cax=cax, label='Cloud-top temperature ($^{\circ}$C)')
    cb.ax.tick_params(labelsize=12)

    cax = f.add_axes([0.87, 0.065, 0.025, 0.175])
    cb = plt.colorbar(mp, cax=cax, label='Wavelet power')
    cb.ax.tick_params(labelsize=12)

    fsiz = 14
    x = 0.02
    plt.annotate('a)', xy=(x, 0.96), xytext=(0, 4), size=fsiz, xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points')
    plt.annotate('b)', xy=(x, 0.51), xytext=(0, 4), size=fsiz, xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points')
    plt.annotate('c)', xy=(x, 0.245), xytext=(0, 4), size=fsiz, xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points')

    plt.show()
    spath = '/users/global/cornkle/C_paper/wavelet/figs/paper/'
    plt.savefig(spath+'/method2.png', dpi=300)

    dic.close()

    plt.close('all')


if __name__ == "__main__":
    files = ua.locate(".nc", '/users/global/cornkle/MCSfiles/WA15_big_-40_15W-20E_size_zR/')
    run(files[130])