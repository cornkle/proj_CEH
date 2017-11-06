import numpy as np

import matplotlib.pyplot as plt
import matplotlib
import matplotlib.patches as mpatches
from utils import u_arrays as ua
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from scipy import ndimage
import pdb
from wavelet import util
import time



def ellipse():
    matplotlib.rc('xtick', labelsize=10)
    matplotlib.rc('ytick', labelsize=10)


    ellipse = np.zeros((100,100))+70
    short = 12
    long = 12

    xcirc, ycirc = ua.draw_ellipse(50,50,short, long)


    ellipse[ycirc,xcirc] = 80
    #ellipse[np.arange(50,54), [56]*4] = -74

    wav = util.waveletLSTA(ellipse, 5, method='dry')
    lab = 'power'
    wll = wav[lab]
    arr = np.round(wav['scales'])
    maxs = np.zeros_like(wll)
    yl = []
    xl = []
    for nb in range(wav[lab].shape[0]):

        orig = float(wav['scales'][nb])
        scale = int(np.round(orig))

        print(np.round(orig))

        wl = wav[lab][nb, :, :]
        # maxout = maxoutt[nb, :, :]

        maxout = (
            wl == ndimage.maximum_filter(wl, (5, 5), mode='constant', cval=np.amax(wl) + 1))  # (np.round(orig / 5))

        try:
            yy, xx = np.where((maxout == 1) & ((wl >= np.percentile(wl[wl >= 0.5], 90)) & (wl > orig*15 )))
        except IndexError:
            continue
        #pdb.set_trace()
        for y, x in zip(yy, xx):
            print(arr[nb],y,x)

            maxs[nb,y,x] = 1
            print('Power value',wll[nb,y,x])
            yl.append(y)
            xl.append(x)

    print('finish loop')

    maxout2 = (
            wll == ndimage.maximum_filter(wll, (5, 5,5), mode='constant', cval=np.amax(wl) + 1))  # (np.round(orig / 5))

    zl, yl, xl = np.where((maxout2==1))
    print('Scales: ', arr[zl])
    print(np.where((maxout2==1)))
    print(arr[np.where((maxout2==1))[0]])


    amax = np.unravel_index(np.argmax(wll), wll.shape)
    print('Totalmax', arr[amax[0]])
    print('Pixelacross', arr[amax[0]]/5)
    print('Power max', np.max(wll))

    f = plt.figure(figsize = (6.5,11), dpi=300)
    #
    gridspec.GridSpec(3,1)
    posi = 56#116 ## 118
    ax1 = plt.subplot2grid((3,1),(0,0),rowspan=2)
    ax3 = plt.subplot2grid((3,1),(2,0))
    #
    lev = np.arange(-90,-39,4)
    mt = ax1.contourf(np.arange(wll.shape[2])*5,np.arange(wll.shape[1])*5,ellipse, cmap='Greys')
    ax1.plot(np.arange(wll.shape[2])*5, [posi*5]*len(np.arange(wll.shape[2])*5), linestyle = '--', linewidth=2, color = 'black')
    ax1.invert_yaxis()
    ax1.legend(loc=4)
    ax1.set_ylabel('Spatial extent (km)')
    #colors = cm.viridis(np.linspace(0, 1, len([0,1, 2,5,10,20,40,60,80,100])))
    #
    mp = ax3.contourf(np.arange(wll.shape[2])*5, arr,wll[:,posi,:], levels=[0,1, 2,5,10,20,40,80,100, 130, 150, 180, 200, 300,400], cmap='viridis')
    maxs = np.mean(maxs[:,posi-1:posi+2, :], 1) # -1, +2
    #ax3.contour(np.arange(wll.shape[2])*5, arr,maxs, cmap='Greys_r')

    ppos = np.where(maxs)

    #for p1, p2 in zip(ppos[1], ppos[0]):
    #    ax3.errorbar((np.arange(wll.shape[2])*5)[p1], arr[p2], xerr=arr[p2]/2, fmt='o', ecolor='white', color='white', capthick=3, ms=3, elinewidth=0.7)
    #ax3.set_xlim(100,700)
    ax3.set_ylim(15, 180)
    ax3.set_xlabel('Spatial extent (km)')
    ax3.set_ylabel('Length scale (km)')

    #plt.tight_layout()

    f.subplots_adjust(right=0.86)

    cax = f.add_axes([0.87, 0.545, 0.025, 0.415])
    cb = plt.colorbar(mt, cax=cax, label='Cloud-top temperature ($^{\circ}$C)')
    cb.ax.tick_params(labelsize=12)

    cax = f.add_axes([0.87, 0.065, 0.025, 0.175])
    cb = plt.colorbar(mp, cax=cax, label='Wavelet power')
    cb.ax.tick_params(labelsize=12)

    plt.show()

    f = plt.figure()
    plt.imshow(ellipse)
    print(xl,yl)
    plt.plot(xl,yl, 'ro')
    plt.show()



def circle():
    matplotlib.rc('xtick', labelsize=10)
    matplotlib.rc('ytick', labelsize=10)


    ellipse = np.zeros((100,100))-70
    short = 6

    xcirc, ycirc = ua.draw_circle(50,50,short)


    ellipse[ycirc,xcirc] = -80
    #ellipse[np.arange(50,54), [56]*4] = -74

    wav = util.waveletT(ellipse, 5)

    wll = wav['t']
    arr = np.round(wav['scales'])
    maxs = np.zeros_like(wll)
    yl = []
    xl = []
    for nb in range(wav['t'].shape[0]):

        orig = float(wav['scales'][nb])
        scale = int(np.round(orig))

        print(np.round(orig))

        wl = wav['t'][nb, :, :]
        # maxout = maxoutt[nb, :, :]

        maxout = (
            wl == ndimage.maximum_filter(wl, (5, 5), mode='constant', cval=np.amax(wl) + 1))  # (np.round(orig / 5))

        try:
            yy, xx = np.where((maxout == 1) & ((wl >= np.percentile(wl[wl >= 0.5], 90)) & (wl > orig*15 )))
        except IndexError:
            continue
        #pdb.set_trace()
        for y, x in zip(yy, xx):
            print(arr[nb],y,x)

            maxs[nb,y,x] = 1
            print('Power value',wll[nb,y,x])
            yl.append(y)
            xl.append(x)

    print('finish loop')

    maxout2 = (
            wll == ndimage.maximum_filter(wll, (5, 5,5), mode='constant', cval=np.amax(wl) + 1))  # (np.round(orig / 5))

    zl, yl, xl = np.where((maxout2==1))
    print('Scales: ', arr[zl])
    print(np.where((maxout2==1)))
    print(arr[np.where((maxout2==1))[0]])


    amax = np.unravel_index(np.argmax(wll), wll.shape)
    print('Totalmax', arr[amax[0]])
    print('Pixelacross', arr[amax[0]]/5)
    print('Power max', np.max(wll))

    f = plt.figure(figsize = (6.5,11), dpi=300)
    #
    gridspec.GridSpec(3,1)
    posi = 56#116 ## 118
    ax1 = plt.subplot2grid((3,1),(0,0),rowspan=2)
    ax3 = plt.subplot2grid((3,1),(2,0))
    #
    lev = np.arange(-90,-39,4)
    mt = ax1.contourf(np.arange(wll.shape[2])*5,np.arange(wll.shape[1])*5,ellipse, cmap='Greys', vmax=-40, levels = lev)
    ax1.plot(np.arange(wll.shape[2])*5, [posi*5]*len(np.arange(wll.shape[2])*5), linestyle = '--', linewidth=2, color = 'black')
    ax1.invert_yaxis()
    ax1.legend(loc=4)
    ax1.set_ylabel('Spatial extent (km)')
    #colors = cm.viridis(np.linspace(0, 1, len([0,1, 2,5,10,20,40,60,80,100])))
    #
    mp = ax3.contourf(np.arange(wll.shape[2])*5, arr,wll[:,posi,:], levels=[0,1, 2,5,10,20,40,80,100, 130, 150, 180, 200, 300,400], cmap='viridis')
    maxs = np.mean(maxs[:,posi-1:posi+2, :], 1) # -1, +2
    #ax3.contour(np.arange(wll.shape[2])*5, arr,maxs, cmap='Greys_r')

    ppos = np.where(maxs)

    #for p1, p2 in zip(ppos[1], ppos[0]):
    #    ax3.errorbar((np.arange(wll.shape[2])*5)[p1], arr[p2], xerr=arr[p2]/2, fmt='o', ecolor='white', color='white', capthick=3, ms=3, elinewidth=0.7)
    #ax3.set_xlim(100,700)
    ax3.set_ylim(15, 180)
    ax3.set_xlabel('Spatial extent (km)')
    ax3.set_ylabel('Length scale (km)')

    #plt.tight_layout()

    f.subplots_adjust(right=0.86)

    cax = f.add_axes([0.87, 0.545, 0.025, 0.415])
    cb = plt.colorbar(mt, cax=cax, label='Cloud-top temperature ($^{\circ}$C)')
    cb.ax.tick_params(labelsize=12)

    cax = f.add_axes([0.87, 0.065, 0.025, 0.175])
    cb = plt.colorbar(mp, cax=cax, label='Wavelet power')
    cb.ax.tick_params(labelsize=12)

    plt.show()

    f = plt.figure()
    plt.imshow(ellipse)
    print(xl,yl)
    plt.plot(xl,yl, 'ro')
    plt.show()