import numpy as np

import matplotlib.pyplot as plt
import matplotlib
from utils import u_arrays as ua
import matplotlib.gridspec as gridspec
from scipy import ndimage
import pdb
from wavelet import util
from utils import u_plot



def ellipse_simple():
    matplotlib.rc('xtick', labelsize=10)
    matplotlib.rc('ytick', labelsize=10)

    ellipse = np.zeros((100,100))

    ycirc, xcirc = ua.draw_ellipse(10,10,4, 4)
    yycirc, xxcirc = ua.draw_ellipse(25, 50, 10, 10)


    ellipse[ycirc,xcirc] = 20
    ellipse[yycirc, xxcirc] = 20


    wav = util.waveletLSTA_power(ellipse, 1)
    lab = 'power'
    wll = wav[lab]

    print('Scales', wav['scales'])


    f = plt.figure()
    plt.imshow(ellipse)

    f = plt.figure()
    plt.imshow(wll[0,:,:])

    print('Small scale max', np.max(np.abs(wll[0,:,:])))

    f = plt.figure()
    plt.imshow(wll[1, :, :])
    print('Mid scale max', np.max(np.abs(wll[1, :, :])))

    f = plt.figure()
    plt.imshow(wll[2, :, :])
    print('Large scale max', np.max(np.abs(wll[2, :, :])))

    f = plt.figure()

    plt.contourf(np.arange(wll.shape[1]),wav['scales'],wll[:,60,:])
    plt.colorbar()

    print((np.max(wll[1, :, :])-np.max(wll[0,:,:])) / np.max(wll[1, :, :]))

def ellipse():
    matplotlib.rc('xtick', labelsize=10)
    matplotlib.rc('ytick', labelsize=10)

    ellipse = np.zeros((100,100))

    ycirc, xcirc = ua.draw_ellipse(50,50,10, 10)
    yycirc, xxcirc = ua.draw_ellipse(25, 50, 4, 4)


    ellipse[ycirc,xcirc] = -50
    ellipse[yycirc, xxcirc] = -50


    wav = util.waveletLSTA_power(ellipse, 1)
    lab = 'power'
    wll = wav[lab]

    pdb.set_trace()

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
            wl == ndimage.maximum_filter(wl, (3, 3), mode='constant', cval=np.amax(wl) + 1))  # (np.round(orig / 5))

        try:
            yy, xx = np.where((maxout == 1)  & (wl > orig**0.5 )) #((wl >= np.percentile(wl[wl >= 0.5], 90))
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
            wll == ndimage.maximum_filter(wll, (3, 3,3), mode='constant', cval=np.amax(wl) + 1))  # (np.round(orig / 5))

    zl, yl, xl = np.where((maxout2==1))
    print('All scales', arr)
    print('Scales with max: ', arr[zl])
    print(np.where((maxout2==1)))
    print(arr[np.where((maxout2==1))[0]])


    amax = np.unravel_index(np.argmax(wll), wll.shape)
    print('Totalmax scale', arr[amax[0]])
    print('MaxScale equivalent pixels', arr[amax[0]]/3)
    print('Circle pixel area', len(ycirc)*3)
    print('Pixel equivalent radius', np.sqrt((len(ycirc)*3)/np.pi))
    print('Power max', np.max(wll))


    f = plt.figure(figsize = (6.5,6), dpi=300)
    #
    gridspec.GridSpec(3,1)
    posi = 56#116 ## 118
    ax1 = plt.subplot2grid((3,1),(0,0),rowspan=2)
    ax3 = plt.subplot2grid((3,1),(2,0))
    #
    lev = np.arange(-90,-39,4)
    mt = ax1.contourf(np.arange(wll.shape[2])*3,np.arange(wll.shape[1])*3,ellipse, cmap='Greys')
    ax1.plot(np.arange(wll.shape[2])*3, [posi*3]*len(np.arange(wll.shape[2])*3), linestyle = '--', linewidth=2, color = 'black')
    ax1.invert_yaxis()
    ax1.legend(loc=4)
    ax1.set_ylabel('Spatial extent (km)')

    mp = ax3.contourf(np.arange(wll.shape[2])*3, arr,wll[:,posi,:], cmap='viridis')
    #levels=[0,1, 2,5,10,20,40,80,100, 130, 150, 180, 200, 300,400]
    maxs = np.mean(maxs[:,posi-1:posi+2, :], 1) # -1, +2
    #ax3.contour(np.arange(wll.shape[2])*5, arr,maxs, cmap='Greys_r')

    ppos = np.where(maxs)

    #for p1, p2 in zip(ppos[1], ppos[0]):
    #    ax3.errorbar((np.arange(wll.shape[2])*5)[p1], arr[p2], xerr=arr[p2]/2, fmt='o', ecolor='white', color='white', capthick=3, ms=3, elinewidth=0.7)
    #ax3.set_xlim(100,700)
    #ax3.set_ylim(15, 180)
    ax3.set_xlabel('Spatial extent (km)')
    ax3.set_ylabel('Length scale (km)')

    plt.tight_layout()

    f.subplots_adjust(right=0.86)

    cax = f.add_axes([0.87, 0.545, 0.025, 0.415])
    cb = plt.colorbar(mt, cax=cax, label='Cloud-top temperature ($^{\circ}$C)')
    cb.ax.tick_params(labelsize=12)

    cax = f.add_axes([0.87, 0.065, 0.025, 0.175])
    cb = plt.colorbar(mp, cax=cax, label='Wavelet power')
    cb.ax.tick_params(labelsize=12)

    path = '/users/global/cornkle/figs/wavelet_examples'
    figname = 'cross'
    filetype = 'png'
    u_plot.savefig(path, figname, filetype)

    # f = plt.figure()
    # plt.imshow(ellipse)
    # print(xl,yl)
    # plt.plot(xl,yl, 'ro')
    # plt.show()
    # try:
    #     f = plt.figure()
    #     plt.imshow(dom)
    # except:
    #     pass
    f = plt.figure()
    plt.imshow(wll[:,posi,:])
    plt.gca().invert_yaxis()

    f = plt.figure()
    plt.title('<=' + str(arr[1]))
    plt.imshow(np.sum(wll[0:1,:,:], axis=0))

    plt.show()

    f = plt.figure()
    plt.title('>='+str(arr[-3]))
    plt.imshow(np.sum(wll[-3::], axis=0))

    plt.show()


def circle():
    matplotlib.rc('xtick', labelsize=10)
    matplotlib.rc('ytick', labelsize=10)


    ellipse = np.zeros((100,100))-70
    short =20

    xcirc, ycirc = ua.draw_circle(50,50,short)


    ellipse[ycirc,xcirc] = -80

    wav = util.waveletT(ellipse, dx=1, dist=1/12., start=15, nb=45)#dx=5, dist=0.08,start=15,nb=15 )

    wll = wav['t']
    arr = np.round(wav['scales'])
    print('AVAIL WAVELET SCALES: ', arr)
    maxs = np.zeros_like(wll)
    yl = []
    xl = []
    for nb in range(wav['t'].shape[0]):

        orig = float(wav['scales'][nb])

        wl = wav['t'][nb, :, :]

        maxout = (
            wl == ndimage.maximum_filter(wl, (5, 5), mode='constant', cval=np.amax(wl) + 1))  # (np.round(orig / 5))

        try:
            yy, xx = np.where((maxout == 1) & (wl > orig ** .5))
        except IndexError:
            continue

        for y, x in zip(yy, xx):
            #print(arr[nb],y,x)

            maxs[nb,y,x] = 1
            #print('Power value',wll[nb,y,x])
            yl.append(y)
            xl.append(x)

    print('finish loop')

    maxout2 = (
            wll == ndimage.maximum_filter(wll, (5, 5,5), mode='reflect', cval=np.amax(wl) + 1))  # (np.round(orig / 5))

    zl, yl, xl = np.where((maxout2==1) & (wll > arr.repeat(100*100,axis=0).reshape((46,100,100)) ** .5))
    wlmax = np.max(wll[zl,yl,xl])
    pl = np.where(wll == wlmax)

    zll, yll, xll = np.where((maxs == 1))
    wllmax = np.max(wll[zll,yll,xll])
    pll = np.where(wll == wllmax)

    print('Max point scales: ', arr[zl])

    amax = np.unravel_index(np.argmax(wll), wll.shape)
    print('Totalmax whole domain', arr[amax[0]])
    print('Totalmax 3d', arr[pl[0]])
    print('Totalmax 2d', arr[pll[0]])
    print('Scale of perfect circle', (2*short+1)*wav['res'])
    print('Pixel-adjusted perfect circle', (2 *short + 1 -2) * wav['res'])
    print('Max scale in pixelacross', arr[amax[0]]/wav['res'])
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

    f = plt.figure()
    pos = np.argmin(np.abs((2*short+1)*wav['res']-arr))

    plt.imshow(wll[pos, :,:])
    print(xl,yl)
    plt.plot(xl,yl, 'ro')
    plt.show()
    f = plt.figure()
    plt.imshow(wll[amax[0], :,:])
    print(xl,yl)
    plt.plot(xl,yl, 'ro')
    plt.show()