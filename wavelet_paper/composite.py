import numpy as np

import matplotlib.pyplot as plt
import matplotlib
import pdb
import seaborn as sns
import pickle as pkl
#pal = sns.color_palette('Blues')
#sns.set_context("paper", font_scale=1.5)
sns.set_style("darkgrid")
matplotlib.rc('xtick', labelsize=8)
matplotlib.rc('ytick', labelsize=8)



out = '/users/global/cornkle/C_paper/wavelet/saves/pandas/'
comp_collect = pkl.load(open(out + 'comp_collect_composite_noR.p','rb'))


siz = 3

keys = comp_collect.keys()
print(keys)
keys = sorted(keys)
keys = np.array(keys)

ranges = [15, 30,  60, 90, 120, 150, 180, 205 ]
outrange = [30,  60, 90, 120, 150, 180, 205 ]

out = []

for id, r in enumerate(ranges):
    if id == 0:
        continue
    precede=ranges[id-1]
    if precede == 15:
        precede = precede -1


    klist = keys[(keys <= r) & (keys>ranges[id-1])]

    print(r, len(klist))
    print(klist)

    bbig = []
    ffin = []
    aarr = []
    taarr = []
    nnz = []

    for k in klist:

        p = np.array(comp_collect[k]['p'])
        t = np.array(comp_collect[k]['t'])

        pos = np.where(t[:, 21, 21] <= -45)
        pos = pos[0]

        big = (p>30)[pos,:,:]
        fin = (np.isfinite(p))[pos,:,:]
        arr = p[pos,:,:]
        nz = (p>0.1)[pos,:,:]
        tarr = t[pos,:,:]

        bbig.append(big)
        ffin.append(fin)
        aarr.append(arr)
        taarr.append(tarr)
        nnz.append(nz)


    try:

        bbig = np.concatenate(bbig, axis=0)
        ffin = np.concatenate(ffin, axis=0)
        aarr = np.concatenate(aarr, axis=0)
        taarr = np.concatenate(taarr, axis=0)
        nnz = np.concatenate(nnz, axis=0)

    except ValueError:
        print('return')



    bla = np.nansum(aarr, 0) / np.nansum(nnz, 0)
    bla1 = np.nansum(ffin, 0)
    blab = np.nansum(bbig, 0) / np.nansum(ffin,0) * 100

    tbla = np.nanmedian(taarr, 0)

    out.append((bla,bla1,blab, tbla))




######### 2d plots


f = plt.figure(figsize=(12, 5), dpi=400)
ll = [30, 60, 90, 120, 180]  # keys

for ind, k in enumerate(ll):

    pos=outrange.index(k)

    outbla = out[pos]

    x = np.linspace(-2, 2, outbla[1].shape[0])
    fos = 12

    print(x)
    print(np.arange(-2,3, 1))

    # ax3 = f.add_subplot(4, 5, 11 + ind)
    # mp3 = ax3.imshow(outbla[0], cmap='viridis', vmin=3, vmax=8)  # vmin=0, vmax=6,
    # plt.title(str(ranges[pos])+'-'+str(k) + ' km', fontsize=fos)
    # ax3.plot(20, 20, 'ro', markersize=siz)
    # ax3.set_xticklabels(np.arange(-3, 3.1, 1))
    # ax3.set_yticklabels(np.arange(-3, 3.1, 1))
    # cbar = plt.colorbar(mp3)
    # cbar.set_label('Average rain (mm h$^{-1}$)', fontsize=fos)

    ax2 = f.add_subplot(2, 5, 6 + ind)
    mp2 = ax2.imshow(outbla[2], aspect=1, cmap='viridis') #x,x,
    #plt.title(str(ranges[pos])+'-'+str(k) + ' km', fontsize=fos)
    ax2.plot(20, 20, 'ro', markersize=siz)
    ax2.set_xticklabels('')
    ax2.set_yticklabels('')
    if ind == 0:
        ax2.set_yticklabels(np.arange(-3,4, 1))
        ax2.set_ylabel('10$^2$ km')

    ax2.set_xticklabels(np.arange(-3,4, 1))
    ax2.set_xlabel('10$^2$ km')


    # ax4 = f.add_subplot(2, 5, 6 + ind)
    # mp4 = ax4.imshow(outbla[2], vmin=0, vmax=4, cmap='viridis')
    # #plt.title(str(ranges[pos])+'-'+str(k) + ' km', fontsize=fos)
    # ax4.plot(20, 20, 'ro', markersize=siz)
    # ax4.set_xticklabels('')
    # ax4.set_yticklabels('')
    # if ind == 0:
    #     ax4.set_yticklabels(np.arange(-3, 4, 1))
    #     ax4.set_ylabel('10$^2$ km')
    # cbar = plt.colorbar(mp4)
    # if ind == 4:
    #     cbar.set_label('Pixel probability (%)', fontsize=fos)

    ax1 = f.add_subplot(2, 5, 1 + ind)
    mp1 = ax1.imshow(outbla[3], vmin=-70, vmax=-55, cmap='inferno', aspect=1)
    ax1.set_xticklabels('')

    plt.title(str(ranges[pos])+'-'+str(k) + ' km', fontsize=13)
    ax1.plot(20, 20, 'ro', markersize=siz)

    ax1.set_yticklabels('')
    if ind == 0:
        ax1.set_yticklabels(np.arange(-3, 4, 1))
        ax1.set_ylabel('10$^2$ km')
plt.tight_layout()

f.subplots_adjust(right=0.91)
cax = f.add_axes([0.92,0.545,0.025,0.37])
cbar = f.colorbar(mp1, cax)
cbar.ax.tick_params(labelsize=12)
cbar.set_label('Cloud top temperature ($\circ$C)', fontsize=12)
cax = f.add_axes([0.92,0.12,0.025,0.37])
cbar = f.colorbar(mp2, cax)
cbar.ax.tick_params(labelsize=12)
cbar.set_label('Pixel probability (%)', fontsize=12)

plt.savefig('/users/global/cornkle/C_paper/wavelet/figs/paper/composite3d_noR.png')
plt.close('all')

col = ['r', 'b', 'g', 'y', 'black']

