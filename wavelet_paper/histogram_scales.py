import seaborn as sns
pal = sns.color_palette('Blues')
sns.set_context("paper", font_scale=1.5)
sns.set_style("ticks")
import ipdb
import matplotlib.cm as cm

from utils import u_graphics as ug

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl

def run_pcp():
    dic = pkl.load(open('/users/global/cornkle/C_paper/wavelet/saves/pandas/precip_3dmax_gt15000_-70_15km.p', 'rb'))
    keys=[]
    for k in dic.keys():
        keys.append(k)
    keys=np.sort(keys)

    for k in keys:
     #   dic[k] = [item for sublist in dic[k] for item in sublist]
        arr = np.array(dic[k])
        arr = arr[(np.isfinite(arr)) & (arr>=0.1)]

        dic[k] = arr

    f = plt.figure()
    ax = f.add_subplot(111)

    colors = cm.rainbow(np.linspace(0,1,len(keys)))

    for k,c in zip(keys[::-1], colors):
        weights = np.ones_like(dic[k]) / float(len(dic[k]))
        hist, h = np.histogram(dic[k], bins=np.arange(0.1,100+1,1), weights=weights, range=(0.1,100))
        print(h)

        line, = ax.semilogy(hist, color=c, lw=2, label=str(k))
        plt.ylabel('ln(normalised frequency of non-zero rain)')
        plt.xlabel('rainfall (mm h-1)')
        plt.title('Sub system features of storms <70degC, >15000km2')

    plt.legend()

    # ax = f.add_subplot(122)
    #
    # colors = cm.rainbow(np.linspace(0, 1, len(keys)))
    #
    # for k, c in zip(keys[::-1], colors):
    #     weights = np.ones_like(dic[k]) / float(len(dic[k]))
    #     hist, h = np.histogram(dic[k], bins=np.arange(0.1, 100 + 1, 1), weights=weights, range=(0.1, 100))
    #     print(h)
    #
    #     plt.plot(hist, color=c, lw=2, label=str(k))
    #     plt.ylabel('normalised frequency of non-zero rain')
    #     plt.xlabel('rainfall (mm h-1)')
    #
    # plt.legend()