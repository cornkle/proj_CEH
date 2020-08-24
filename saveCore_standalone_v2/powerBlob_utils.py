# -*- coding: utf-8 -*-
import numpy as np
from scipy import ndimage
from scipy.ndimage.measurements import label
import ipdb

def filter_img(inarr, data_resolution):
    """
    Filters clouds of set area threshold and prepares image for wavelet analysis via adjusting background temperature
    and smoothing cloud edges.
    Current default values:
    Cloud size cutoff threshold: -40C
    Cloud edge cutoff threshold: -50C

    :param inarr: numpy array of cloud top temperatures
    :return: cloud top temperatures with smoothed cloud edges and adjusted background temperature and set tresholds
    """
    outt = inarr.copy()
    print('outmin', np.nanmin(outt), np.nanmax(outt))

    t_thresh_size = -40 # cloud max temperature
    t_thresh_cut = -50 # wavelet edge temperature

    outt[outt >= t_thresh_size] = 0
    outt[np.isnan(outt)] = 0

    labels, numL = label(outt)

    u, inv = np.unique(labels, return_inverse=True)
    n = np.bincount(inv)

    # Min points suited for wavelets ~ 3 pixels across + needed padding = min. of 5 pixels across (25km at 5km res)
    min_diameter_cloud = 5 * data_resolution
    # min number of pixels in circular cloud
    pix_nb = (np.pi * min_diameter_cloud**2)/data_resolution**2  # ~ 500km2 cloud = 20 pixel at 5km res

    badinds = u[(n < pix_nb)]
    goodinds = u[n >= pix_nb]

    area_img = np.zeros_like(outt)
    for bi in badinds:
        inds = np.where(labels == bi)
        outt[inds] = 0
    for bi in goodinds:
        inds = np.where(labels==bi)
        area_img[inds]= np.float(len(inds[0]))

    outt[outt >= t_thresh_cut] = 150

    grad = np.gradient(outt)
    outt[outt == 150] = np.nan

    nogood = np.isnan(outt)  # filters edge maxima later, no maxima in -40 edge area by definition!

    # tdiff = np.nanmax(outt) - np.nanmin(outt)  # define background temperature for image
    # if tdiff > 28:  # temp difference of 28 degrees
    #     xmin = 15
    # else:
    #     xmin = 10


    #### setting xmin to something bigger than 0 and edge smoothing is optional, often improves detection of cores
    #### within big storms and reduces pure edge detection. However, it does reduce the number of detected smaller
    #### storms. For 'popcorn detection', test with xmin = 0 and no edge smoothing.

    xmin = 12
    outt[nogood] = t_thresh_cut - xmin
    nok = np.where(abs(grad[0]) > 80)
    d = 2
    i = nok[0]
    j = nok[1]
    # edge smoothing for wavelet application
    for ii, jj in zip(i, j):
        kern = outt[ii - d:ii + d + 1, jj - d:jj + d + 1]
        outt[ii - d:ii + d + 1, jj - d:jj + d + 1] = ndimage.gaussian_filter(kern, 1, mode='nearest')

    return outt, nogood, t_thresh_size, t_thresh_cut, pix_nb, area_img



def find_dominant_power(wav, no_good, area, data_resolution, dataset=None):
    """
    This routine sums up power values of all available scales and identifies areas of dominant power.
    :param wav: wavelet dictionary, output from standard wavelet routine
    :param no_good: mask indicating cloud areas that are accepted for dominant power detection
    :param area: 2d array indicating the number of pixels per MCS
    :param dataset: string to define input dataset for threshold setting
    :return: 2d array of dominant power areas, negative values indicate max power centres (-999)
    The power values for different datasets are not directly comparable. They would have to be normalised.
    Can directly used for frequency analysis though.
    """
    dataset_dic = {

        'MFG' : -17,   # purely empirical, sorry
        'MSG' : -8,
        'GRIDSAT' : -20,
        'neutral' : 0
    }

    wll = wav['power']

    power_img = np.sum(wll, axis=0)
    power_img[no_good] = 0

    smaller = dataset_dic[dataset]
    thresh_p = np.zeros_like(power_img)+np.sum((wav['scales'] + smaller) ** .5) # set different power threshold adjustments to datasets

    #thresh_p[area>(wav['scales'][-3])**2*4] = np.sum((wav['scales']) ** .5+5)
    #ipdb.set_trace()
    try:
        power_img[(power_img < np.percentile(power_img[power_img > 1], 25)) | (power_img < (thresh_p))] = 0
    except IndexError:
        return power_img * 0

    labels, numL = label(power_img)

    u, cnt = np.unique(labels, return_counts=True)

    for inds, c in zip(u, cnt):
        if inds == 0:
            continue



        # remove cores that don't reach minimum wavelet scale representing noise.
        if c*data_resolution**2 < (np.pi * (int(wav['scales'][0])**2))/4:
            power_img[np.where(labels==inds)] = 0
            continue

        arr = power_img.copy()
        arr[np.where(labels != inds)] = 0
        pos = np.argmax(arr)
        power_img.flat[pos] = area.flat[pos]*(-1)  # maximum power in core replaced by MCS area

    return power_img