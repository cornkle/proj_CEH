# -*- coding: utf-8 -*-
import numpy as np
from scipy import ndimage
from scipy.ndimage.measurements import label


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

    t_thresh_size = -40
    t_thresh_cut = -50

    outt[outt >= t_thresh_size] = 0
    outt[np.isnan(outt)] = 0

    labels, numL = label(outt)

    u, inv = np.unique(labels, return_inverse=True)
    n = np.bincount(inv)

    pix_nb = 700/data_resolution**2

    badinds = u[(n < pix_nb)]
    # all blobs with more than 1000 pixels = 25,000km2 (meteosat regridded 5km), 200pix = 5000km2, 8pix = 200km2
    # scale 30km, radius 15km ca. 700km2 circular area equals 28 pix

    for bi in badinds:
        inds = np.where(labels == bi)
        outt[inds] = 0

    outt[outt >= t_thresh_cut] = 150

    grad = np.gradient(outt)
    outt[outt == 150] = np.nan

    nogood = np.isnan(outt)  # filters edge maxima later, no maxima in -40 edge area by definition!

    # tdiff = np.nanmax(outt) - np.nanmin(outt)  # define background temperature for image
    # if tdiff > 28:  # temp difference of 28 degrees
    #     xmin = 15
    # else:
    #     xmin = 10

    xmin = 10
    outt[nogood] = t_thresh_cut - xmin
    nok = np.where(abs(grad[0]) > 80)
    d = 2
    i = nok[0]
    j = nok[1]
    # edge smoothing for wavelet application
    for ii, jj in zip(i, j):
        kern = outt[ii - d:ii + d + 1, jj - d:jj + d + 1]
        outt[ii - d:ii + d + 1, jj - d:jj + d + 1] = ndimage.gaussian_filter(kern, 3, mode='nearest')

    return outt, nogood, t_thresh_size, t_thresh_cut, pix_nb



def find_scales_dominant(wav, no_good, dataset=None):
    """
    This routine sums up power values of all available scales and identifies areas of dominant power.
    :param wav: wavelet dictionary, output from standard wavelet routine
    :param core_min:
    :param no_good: mask indicating cloud areas that are accepted for dominant power detection
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

    wll = wav['t']

    power_img = np.sum(wll, axis=0)
    power_img[no_good] = 0


    smaller = dataset_dic[dataset]
    thresh_p = np.sum((wav['scales'] + smaller) ** .5) # set different power threshold adjustments to datasets
    try:
        power_img[(power_img < np.percentile(power_img[power_img > 1], 25)) | (power_img < (thresh_p))] = 0
    except IndexError:
        return

    labels, numL = label(power_img)
    u, inv = np.unique(labels, return_inverse=True)

    for inds in u:
        if inds == 0:
            continue

        arr = power_img.copy()
        arr[np.where(labels != inds)] = 0
        power_img.flat[np.argmax(arr)] = -999

    return power_img