# -*- coding: utf-8 -*-
import numpy as np
from scipy import ndimage
from scipy.ndimage.measurements import label
import ipdb
import matplotlib.pyplot as plt

def find_power_sum(coreObj):
    """
    This routine sums up power values of all available scales and identifies areas of dominant power.
    :param wav: wavelet dictionary, output from standard wavelet routine
    :param no_good: mask indicating cloud areas that are accepted for dominant power detection
    :param area: 2d array indicating the number of pixels per MCS
    :return: 2d array of dominant power areas, negative values indicate max power centres (-999)
    The power values for different datasets are not directly comparable. They would have to be normalised.
    Can directly used for frequency analysis though.
    """

    power_img = coreObj.power
    for inds in range(power_img.shape[0]):


        slice = power_img[inds,:,:]
        if coreObj.scales[inds] > 35:
            tsp = 5
        else:
            tsp = 2
        thresh_p = (coreObj.scales[inds]) ** .5  * tsp

        try:
            slice[(slice < np.percentile(slice[slice > 1], 25)) | (slice < (thresh_p))] = 0
        except IndexError:
            return slice * 0

    power_img = np.sum(power_img, axis=0)
    try:
        power_img[coreObj.invalid] = 0
    except:
        ipdb.set_trace()


    labels, numL = label(power_img)

    u, cnt = np.unique(labels, return_counts=True)

    for inds, c in zip(u, cnt):
        if inds == 0:
            continue

        # remove cores that don't reach minimum wavelet scale representing noise.
        if c*coreObj.res**2 < (np.pi * (int(coreObj.scales[0])**2))/4:
            power_img[np.where(labels==inds)] = 0
            continue

        arr = coreObj.image.copy()
        arr[np.where(labels != inds)] = 0
        pos = np.argmin(arr)
        power_img.flat[pos] = coreObj.area.flat[pos] * (-1)  # centre at minT location

        del arr

    return power_img



def find_power_individual(coreObj):
    """
    This routine sums up power values of all available scales and identifies areas of dominant power.
    :param wav: wavelet dictionary, output from standard wavelet routine
    :param no_good: mask indicating cloud areas that are accepted for dominant power detection
    :param area: 2d array indicating the number of pixels per MCS
    :return: 2d array of dominant power areas, negative values indicate max power centres (-999)
    The power values for different datasets are not directly comparable. They would have to be normalised.
    Can directly used for frequency analysis though.
    """


    wll = coreObj.power

    for ids in range(wll.shape[0]):
        arr = wll[ids, :, :]
        out = arr / np.std(arr)
        wll[ids, :, :] = out

    small = (coreObj.scales>12) & (coreObj.scales<=35)
    medium = (coreObj.scales>35) & (coreObj.scales<65)
    large = coreObj.scales>=65

    psmall = np.sum(wll[small,:,:], axis=0)
    pmed = np.sum(wll[medium,:,:], axis=0)
    plarge = np.sum(wll[large,:,:], axis=0)

    scalist = [coreObj.scales[small], coreObj.scales[medium], coreObj.scales[large]]

    thresh_ls = np.sum(coreObj.scales[large]) ** .5  *0.9#* len(coreObj.scales[large])*1.5
    thresh_ss = np.sum(coreObj.scales[small]) ** .5 * len(coreObj.scales[small])*0.75 #*20
    thresh_sl = np.sum(coreObj.scales[small]) ** .5 * 0.01

    thresh_mm = np.sum(coreObj.scales[medium]) ** .5 * len(coreObj.scales[medium])*0.75 #*20
    thresh_ml = np.sum(coreObj.scales[medium]) ** .5 * 0.01

    #maskout =  np.sum(coreObj.power[0:np.min(pos), :, :], axis=0) < 1.5*nbl #np.sum(coreObj.power, axis=0)*0.025#1  #np.sum(coreObj.power[0:np.min(pos), :, :]<1, axis=0) > 0.8*nbl
    ls = (plarge > thresh_ls)
    ss = (psmall > thresh_ss)
    sl = (psmall > thresh_sl)

    mm = (pmed > thresh_mm)
    ml = (pmed > thresh_ml)

    mask = (ls & sl) | ss  #ss |
    med_mask = (ls & ml) | mm

    try:
        psmall[((psmall < np.percentile(psmall[psmall > 1], 25)) | ~mask)] = 0  #| (power_img < (thresh_p))
        psmall[coreObj.invalid] = 0
    except IndexError:
        return


    try:
        pmed[((pmed < np.percentile(pmed[pmed > 1], 25)) | ~med_mask)] = 0  #| (power_img < (thresh_p))
        pmed[coreObj.invalid] = 0
    except IndexError:
        return


    try:
        plarge[((plarge < np.percentile(plarge[plarge > 1], 25)) | ~ls)] = 0  #| (power_img < (thresh_p))
        plarge[coreObj.invalid] = 0
    except IndexError:
        return


    power_img = np.stack([psmall, pmed, plarge], axis=0)

    for idds, pi in enumerate(power_img):
        try:
            pi[coreObj.invalid] = 0
        except:
            ipdb.set_trace()

        labels, numL = label(pi)

        u, cnt = np.unique(labels, return_counts=True)

        for inds, c in zip(u, cnt):
            if inds == 0:
                continue


            # remove cores that don't reach minimum wavelet scale representing noise.
            if c*coreObj.res**2 < (np.pi * (int(np.min(scalist[idds]))**2))/4:
                pi[np.where(labels==inds)] = 0
                continue

            # if c*coreObj.res**2 > (np.pi * (int(np.max(scalist[idds]))**2))/4:
            #     pi[np.where(labels==inds)] = 0
            #     continue

            arr = coreObj.image.copy()
            arr[np.where(labels != inds)] = 0
            pos = np.argmin(arr)
            pi.flat[pos] = coreObj.area.flat[pos] * (-1)

            del arr

    return power_img



def find_power_dominant(coreObj):
    """
    This routine identifies dominant scales and identifies areas of dominant power.
    :param wav: wavelet dictionary, output from standard wavelet routine
    :param no_good: mask indicating cloud areas that are accepted for dominant power detection
    :param area: 2d array indicating the number of pixels per MCS
    :param dataset: string to define input dataset for threshold setting
    :return: 2d array of dominant power areas, negative values indicate max power centres (-999)
    The power values for different datasets are not directly comparable. They would have to be normalised.
    Can directly used for frequency analysis though.
    """


    ss = coreObj.scales<150
    wll = coreObj.power[ss,:,:]
    scales = coreObj.scales[ss]



    for ids in range(wll.shape[0]):
        arr = wll[ids, :, :]
        out = arr / np.std(arr)
        wll[ids, :, :] = out

    power_img = np.sum(wll, axis=0)*0


    maxoutt = (
            wll == ndimage.maximum_filter(wll, (len(scales), coreObj.res, coreObj.res), mode='reflect',  # 5,4,4
                                          cval=np.amax(wll) + 1))  # (np.round(orig / 5))

    for nb in range(scales.size): #[::-1]

        orig = float(scales[nb])

        scale = int(np.round(orig))

        print(np.round(orig))

        wl = wll[nb, :, :]
        wl[coreObj.invalid] = 0
        maxout = maxoutt[nb, :, :]

        vals = wl[wl > 0.01]

        try:
            yy, xx = np.where((maxout == 1)  & (wl > 90/orig*  np.std(vals)))  # &  (wl > orig**.5)) #  #  &

        except IndexError:
            continue



        #wl[(wl < orig ** .5) | (wl - np.mean(vals) < np.std(vals))] = 0#1.3*
        wl[(wl < 1) | (wl - np.mean(vals) <  np.std(vals))] = 0  # 1.3*
        labels, numL = label(wl)

        for y, x in zip(yy, xx):

            inds = labels[y,x]

            if inds == 0:
                continue

            # remove cores that don't reach minimum wavelet scale representing noise.
            if np.sum(labels == inds) * coreObj.res ** 2 < (np.pi * (int(scales[0]) ** 2)) / 4:
                power_img[np.where(labels == inds)] = 0
                continue

            # arr = t.copy()
            # arr[np.where(labels != inds)] = 0
            #power_img[arr<0] = scale
            power_img[labels==inds] = wl[labels==inds] #scale

            arr = coreObj.image.copy() # wl.copy()
            arr[np.where(labels != inds)] = 0
            pos = np.argmin(arr)
            power_img.flat[pos] = scale * (-1) #coreObj.area.flat[pos] * (-1)



    # labels, numL = label(power_img)
    # inds = np.unique(labels)
    # #
    # for ind in inds:
    #     if ind == 0:
    #         continue
    #     arr = coreObj.image.copy()
    #     arr[np.where(labels != ind)] = 0
    #     #ipdb.set_trace()
    #     counts = np.bincount(power_img[labels==ind].astype(int))
    #     maj = np.argmax(counts)
    #     power_img[labels==ind] = maj
    #     pos = np.argmin(arr)
    #     power_img.flat[pos] = coreObj.area.flat[pos] * (-1)
    #
    #     del arr

    return power_img



def find_power_nflics(coreObj):
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


    power_img = np.sum(coreObj.power, axis=0)
    power_img[coreObj.invalid] = 0

    try:
        smaller = dataset_dic[coreObj.data_tag]
    except KeyError:
        print('Data tag not found. Please choose from '+str(dataset_dic.keys()))
        return

    thresh_p = np.sum((coreObj.scales + smaller) ** .5) # set different power threshold adjustments to datasets
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
        pos = np.argmax(arr)
        power_img.flat[pos] = coreObj.area.flat[pos]*(-1)

    return power_img




def find_power_nflicsv2(coreObj):
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

    large_scale = 100
    small_scale = 65
    pos = np.where(coreObj.scales > large_scale) # position of 'large scales'
    mpos = np.where(coreObj.scales <= small_scale)
    #nbl = np.sum(coreObj.scales < large_scale) # number of scales above large scale definition

    if np.max(coreObj.scales > large_scale):  # large scale threshold adjustment

        power_img = np.sum(coreObj.power[0:np.max(mpos), :, :], axis=0)

        thresh_ls = np.sum((coreObj.scales[np.min(pos)::])) ** .5 * len(pos[0])*1.5
        thresh_ss = np.sum((coreObj.scales[0:np.max(mpos)])) ** .5 * len(mpos[0])
        thresh_sl = np.sum((coreObj.scales[0:np.max(mpos)])) ** .5 * 0.5

        #maskout =  np.sum(coreObj.power[0:np.min(pos), :, :], axis=0) < 1.5*nbl #np.sum(coreObj.power, axis=0)*0.025#1  #np.sum(coreObj.power[0:np.min(pos), :, :]<1, axis=0) > 0.8*nbl
        ls = (np.sum(coreObj.power[np.min(pos)::, :, :], axis=0) > thresh_ls)
        ss = (np.sum(coreObj.power[0:np.max(mpos), :, :], axis=0) > thresh_ss)
        sl = (np.sum(coreObj.power[0:np.max(mpos), :, :], axis=0) > thresh_sl)

        mask = (ls & sl) | ss  #ss |

    else:
        power_img = np.sum(coreObj.power, axis=0)
        thresh_all = np.sum((coreObj.scales) ** .5) * len(mpos[0])
        mask = power_img > thresh_all


    power_img[coreObj.invalid] = 0

    try:
        power_img[((power_img < np.percentile(power_img[power_img > 1], 25)) | ~mask)] = 0
    except IndexError:
        return

    labels, numL = label(power_img)
    u, inv = np.unique(labels, return_inverse=True)


    for inds in u:
        arr = coreObj.image.copy()
        if inds == 0:
            continue

        if np.sum(labels == inds) * coreObj.res**2 < (np.pi * (int(coreObj.scales[0])**2))/4:
            power_img[np.where(labels==inds)] = 0
            continue

        arr[np.where(labels != inds)] = 0
        pos = np.argmin(arr)
        power_img.flat[pos] = coreObj.area.flat[pos]*(-1)

    del arr

    return power_img




