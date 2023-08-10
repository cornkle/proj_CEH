"""

@author: C. Klein
"""
import numpy as np
from land_wavelet import twod as w2d
import ipdb


class wavelet(object):


    def __init__(self, res, dist, nb, mother2d = w2d.Mexican_hat(), start=None):

        """
        2D continuous wavelet analysis initialisation. This only supports dx == dy.
        Initialisation sets the scales we want to decompose into.
        From Torrence and Compo: Mexican Hat period, in Fourier sense, is 4 * wavelet scale
        :param res: pixel resolution of prospective input data (e.g. in km)
        :param dist: exponential factor for calculation of distance between decomposition scales, check resulting scales!
        :param start: smallest decomposition scale, smallest resolvable scale is 2*res (== 2*dx)
        :param nb: the number of scales the data is decomposed into
        :param mother2d: a wavelet object, by default Mexican hat
        """
        if start:
            s0 = 2 * start / mother2d.flambda()  # user-defined start scale
        else:
            print('No start scale given, set to 2*dx')
            s0 = 4 * res / mother2d.flambda()  # default start scale: 2 * dx

        a = s0 * 2. ** (np.arange(0, nb + 1) * dist)  # The scales in wavelet space ('wavelet scale')
        freqs = 1. / (mother2d.flambda() * a)  # As of Mallat 1999

        period = 1. / freqs # scales of full period, better represent scale at which >0 values correspond to real scale.
        rscales = period/2. # real scale approximation, for detection of dominant scale from max. coefficient / power values

        self.scale_dist = dist  # exponential factor for calculation of distance between decomposition scales, check resulting scales!
        self.scale_start = s0 # smallest decomposition scale
        self.scale_number = nb # the number of scales the data is decomposed into
        self.res = res # pixel resolution (e.g. in km)
        self.scales = rscales # scales in unit of given pixel resolution
        self.wav_scales = a # scales in wavelet space, can be used to normalise power spectrum
        self.period = period # period scales in unit of given pixel resolution



    def calc_coeffs(self, data, le_thresh=None, ge_thresh=None, fill=0, power_normed='stddev'):
        """
        Calculate pos/neg wavelet coefficients and scale-normalised (always positive) wavelet powers
        :param data: 2d array to decompose into scales
        :param le_thresh: less or equal threshold for wavelet coefficients to be filled with fill value
        :param ge_thresh: greater or equal threshold for wavelet coefficients to be filled with fill value
        :param fill:  fill value for greater / less equal option
        :return: wav_coeffs: positive and negative wavelet coefficients
                 norm_power: normalised wavelet power spectrum
        """

        wav_coeffs = w2d.cwt2d(data, self.res, self.wav_scales)

        wav_coeffs_pure = np.real(wav_coeffs.copy())
        if le_thresh!=None:
            wav_coeffs[np.real(wav_coeffs <= le_thresh)] = fill

        if ge_thresh!=None:
            wav_coeffs[np.real(wav_coeffs >= ge_thresh)] = fill

        norm_power = (np.abs(wav_coeffs)) * (np.abs(wav_coeffs))  # squared wavelet coefficients = power
        scale_dummy = np.reshape(self.wav_scales, (len(self.wav_scales), 1, 1))
        if power_normed == 'scale':
            norm_power = norm_power / (scale_dummy * scale_dummy) # Normalized wavelet power spectrum by squared scale
            # Note: Liu et al 2007 JOAT suggest dividing by wavelet scale only - we emphasize small scales more.
        if power_normed == 'stddev':

            for ids in range(norm_power.shape[0]):
                arr = norm_power[ids,:,:]
                out = arr / np.std(arr)
                norm_power[ids,:,:] = out

            for ids in range(wav_coeffs.shape[0]):
                arr = wav_coeffs[ids,:,:]
                out = arr / np.std(arr)
                wav_coeffs[ids,:,:] = out

        if power_normed == 'none':
            norm_power = norm_power

        return wav_coeffs_pure, norm_power



    def calc_coeffs_inverse(self, data, per_scale=False):
        """
        Calculate pos/neg wavelet coefficients and scale-normalised (always positive) wavelet powers
        :param data: 2d array wavelet coefficients

        :return: array of reconstructed input signal, either as sum over all scales (2d) or for each scale (3d)
        """


        variable = w2d.icwt2d(data, self.res, self.wav_scales, individual=per_scale)

        return variable
