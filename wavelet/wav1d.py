"""

@author: C. Klein
"""
import numpy as np
from wavelet import oned as w1d


class wavelet(object):


    def __init__(self, res, dist, nb, mother1d = w1d.Mexican_hat(), start=None):

        """
        1D continuous wavelet analysis initialisation.
        Initialisation sets the scales we want to decompose into.

        :param res: pixel resolution of prospective input data (e.g. in km)
        :param dist: exponential factor for calculation of distance between decomposition scales, check resulting scales!
        :param start: smallest decomposition scale, smallest resolvable scale is 2*res (== 2*dx)
        :param nb: the number of scales the data is decomposed into
        :param mother1d: a wavelet object, by default Mexican hat
        """
        if start:
            s0 = 2 * start / mother1d.flambda()  # user-defined start scale
        else:
            print('No start scale given, set to 2*dx')
            s0 = 4 * res / mother1d.flambda()  # default start scale: 2 * dx

        a = s0 * 2. ** (np.arange(0, nb + 1) * dist)  # The scales in wavelet space ('wavelet scale')
        freqs = 1. / (mother1d.flambda() * a)  # As of Mallat 1999
        period = 1. / freqs
        scales = period/2. #(period/3)*2 # # 'real' scale approximation, alternative: (period/3)*2

        self.scale_dist = dist  # exponential factor for calculation of distance between decomposition scales, check resulting scales!
        self.scale_start = s0 # smallest decomposition scale
        self.scale_number = nb # the number of scales the data is decomposed into
        self.res = res # pixel resolution (e.g. in km)
        self.scales = scales # scales in unit of given pixel resolution
        self.norm_scales = a # wavelet scales for normalising power spectrum
        self.mother = mother1d


    def significance(self, signal3d, sig3d, direction='x', mask=None):

        out = np.zeros((len(self.scales), signal3d.shape[1], signal3d.shape[2]))

        if direction == 'x':
            for row in range(signal3d.shape[1]):

                data = signal3d[:,row,:]
                power = sig3d[:,row,:]

                for ids in range(signal3d.shape[0]):
                    pscale = power[ids,:]
                    if mask is None:
                        inn = pscale
                    else:
                        inn = pscale[mask[row,:]]
                    s = np.std(inn)
                    sigout = data[ids,:] > 10*s
                    out[ids, row, :] = sigout

        if direction == 'y':
            for column in range(signal3d.shape[2]):

                data = signal3d[:,:,column]
                power = sig3d[:,:,column]

                for ids in range(signal3d.shape[0]):
                    pscale = power[ids,:]
                    if mask is None:
                        inn = pscale
                    else:
                        inn = pscale[mask[:,column]]
                    s = np.std(inn)
                    sigout = data[ids, :] > 10*s
                    out[ids, :, column] = sigout


        return out

        # # Calculates the global wavelet spectrum and determines its significance level.
        # glbl_power = std2 * power.mean(axis=1)
        # dof = data.size - self.scales  # Correction for padding at edges
        # glbl_signif, tmp = wvt.significance(std2, self.res, self.scales, 1, alpha,
        #                                             dof=dof, wavelet=self.mother)



    def calc_coeffs(self, data, le_thresh=None, ge_thresh=None, fill=0, direction='x'):
        """
        Calculate pos/neg wavelet coefficients and scale-normalised (always positive) wavelet powers
        :param data: 2d array to decompose into scales
        :param le_thresh: less or equal threshold for wavelet coefficients to be filled with fill value
        :param ge_thresh: greater or equal threshold for wavelet coefficients to be filled with fill value
        :param fill:  fill value
        :return: wav_coeffs: positive and negative wavelet coefficients
                 norm_power: normalised wavelet power spectrum
        """

        wav_coeffs = w1d.cwt1d(data, self.res, dj=self.scale_dist, s0=self.scale_start, J=self.scale_number, direction=direction)


        wav_coeffs_pure = np.real(wav_coeffs.copy())
        if le_thresh!=None:
            wav_coeffs[np.real(wav_coeffs <= le_thresh)] = fill

        if ge_thresh!=None:
            wav_coeffs[np.real(wav_coeffs >= ge_thresh)] = fill

        norm_power = (np.abs(wav_coeffs)) * (np.abs(wav_coeffs))  # squared wavelet coefficients
        scale_dummy = np.reshape(self.norm_scales, (len(self.norm_scales), 1, 1))

        for ids in range(norm_power.shape[0]):
            arr = norm_power[ids,:,:]
            out = arr / np.std(arr)
            norm_power[ids,:,:] = out

        #norm_power = norm_power / (scale_dummy * scale_dummy) # Normalized wavelet power spectrum
        # Note: Liu et al 2007 JOAT suggest dividing by wavelet scale only - we emphasize small scales more.

        return wav_coeffs_pure, norm_power
