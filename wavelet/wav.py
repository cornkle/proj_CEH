"""

@author: C. Klein
"""
import numpy as np
from wavelet import twod as w2d
from scipy import ndimage
import pdb
import matplotlib.pyplot as plt
from wavelet import twod as w2d


class method(object):

    pass

    # def
    #
    # all_methods = {
    #
    #
    #
    #
    #
    # }

    def __init__(self, res, dist, start, nb, mother2d = w2d.Mexican_hat()):

        # 2D continuous wavelet analysis, this only supports dx == dy. By default, uses the Mexican hat wavelet
        #
        # dj: exponential distance between scales according to s0 * 2 ** (j * dj)
        # s0: start scale, approx 2*3*pixel scale (3 pixels necessary for signal detection)
        # j: number of scales

        s0 = 2 * start / mother2d.flambda()

        a = s0 * 2. ** (np.arange(0, nb + 1) * dist)  # The scales
        freqs = 1. / (mother2d.flambda() * a)  # As of Mallat 1999
        period = 1. / freqs
        scales = period/2.

        self.scale_dist = dist
        self.scale_start = s0
        self.scale_number = nb
        self.pixel = res
        self.scales = scales
        self.norm_scales = a

    def calc_coeff(self, data):

        p_spec = w2d.cwt2d(data, self.pixel, self.pixel, dj=self.scale_dist, s0=self.scale_start, J=self.scale_number)
