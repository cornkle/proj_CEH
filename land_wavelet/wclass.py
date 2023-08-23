from land_wavelet import constants, wav, wav1d
import numpy as np
from scipy.ndimage.measurements import label
from scipy import ndimage
import xarray as xr
import ipdb
import os

class landwav(object):

    def __init__(self, dataname):

        if dataname in constants.NAMES:
            dic = constants.NAMES[dataname]
        else:
            print('Dataset not found, ERROR')
            return

        self.name = dataname
        self.res = dic['dx']
        self.dist = dic['dist']
        self.nb = dic['nb']
        self.start = dic['start']


        wobj = wav.wavelet(self.res, self.dist, self.nb, start=self.start)
        wobj_1d = wav1d.wavelet_1d(self.res, self.dist, self.nb, start=self.start)
        self.scales = wobj.scales # actual scales at across-scale power maximum
        self.period = wobj.period # period scale better representative of coeff > 0 area.
        self.wobj = wobj
        self.wobj_1d = wobj_1d

        print('Initialised wavelet with scales: ', self.scales)

    def __repr__(self):
        return f"<landwav name attributes: name:{self.name} res:{self.res} dist:{self.dist} nb:{self.nb} start:{self.start}, period | " \
               f"image routine: image, original, lon, lat | wav routine: power, coeffs>"



    def read_img(self, image, lon, lat, vmin=None, vmax=None, nanfill=0):
        """
        Filters clouds of set area threshold and prepares image for wavelet analysis via adjusting background temperature
        and smoothing cloud edges.
        t: numpy array, cloud top temperature data
        lon: 1d numpy array, longitude or x
        lat: 1d numpy array, latitude or y
        edge_smoothing: optional cloud edge smoothing via gaussian filter - can help in case of excessive core
                        identification at cloud edges (default: False)
        dynamic_background: optional dynamical background temperature according to coldest pixel in image -
                            can help in case of excessive core identification at cloud edges (default: False)
        min_area: optional minimum area threshold for identified clouds. If false, minimum is defined by the minimum
                  core scale (default: False)

        :return: filtered cloud top temperatures with adjusted background temperature
        """

        londiff = lon[0:-1]-lon[1::]
        if not np.allclose(londiff, np.zeros_like(londiff)+londiff[0]):
            print('Mean res is', str(np.abs(np.mean(londiff))), ' I found that grid is not regular. If in doubt, please check.')

        if not isinstance(image, np.ndarray):
            print('ERROR: Input needs to be a numpy array, please check.')
            return


        t = image.copy()

        if vmin is not None:
            t[t <= vmin] = 0

        if vmax is not None:
            t[t >= vmin] = 0


        self.invalid = np.isnan(t)
        t[np.isnan(t)] = nanfill

        self.imean = np.mean(t)

        t = t - self.imean

        self.image = t
        self.original = image
        self.lon = lon
        self.lat = lat


    def applyWavelet(self, ge_thresh=None, fill=0.01, le_thresh=None, normed='none'):
        """
        Applies the wavelet functions and handles wavelet coefficient filtering.
        :param ge_thresh: greater-equal threshold for coefficient filtering.
        :param fill: fill value for filtering thresholds
        :param le_thresh: less-equal threshold for coefficient filtering.
        :return: Wavelet coefficient and wavelet power attributes of the wavelet object.
        """

        try:
            data = self.image.copy()
        except NameError:
            print('No image found to apply wavelet. Please read in an image first.')
            return

        print('Wavelet coeffs (none or stddev) and power (none, stddev or scale) normed by:', normed, 'Please note: Choose none if value reconstruction is intended.')

        #obj = wav.wavelet(self.res, self.dist, self.nb, start=self.start)

        coeffs, power = self.wobj.calc_coeffs(data, ge_thresh=ge_thresh, fill=fill, le_thresh=le_thresh, power_normed=normed)

        self.power = power
        self.coeffs = coeffs

        del data

        return coeffs, power, self.scales, self.period

    def applyWavelet_1d(self, ge_thresh=None, fill=0.01, le_thresh=None, normed='none', direction='x'):
        """
        Applies the wavelet functions and handles wavelet coefficient filtering.
        :param ge_thresh: greater-equal threshold for coefficient filtering.
        :param fill: fill value for filtering thresholds
        :param le_thresh: less-equal threshold for coefficient filtering.
        :return: Wavelet coefficient and wavelet power attributes of the wavelet object.
        """

        try:
            data = self.image.copy()
        except NameError:
            print('No image found to apply wavelet. Please read in an image first.')
            return

        print('Wavelet coeffs (none or stddev) and power (none, stddev or scale) normed by:', normed, 'Please note: Choose none if value reconstruction is intended.')

        #obj = wav.wavelet(self.res, self.dist, self.nb, start=self.start)

        coeffs, power = self.wobj_1d.calc_coeffs(data, ge_thresh=ge_thresh, fill=fill, le_thresh=le_thresh, normed=normed, direction=direction)

        self.power = power
        self.coeffs = coeffs

        del data

        return coeffs, power, self.scales, self.period


    def applyInverseWavelet(self, per_scale=False, anomaly_out=True):
        """
        Applies the wavelet functions and handles wavelet coefficient filtering.
        :param ge_thresh: greater-equal threshold for coefficient filtering.
        :param fill: fill value for filtering thresholds
        :param le_thresh: less-equal threshold for coefficient filtering.
        :return: Wavelet coefficient and wavelet power attributes of the wavelet object.
        """

        #obj = wav.wavelet(self.res, self.dist, self.nb, start=self.start)
        variable = self.wobj.calc_coeffs_inverse(self.coeffs, per_scale=per_scale)

        if anomaly_out==False:
            variable = variable + self.imean

        return variable, self.scales











