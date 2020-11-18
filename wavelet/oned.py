"""
Continuous wavelet transform module for Python. Includes a collection
of routines for wavelet transform and statistical analysis via FFT
algorithm. This module references to the numpy, scipy and pylab Python
packages.
"""
from __future__ import division

import numpy as np
import numpy.fft as fft
from numpy.lib.polynomial import polyval
from scipy.special import gamma
from scipy.special.orthogonal import hermitenorm


class DOG:
    """
    Implements the derivative of a Guassian wavelet class.
    Note that the input parameter f is the angular frequency and that
    for m=2 the DOG becomes the Mexican hat wavelet and that
    the default order for this wavelet is m=6.
    #TODO: Implenment arbitarty order
    """
    def __init__(self, m=6):
        self._set_m(m)
        self.name = 'DOG'

    def psi_ft(self, f):
        """Fourier transform of the DOG wavelet."""
        return (- 1j ** self.m / np.sqrt(gamma(self.m + 0.5)) * f ** self.m *
                np.exp(- 0.5 * f ** 2))

    def psi(self, t):
        """DOG wavelet as described in Torrence and Compo (1998)
        The derivative of a Gaussian of order n can be determined using
        the probabilistic Hermite polynomials. They are explicitly
        written as:
            Hn(x) = 2 ** (-n / s) * n! * sum ((-1) ** m) * (2 ** 0.5 *
                x) ** (n - 2 * m) / (m! * (n - 2*m)!)
        or in the recursive form:
            Hn(x) = x * Hn(x) - nHn-1(x)
        Source: http://www.ask.com/wiki/Hermite_polynomials
        """
        p = hermitenorm(self.m)
        return ((-1) ** (self.m + 1) * polyval(p, t) * np.exp(-t ** 2 / 2) /
                np.sqrt(gamma(self.m + 0.5)))

    def flambda(self):
        """Fourier wavelength as of Torrence and Compo (1998)."""
        return (2 * np.pi / np.sqrt(self.m + 0.5))

    def coi(self):
        """e-Folding Time as of Torrence and Compo (1998)."""
        return 1. / np.sqrt(2.)

    def sup(self):
        """Wavelet support defined by the e-Folding time."""
        return 1. / self.coi

    def _set_m(self, m):
        # Sets the m derivative of a Gaussian, the degrees of freedom and the
        # empirically derived factors for the wavelet bases C_{\delta}, \gamma,
        # \delta j_0 (Torrence and Compo, 1998, Table 2)
        self.m = m               # m-derivative
        self.dofmin = 1          # Minimum degrees of freedom
        if self.m == 2:
            self.cdelta = 3.541  # Reconstruction factor
            self.gamma = 1.43    # Decorrelation factor for time averaging
            self.deltaj0 = 1.40  # Factor for scale averaging
        elif self.m == 6:
            self.cdelta = 1.966
            self.gamma = 1.37
            self.deltaj0 = 0.97
        else:
            self.cdelta = -1
            self.gamma = -1
            self.deltaj0 = -1


class Mexican_hat(DOG):
    """
    Implements the Mexican hat wavelet class.
    This class inherits the DOG class using m=2.
    """
    def __init__(self):
        self._set_m(2)
        self.name = 'Mexican Hat'

def rect(x, normalize=False) :
    if type(x) in [int, float]:
        shape = [x, ]
    elif type(x) in [list, dict]:
        shape = x
    elif type(x) in [np.ndarray, np.ma.core.MaskedArray]:
        shape = x.shape
    X = np.zeros(shape)
    X[0] = X[-1] = 0.5
    X[1:-1] = 1

    if normalize:
        X /= X.sum()

    return X


def fftconv(x, y):
    """ Convolution of x and y using the FFT convolution theorem. """
    N = len(x)
    n = int(2 ** np.ceil(np.log2(N))) + 1
    X, Y, x_y = fft(x, n), fft(y, n), []
    for i in range(n):
        x_y.append(X[i] * Y[i])

    # Returns the inverse Fourier transform with padding correction
    return fft.ifft(x_y)[4:N+4]


def cwt(signal, dt, dj=1./12, s0=-1, J=-1, wavelet=Mexican_hat()):
    """
    Continuous wavelet transform of the signal at specified scales.
    Parameters
    ----------
        signal : numpy.ndarray, list
            Input signal array
        dt : float
            Sample spacing.
        dj : float, optional
            Spacing between discrete scales. Default value is 0.25.
            Smaller values will result in better scale resolution, but
            slower calculation and plot.
        s0 : float, optional
            Smallest scale of the wavelet. Default value is 2*dt.
        J : float, optional
            Number of scales less one. Scales range from s0 up to
            s0 * 2**(J * dj), which gives a total of (J + 1) scales.
            Default is J = (log2(N*dt/so))/dj.
        wavelet : instance of a wavelet class, optional
            Mother wavelet class. Default is Morlet wavelet.
    Returns
    -------
        W  : numpy.ndarray
            Wavelet transform according to the selected mother wavelet.
            Has (J+1) x N dimensions.
        sj : numpy.ndarray
            Vector of scale indices given by sj = s0 * 2**(j * dj),
            j={0, 1, ..., J}.
        freqs : array like
            Vector of Fourier frequencies (in 1 / time units) that
            corresponds to the wavelet scales.
        coi : numpy.ndarray
            Returns the cone of influence, which is a vector of N
            points containing the maximum Fourier period of useful
            information at that particular time. Periods greater than
            those are subject to edge effects.
        fft : numpy.ndarray
            Normalized fast Fourier transform of the input signal.
        fft_freqs : numpy.ndarray
            Fourier frequencies (in 1/time units) for the calculated
            FFT spectrum.
    Example
    -------
        mother = wavelet.Morlet(6.)
        wave, scales, freqs, coi, fft, fftfreqs = wavelet.cwt(var,
            0.25, 0.25, 0.5, 28, mother)
    """
    n0 = len(signal)                              # Original signal length.
    if s0 == -1: s0 = 2 * dt / wavelet.flambda()  # Smallest resolvable scale:
    if J == -1: J = int(np.log2(n0 * dt / s0) / dj)  # Number of scales
    N = 2 ** (int(np.log2(n0)) + 1)                  # Next higher power of 2.
    signal_ft = fft.fft(signal, N)                    # Signal Fourier transform
    ftfreqs = 2 * np.pi * fft.fftfreq(N, dt)             # Fourier angular frequencies

    sj = s0 * 2 ** (np.arange(0, J+1) * dj)         # The scales
    freqs = 1 / (wavelet.flambda() * sj)         # As of Mallat 1999

    # Creates an empty wavlet transform matrix and fills it for every discrete
    # scale using the convolution theorem.
    W = np.zeros((len(sj), N), 'complex')
    for n, s in enumerate(sj):
        psi_ft_bar = ((s * ftfreqs[1] * N) ** .5 *
            np.conjugate(wavelet.psi_ft(s * ftfreqs)))
        W[n, :] = fft.ifft(signal_ft * psi_ft_bar, N)

    # Checks for NaN in transform results and removes them from the scales,
    # frequencies and wavelet transform.
    sel = np.logical_not(np.isnan(W).all(axis=1))
    sj = sj[sel]
    freqs = freqs[sel]
    W = W[sel, :]

    # Determines the cone-of-influence. Note that it is returned as a function
    # of time in Fourier periods. Uses triangualr Bartlett window with non-zero
    # end-points.
    coi = (n0 / 2. - abs(np.arange(0, n0) - (n0 - 1) / 2))
    coi = wavelet.flambda() * wavelet.coi() * dt * coi
    #
    return W[:, :n0]


def cwt1d(signal2d, dt, dj=1./12, s0=-1, J=-1, wavelet=Mexican_hat(), direction='x'):


    out = np.zeros((J+1, signal2d.shape[0], signal2d.shape[1]))

    if direction == 'x':
        for row in range(signal2d.shape[0]):

            coeffs = cwt(signal2d[row,:], dt, dj=dj, s0=s0, J=J, wavelet=wavelet)
            out[:, row, :] = coeffs

    if direction == 'y':
        for column in range(signal2d.shape[1]):
            coeffs = cwt(signal2d[:, column], dt, dj=dj, s0=s0, J=J, wavelet=wavelet)
            out[:, :, column] = coeffs

    return out