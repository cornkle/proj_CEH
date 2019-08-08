"""
Bi-dimensional continuous wavelet transform module for Python. Includes a 
collection of routines for wavelet transform and statistical analysis via
FFT algorithm. This module references to the numpy, scipy and pylab Python
packages.
DISCLAIMER
    This software may be used, copied, or redistributed as long as it
    is not sold and this copyright notice is reproduced on each copy
    made. This routine is provided as is without any express or implied
    warranties whatsoever.
AUTHOR
    Sebastian Krieger
    email: sebastian@nublia.com
REVISION
    1 (2011-04-30 19:48 -3000)
    2 (2016-02-04: DB introduced scale/period relationship as in 1D case)
REFERENCES
    [1] Wang, Ning and Lu, Chungu (2010). Two-dimensional continuous
        wavelet analysis and its application to meteorological data
"""

__version__ = '$Revision: 1 $'
# $Source$

from numpy import (arange, ceil, concatenate, conjugate, cos, exp, floor, 
                   isnan, log, log2, meshgrid, ones, pi, prod, real, sqrt,
                   zeros, polyval)
from numpy.fft import ifft2, fftfreq, fft2
import pdb
from pylab import find


class Mexican_hat():
    """Implements the Mexican hat wavelet class."""

    name = 'Mexican hat'
    
    def __init__(self):
        # Reconstruction factor $C_{\psi, \delta}$
        self.cpsi = 1. # pi

    def psi_ft(self, k, l):
        """
        Fourier transform of the Mexican hat wavelet as in Wang and
        Lu (2010), equation [15].
 
        """
        K, L = meshgrid(k, l)
        return (K ** 2. + L ** 2.) * exp(-0.5 * (K ** 2. + L ** 2.))

    def psi(self, x, y):
        """Mexican hat wavelet as in Wang and Lu (2010), equation [14]."""
        X, Y = meshgrid(x, y)
        return (2. - (X ** 2. + Y ** 2.)) * exp(-0.5 * (X ** 2. + Y ** 2.))
        
    def flambda(self):
        """Fourier wavelength as of Torrence and Compo (1998)."""
        return (2 * pi / sqrt(2.5))


def cwt2d(f, dx, dy, dj=1./12, s0=-1, J=-1, wavelet=Mexican_hat()):
    """
    Bi-dimensional continuous wavelet transform of the signal at 
    specified scale a.
    PARAMETERS
        f (array like):
            Input signal array.
        dx, dy (float):
            Sample spacing for each dimension.
        dj (float, optional) :
            Spacing between discrete scales. Default value is 0.25.
            Smaller values will result in better scale resolution, but
            slower calculation and plot.
        s0 (float, optional) :
            Smallest scale of the wavelet. Default value is 2*dt.
        J (float, optional) :
            Number of scales less one. Scales range from s0 up to
            s0 * 2**(J * dj), which gives a total of (J + 1) scales.
            Default is J = (log2(N*dt/so))/dj.
        a (array like, optional):
            Scale parameter array.
        wavelet (class, optional) :
            Mother wavelet class. Default is Mexican_hat()
    RETURNS
        Wf (array like) :
            2D wavelet transform according to the selected mother wavelet.
            Has (J+1) x N x M dimensions.
        a (array like) :
            Vector of scale indices given by a = s0 * 2**(j * dj),
            j={0, 1, ..., J}.
        freqs (array like) :
            Vector of Fourier frequencies (in 1 / time units) that
            corresponds to the wavelet scales.
    EXAMPLE
        wave, scales, freqs = twod.cwt2d(var, 9., 9., 1./12, -1, -1)
    """
    #**************DB NOTE: new scales work ONLY for dx = dy!!!!!********************
    # Determines the shape of the arrays and the discrete scales.

    n0, m0 = f.shape
    if s0 == -1: s0 = 2 * max(dx,dy) / wavelet.flambda()  # Smallest resolvable scale
    if J == -1: J = int(log2(max(n0,m0) * max(dx,dy) / s0) / dj)  # Number of scales
    N, M = 2 ** int(ceil(log2(n0))), 2 ** int(ceil(log2(m0)))   # Next higher power of 2
    
    a = s0 * 2. ** (arange(0, J+1) * dj)         # The scales
    A = len(a)
    # Calculates the zonal and meridional wave numbers.
    l, k = 2 * pi * fftfreq(N, dy), 2 * pi * fftfreq(M, dx)
    # Calculates the Fourier transform of the input signal.
    f_ft = fft2(f, s=(N, M))
    # Creates empty wavelet transform array and fills it for every discrete
    # scale using the convolution theorem.
    Wf = zeros((A, N, M), 'complex')
    for i, an in enumerate(a):
        psi_ft_bar = an * wavelet.psi_ft(an * k, an * l)
        Wf[i, :, :] = ifft2(f_ft * psi_ft_bar, s=(N, M))

    return Wf[:, :n0, :m0]


def icwt2d(W, a, dx=0.25, dy=0.25, da=0.25, wavelet=Mexican_hat()):
    """
    Inverse bi-dimensional continuous wavelet transform as in Wang and
    Lu (2010), equation [5].
    PARAMETERS
        W (array like):
            Wavelet transform, the result of the cwt2d function.
        a (array like, optional):
            Scale parameter array.
        w (class, optional) :
            Mother wavelet class. Default is Mexican_hat()
    RETURNS
        iW (array like) :
            Inverse wavelet transform.
    EXAMPLE
    """
    m0, l0, k0 = W.shape
    if m0 != a.size:
        raise Warning('Scale parameter array shape does not match wavelet' \
                       ' transform array shape.')
    # Calculates the zonal and meridional wave numters.
    L, K = 2 ** int(ceil(log2(l0))), 2 ** int(ceil(log2(k0)))
    # Calculates the zonal and meridional wave numbers.
    l, k = fftfreq(L, dy), fftfreq(K, dx)
    # Creates empty inverse wavelet transform array and fills it for every 
    # discrete scale using the convolution theorem.
    iW = zeros((m0, L, K), 'complex')
    for i, an in enumerate(a):
        psi_ft_bar = an * wavelet.psi_ft(an * k, an * l)
        W_ft = fft2(W[i, :, :], s=(L, K))
        iW[i, :, :] = ifft2(W_ft * psi_ft_bar, s=(L, K)) * da / an ** 2.

    return iW[:, :l0, :k0].real.sum(axis=0) / wavelet.cpsi