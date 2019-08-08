from __future__ import division

import numpy as np
import kpywavelet as kpy
import matplotlib.pyplot as plt
import xarray as xr
import ipdb

da = xr.open_dataset('/home/ck/DIR/mymachine/cores_fromMeteosat/cores/coresPower_MSG_-40_700km2_-50points_dominant_2012_04.nc')
datat = da['tir'].values[10,:,:]/100.

dt=5

out = np.zeros((10,datat.shape[0], datat.shape[1]))

for row, data in enumerate(datat):

    try:
        if data.shape[1] == 2: # For the last two data sets which contain time and data
            data = np.asarray(zip(*data)[1])
    except:
        pass

    std2 = data.std() ** 2
    data = (data - data.mean()) / data.std() # Calculating anomaly and normalizing
    time = np.arange(0, data.size)  # Time array in time units of your choice
    #alpha, _, _ = kpy.wavelet.ar1(data) # Lag-1 autocorrelation for white noise
    mother = kpy.wavelet.Mexican_hat() # Morlet mother wavelet with wavenumber=6

    wave, scales, freqs, coi, dj, s0, J = kpy.wavelet.cwt(data, dt, dj=0.5, s0=15, J=9, wavelet=mother)
    print('Scales', scales, len(scales))

    power = (np.abs(wave)) ** 2 # Normalized wavelet power spectrum

    period = 1. / freqs


    # Calculates the global wavelet spectrum and determines its significance level.
    glbl_power = std2  * power.mean(axis=1)
    dof = data.size - scales # Correction for padding at edges


    out[:,row,:] = power

plt.figure()
plt.contourf(datat)
plt.contour(out[5,:,:])

plt.figure()
plt.pcolormesh(out[3,:,:])


plt.figure()
plt.pcolormesh(out[:,50,:])