#!/bin/env python

import matplotlib.pyplot as plt
import numpy as np

from wav_mixed.old_dani import twod_dani_old as w2d

x1 = np.arange(300)
y1 = np.arange(300)
x,y = np.meshgrid(x1,y1)
g = np.exp(-( (x-150)**2/100. + (y-150)**2/3600. ))

g = g - np.mean(g)
mother2d = w2d.Mexican_hat()
wavel2d, scales2d, freqs2d = w2d.cwt2d(g, 1, 1, dj=1./12, s0=4./mother2d.flambda(), J=104)

plt.figure()
plt.plot(g[150,:])
plt.plot(wavel2d[20,150,:])
plt.title('wavel2d')

wavel2d[np.real(wavel2d<=0)] = -0.01

power2d = (np.abs(wavel2d)) ** 2 # Normalized wavelet power spectrum
period2d = 1. / freqs2d
print('Old Scales', period2d/2)
scales2d.shape = (len(scales2d),1,1)
power2ds = power2d / (scales2d**2)

plt.figure()
plt.pcolormesh(g)

plt.figure()
plt.plot(g[150,:])
plt.plot(wavel2d[20,150,:])
plt.title('wavel2d after removal of neg val')

plt.figure()
plt.plot(power2d[20,150,:])
plt.title('power2d')

plt.figure()
plt.plot(power2ds[20,150,:])
plt.title('power2d scale normalised')