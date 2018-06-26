#!/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from wavelet import wav
from wavelet import twod as w2d


x1 = np.arange(300)
y1 = np.arange(300)
x,y = np.meshgrid(x1,y1)
g = np.exp(-( (x-150)**2/100. + (y-150)**2/3600. ))

g = g - np.mean(g)
mother2d = w2d.Mexican_hat()
obj = wav.wavelet(1, 1./12, 104)
# TIR
# wg, powerg = obj.calc_coeffs(g)
#
# plt.figure()
# plt.plot(wg[20,150,:])
# plt.title('wavel2d')
#
# plt.figure()
# plt.plot(powerg[20,150,:])
# plt.title('wavel2d power, all coeffs stay')

wg, powerg = obj.calc_coeffs(g, le_thresh=0, fill=-0.01)

plt.figure()
plt.plot(wg[20,150,:])
plt.title('wavel2d')

plt.figure()
plt.plot(powerg[20,150,:])
plt.title('power with thresh -0.01')



# wg, powerg = obj.calc_coeffs(g, le_thresh=0, fill=0)
#
# plt.figure()
# plt.plot(wg[20,150,:])
# plt.title('wavel2d')
#
# plt.figure()
# plt.plot(powerg[20,150,:])
# plt.title('wavel2d with fill 0')