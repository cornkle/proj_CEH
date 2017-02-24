import numpy as np
import scipy
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
import matplotlib.pyplot as plt

fname = np.empty((10,10))*0
putone = np.arange(0,10,1)
print(putone)
fname[putone, putone]=5


neighborhood_size = 5
threshold = 1

data_max = filters.maximum_filter(fname, (2, 2))

data_min = filters.minimum_filter(fname, (2, 2))

#data_max = filters.maximum_filter(data, neighborhood_size)
maxima = (fname == data_max)
# data_min = filters.minimum_filter(data, neighborhood_size)
diff = ((data_max - data_min) > threshold)

print((data_max - data_min))
maxima[diff == 0] = 0
#
labeled, num_objects = ndimage.label(maxima)
slices = ndimage.find_objects(labeled)
xx, yy = [], []
for dy, dx in slices:
     x_center = (dx.start + dx.stop - 1) / 2
     xx.append(x_center)
     y_center = (dy.start + dy.stop - 1) / 2
     yy.append(y_center)

#yy, xx = np.where((data_max == 1))
print(yy,xx)
plt.pcolormesh(fname)


plt.autoscale(False)
plt.plot(xx, yy, 'yo')
plt.show()