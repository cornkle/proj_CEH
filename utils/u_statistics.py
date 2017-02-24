import numpy as np




def histo_frequency(data):
    weights = np.ones_like(data) / float(len(data))
    hist, h = np.histogram(data, bins=np.arange(0.1, 100 + 1, 1), weights=weights, range=(0.1, 100))
