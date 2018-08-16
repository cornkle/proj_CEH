import numpy as np



def da_minmax(da):

    return (da - da.min()) / (da.max() - da.min())


def minmax(arr):

    return (arr - np.nanmin(arr)) / (np.nanmax(arr) - np.nanmin(arr))