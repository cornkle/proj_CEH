import numpy as np
from numpy.testing import assert_array_equal
import unittest
from scipy import ndimage



class TestComposite(unittest.TestCase):
    def test_locmin(self):
        arr = np.ones((9, 12))
        arr[4, 3] = -2
        arr[4, 7] = -5

        filter = ndimage.minimum_filter(arr, 8, mode='constant', cval=np.min(arr) - 1)
        filter[filter==np.max(arr)] = -999
        maxoutt = (arr==filter)
        yy, xx = np.where((maxoutt == 1))

        assert (yy, xx) == (4, 7)

        filter = ndimage.minimum_filter(arr, 3, mode='constant', cval=np.min(arr) - 1)
        filter[filter == np.max(arr)] = -999
        maxoutt = (arr == filter)
        yy, xx = np.where((maxoutt == 1))

        assert_array_equal(np.array([4,4]),yy)
        assert_array_equal(np.array([3, 7]), xx)


