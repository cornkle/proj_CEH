from utils import u_arrays as ua
import numpy as np
from numpy.testing import assert_array_equal
import unittest


class TestUtils(unittest.TestCase):
    def test_distance(self):
        data = [(2, 4, 2, 4), (2, 4, 2, 6),
                (np.array([2, 2, 1]), np.array([4, 4, 1]), np.array([2, 2, 1]), np.array([4, 6, 4]))]
        res = [np.array(0), np.array(2), np.array([0, 2, 3])]

        for d, r in zip(data, res):
            dist = ua.distance(d[0], d[1], d[2], d[3])
            assert_array_equal(dist, r)

    def test_closest_point(self):
        data = np.array([(1, 2), (3, 4), (5, 1), (7, 1), (0, 0)])
        #data2 = np.array([(1,3,5,7,0), (2,4,1,1,0)])
        res = 3
        point = np.array((7, 2))

        assert ua.closest_point(point, data) == res
