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

    def test_unique_of_tuple(self):

        x = np.array([5, 6, 3, 0, 1, 3, 0, 5])
        y = np.array([8, 6, 2, 0, 0, 2, 1, 8])

        xrand = np.arange(100)
        yrand = np.arange(100)
        np.random.shuffle(xrand)
        np.random.shuffle(yrand)

        uni = ua.unique_of_pair(x,y)
        unir = ua.unique_of_pair(xrand,yrand)

        assert np.unique(uni).size == x.size-2
        assert np.unique(unir).size == xrand.size
