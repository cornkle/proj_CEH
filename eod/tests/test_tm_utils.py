from eod import tm_utils
import numpy as np
from numpy.testing import assert_array_equal, assert_allclose
import xarray as xr
import unittest

class TestTMUtils(unittest.TestCase):

    def test_minute_delta(self):

        tt = [59, 0, 15, 13, 9, 2, 31, 27, 43, 45]
        result = [1, 0, 0, 2, 6, -2, -1, 3, 2, 0]

        for min, res in zip(tt, result):
            assert tm_utils.minute_delta(min, 15) == res



