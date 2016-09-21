from eod import tm_utils, trmm
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

    def test_unique_of_tuple(self):

        x = np.array([5, 6, 3, 0, 1, 3, 0, 5])
        y = np.array([8, 6, 2, 0, 0, 2, 1, 8])

        xrand = np.arange(100)
        yrand = np.arange(100)
        np.random.shuffle(xrand)
        np.random.shuffle(yrand)

        uni = tm_utils.unique_of_pair(x,y)
        unir = tm_utils.unique_of_pair(xrand,yrand)

        assert np.unique(uni).size == x.size-2
        assert np.unique(unir).size == xrand.size

    def test_ll_to_MSG(self):

        zp = 1856

        #  1 x/y ca 0.025 lon/lat = ca. 3km

        test = [(0,0), (-0.01, 0), (-0.025, 0), (0,0.025), (0,-0.1)]
        result = [(zp, zp), (zp, zp), (zp+1, zp), (zp,zp+1), (zp,zp-4)]

        for t, r in zip(test, result):
            dir = tm_utils.ll_toMSG(t[0], t[1])
            assert (dir['x'], dir['y']) == r

    def test_ll_to_MSG_rev(self):

        zp = 1856

        #  1 x/y ca 0.025 lon/lat = ca. 3km

        test = [(0,0), (-0.01, 0), (-0.025, 0), (0,0.025), (0,-0.1)]
        result = [(zp, zp), (zp, zp), (zp+1, zp), (zp,zp+1), (zp,zp-4)]

        for t, r in zip(test, result):
            dir = tm_utils.ll_toMSG_rev(t[0], t[1])
            assert (dir['x'], dir['y']) == r

    def test_ll_to_MSG_TRMM(self):

        test_dir = '/users/global/cornkle/data/pythonWorkspace/proj_CEH/eod/tests/test_files/trmm'
        obj = trmm.ReadWA(test_dir)
        dat = obj.get_data(obj.fpaths[0], cut=[3,4])
        lon = dat['lon'].values
        lat = dat['lat'].values

        dir = tm_utils.ll_toMSG(lon, lat)
        assert np.unique(tm_utils.unique_of_pair(dir['x'],dir['y'])).size == lon.size

        #[item for item, count in Counter(a.flatten()).items() if count > 1]

       # In[112]: nnok
       # Out[112]: (array([1883, 1932]),)
       # In[113]: lon.flatten()[nnok]
       # Out[113]: array([15.71, 15.73])
       # In[114]: lat.flatten()[nnok]
       # Out[114]: array([2.36, 2.34])

        #Out[114]: array([2.36, 2.34])
        #In[115]: tm_utils.ll_toMSG(15.71, 2.36)
        #Out[115]: {'x': 1285.0, 'y': 1942.0}
        #In[116]: tm_utils.ll_toMSG(15.73, 2.34)
        #Out[116]: {'x': 1285.0, 'y': 1942.0}

        check why these are the same!! and how far they are apart..(3km distance)














